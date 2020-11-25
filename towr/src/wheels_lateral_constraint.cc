/*
 * wheels_lateral_constraint.cc
 *
 *  Created on: Feb 13, 2020
 *      Author: vivian
 */

#include <iostream>
#include <towr/constraints/wheels_lateral_constraint.h>
#include <towr/variables/variable_names.h>

namespace towr {

WheelsLateralConstraint::WheelsLateralConstraint (const HeightMap::Ptr& terrain,
												  const EE& ee,
												  const SplineHolder& spline_holder)
	: ConstraintSet(kSpecifyLater, "wheels-nhc-" + std::to_string(ee))
{
  ee_ = ee;
  base_linear_  = spline_holder.base_linear_;
  base_angular_ = EulerConverter(spline_holder.base_angular_);
  ee_spline_ = spline_holder.ee_motion_.at(ee_);

  poly_durations_ = ee_spline_->GetPolyDurations();
  T_ = GetVectorGlobalTimeAtNodes(poly_durations_);

  terrain_ = terrain;

  n_constraints_per_node_ = 2;  // lateral velocity

}

void
WheelsLateralConstraint::InitVariableDependedQuantities (const VariablesPtr& x)
{
  ee_motion_ = x->GetComponent<NodesVariablesPhaseBased>(id::EEMotionNodes(ee_));

  for (int id=0; id<ee_motion_->GetNodes().size(); ++id) {
	  if (ee_motion_->IsConstantNode(id))
	  //if (ee_motion_->AdjacentPolyIdsConstant(id))
		  pure_stance_node_ids_.push_back(id);
  }

//  std::cout << "WheelsLateralConstraint" << std::endl;
//  std::cout << "Number of Nodes " << ee_motion_->GetNodes().size() << std::endl;
//  std::cout << "T vector size " << T_.size() << std::endl;
//  std::cout << "ee" << ee_ << " stance nodes: ";
//  for (int id : pure_stance_node_ids_)
//	  std::cout << "(" << id << ", " << T_.at(id) << ") ";
//  std::cout << std::endl;

  int constraint_count = pure_stance_node_ids_.size()*n_constraints_per_node_;
  SetRows(constraint_count);
}

std::vector<double>
WheelsLateralConstraint::GetVectorGlobalTimeAtNodes(VecDurations poly)
{
	std::vector<double> global_time;

	double sum_time = 0.0;
	for (auto t : poly) {
		global_time.push_back(sum_time);
		sum_time += t;
	}
	global_time.push_back(sum_time);

	return global_time;
}

WheelsLateralConstraint::VectorXd
WheelsLateralConstraint::GetValues () const
{
  VectorXd g(GetRows());

  int row = 0;
  auto nodes = ee_motion_->GetNodes();
  for (int id : pure_stance_node_ids_) {
	Vector3d ee_pos_w = nodes.at(id).p();
	Vector3d ee_vel_w = nodes.at(id).v();
	Vector3d ee_acc_w = ee_spline_->GetPoint(T_.at(id)).a();

	double t = T_.at(id);

	Eigen::Matrix3d w_R_b = base_angular_.GetRotationMatrixBaseToWorld(t);
	Eigen::Matrix3d b_R_w = w_R_b.transpose();

	// rotation matrix from contact frame to world frame
	Eigen::Matrix3d c_R_w;
	Vector3d n = terrain_->GetNormalizedBasis(HeightMap::Normal, ee_pos_w.x(), ee_pos_w.y());
	Vector3d heading = w_R_b * Vector3d::UnitX();
	Vector3d ey = n.cross(heading);
	Vector3d ex = ey.cross(n);
	Vector3d ez = n;
	c_R_w << ex, ey, ez;

	Vector3d ee_vel_c = c_R_w * ee_vel_w;
	Vector3d ee_acc_c = c_R_w * ee_acc_w;

    g(row++) = ee_vel_c(Y);
    g(row++) = ee_acc_c(Y);
  }

  return g;
}

WheelsLateralConstraint::VecBound
WheelsLateralConstraint::GetBounds () const
{
  return VecBound(GetRows(), ifopt::BoundZero);
}

WheelsLateralConstraint::Jacobian
WheelsLateralConstraint::SkewSymmetricMatrix (const Vector3d in) const
{
  Jacobian out(3,3);

  out.coeffRef(0,1) = -in(2); out.coeffRef(0,2) =  in(1);
  out.coeffRef(1,0) =  in(2); out.coeffRef(1,2) = -in(0);
  out.coeffRef(2,0) = -in(1); out.coeffRef(2,1) =  in(0);

  return out;
}

Eigen::Matrix3d
WheelsLateralConstraint::GetJacobianTerrainBasis(HeightMap::Direction basis, double x, double y) const
{
  Eigen::Matrix3d jac = Eigen::Matrix3d::Ones() * 1e-10;

  jac.col(X_) += terrain_->GetDerivativeOfNormalizedBasisWrt(basis, X_, x, y);
  jac.col(Y_) += terrain_->GetDerivativeOfNormalizedBasisWrt(basis, Y_, x, y);

  return jac;
}

void
WheelsLateralConstraint::FillJacobianBlock(std::string var_set, Jacobian& jac) const
{

  if (var_set == id::EEMotionNodes(ee_)) {
	int row = 0;
	auto nodes = ee_motion_->GetNodes();
	for (int id : pure_stance_node_ids_) {
		Vector3d ee_pos_w = nodes.at(id).p();
		Vector3d ee_vel_w = nodes.at(id).v();

		double t = T_.at(id);

		Vector3d ee_acc_w = ee_spline_->GetPoint(t).a();

		Vector3d n = terrain_->GetNormalizedBasis(HeightMap::Normal, ee_pos_w.x(), ee_pos_w.y());
		Eigen::Matrix3d jac_n = GetJacobianTerrainBasis(HeightMap::Normal, ee_pos_w.x(), ee_pos_w.y());

		Eigen::Matrix3d w_R_b = base_angular_.GetRotationMatrixBaseToWorld(t);
		Vector3d heading = w_R_b * Vector3d::UnitX();
		Vector3d ey = n.cross(heading);
		Vector3d ex = ey.cross(n);
		Vector3d ez = n;
		Eigen::Matrix3d Dey = - (SkewSymmetricMatrix(heading) * jac_n);
		Eigen::Matrix3d Dex = - (SkewSymmetricMatrix(n) * Dey) + (SkewSymmetricMatrix(ey) * jac_n);
		Eigen::Matrix3d Dez = jac_n;

		Jacobian jac1 = ee_vel_w.x() * (Dex * ee_spline_->GetJacobianWrtNodes(t, kPos)).sparseView().row(Y);
		Jacobian jac2 = ee_spline_->GetJacobianWrtNodes(t, kVel).row(X) * ex.y();
		Jacobian jac3 = ee_vel_w.y() * (Dey * ee_spline_->GetJacobianWrtNodes(t, kPos)).sparseView().row(Y);
		Jacobian jac4 = ee_spline_->GetJacobianWrtNodes(t, kVel).row(Y) * ey.y();
		Jacobian jac5 = ee_vel_w.z() * (Dez * ee_spline_->GetJacobianWrtNodes(t, kPos)).sparseView().row(Y);
		Jacobian jac6 = ee_spline_->GetJacobianWrtNodes(t, kVel).row(Z) * ez.y();
		jac.row(row++) = (jac1 + jac2 + jac3 + jac4 + jac5 + jac6);


		jac1 = ee_acc_w.x() * (Dex * ee_spline_->GetJacobianWrtNodes(t, kPos)).sparseView().row(Y);
		jac2 = ee_spline_->GetJacobianWrtNodes(t, kAcc).row(X) * ex.y();
		jac3 = ee_acc_w.y() * (Dey * ee_spline_->GetJacobianWrtNodes(t, kPos)).sparseView().row(Y);
		jac4 = ee_spline_->GetJacobianWrtNodes(t, kAcc).row(Y) * ey.y();
		jac5 = ee_acc_w.z() * (Dez * ee_spline_->GetJacobianWrtNodes(t, kPos)).sparseView().row(Y);
		jac6 = ee_spline_->GetJacobianWrtNodes(t, kAcc).row(Z) * ez.y();
		jac.row(row++) = (jac1 + jac2 + jac3 + jac4 + jac5 + jac6);
	}
//	  jac.coeffRef(row++,col) = dR.row(Y) * ee_acc_w;
  }

  if (var_set == id::base_ang_nodes) {
	int row = 0;
	auto nodes = ee_motion_->GetNodes();
	for (int id : pure_stance_node_ids_) {
		Vector3d ee_pos_w = nodes.at(id).p();
		Vector3d ee_vel_w = nodes.at(id).v();

		double t = T_.at(id);

		Vector3d ee_acc_w = ee_spline_->GetPoint(t).a();

		Vector3d n = terrain_->GetNormalizedBasis(HeightMap::Normal, ee_pos_w.x(), ee_pos_w.y());

		Eigen::Matrix3d w_R_b = base_angular_.GetRotationMatrixBaseToWorld(t);
		Jacobian Dheading = base_angular_.DerivOfRotVecMult(t,Vector3d::UnitX(),false);
		Jacobian Dey =   SkewSymmetricMatrix(n) * Dheading;
		Jacobian Dex = - SkewSymmetricMatrix(n) * Dey;

		jac.row(row++) = ee_vel_w.x() * Dex.row(Y) + ee_vel_w.y() * Dey.row(Y);
		jac.row(row++) = ee_acc_w.x() * Dex.row(Y) + ee_acc_w.y() * Dey.row(Y);
	}
//    jac.row(row++) = base_angular_.DerivOfRotVecMult(t,ee_acc_w,true).row(Y);
  }
}

} /* namespace towr */



