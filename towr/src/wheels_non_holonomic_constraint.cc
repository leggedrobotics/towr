/*
 * wheels_non_holonomic_constraint.cc
 *
 *  Created on: Apr 3, 2019
 *      Author: vivian
 */
#include <towr/constraints/wheels_non_holonomic_constraint.h>

namespace towr {

WheelsNonHolonomicConstraint::WheelsNonHolonomicConstraint (const HeightMap::Ptr& terrain,
															double T, double dt, const EE& ee,
															const SplineHolder& spline_holder)
	:TimeDiscretizationConstraint(T, dt, "wheels-nhc-" + std::to_string(ee))
{
  ee_ = ee;
  base_linear_  = spline_holder.base_linear_;
  base_angular_ = EulerConverter(spline_holder.base_angular_);
  ee_wheels_motion_ = spline_holder.ee_motion_.at(ee_);
  durations_ = spline_holder.phase_durations_.at(ee_);

  terrain_ = terrain;

  n_constraints_per_node_ = 1;  // lateral velocity and acceleration

  SetRows(GetNumberOfNodes()*n_constraints_per_node_);
}

void
WheelsNonHolonomicConstraint::UpdateConstraintAtInstance (double t, int k, VectorXd& g) const
{
  int row = k*n_constraints_per_node_;

  Vector3d ee_pos_w = ee_wheels_motion_->GetPoint(t).p();
  Vector3d ee_vel_w = ee_wheels_motion_->GetPoint(t).v();
  Vector3d ee_acc_w = ee_wheels_motion_->GetPoint(t).a();

  Eigen::Matrix3d w_R_b = base_angular_.GetRotationMatrixBaseToWorld(t);
  //Eigen::Matrix3d b_R_w = w_R_b.transpose();

  // rotation matrix from contact frame to world frame
  Eigen::Matrix3d c_R_w;
  Vector3d n = terrain_->GetNormalizedBasis(HeightMap::Normal, ee_pos_w.x(), ee_pos_w.y());
  Vector3d heading = w_R_b * Vector3d::UnitX();
  Vector3d ey = n.cross(heading);
  Vector3d ex = ey.cross(n);
  Vector3d ez = n;
  c_R_w << ex, ey, ez;

  Vector3d ee_vel_c = c_R_w * ee_vel_w;
  //Vector3d ee_acc_c = c_R_w * ee_acc_w;

//  std::cout << ee_ << ": vel_y = " << ee_vel_c(Y) << std::endl;

  g(row++) = ee_vel_c(Y);
//  g(row++) = ee_acc_c(Y);
}

void WheelsNonHolonomicConstraint::UpdateBoundsAtInstance(double t, int k, VecBound& bounds) const {
  int row = k * n_constraints_per_node_;
  double bound_if_not_contact = 1e20;
  if ( durations_->IsContactPhase(t)){
    bounds.at(row++) = ifopt::Bounds(-0.01, 0.01);

  }else{
    bounds.at(row++) = ifopt::Bounds(-bound_if_not_contact, bound_if_not_contact);

  }


  //  bounds.at(row++) = ifopt::BoundZero;
}

WheelsNonHolonomicConstraint::Jacobian
WheelsNonHolonomicConstraint::SkewSymmetricMatrix (const Vector3d in) const
{
  Jacobian out(3,3);

  out.coeffRef(0,1) = -in(2); out.coeffRef(0,2) =  in(1);
  out.coeffRef(1,0) =  in(2); out.coeffRef(1,2) = -in(0);
  out.coeffRef(2,0) = -in(1); out.coeffRef(2,1) =  in(0);

  return out;
}

Eigen::Matrix3d
WheelsNonHolonomicConstraint::GetJacobianTerrainBasis(HeightMap::Direction basis, double x, double y) const
{
  Eigen::Matrix3d jac = Eigen::Matrix3d::Ones() * 1e-10;

  jac.col(X_) += terrain_->GetDerivativeOfNormalizedBasisWrt(basis, X_, x, y);
  jac.col(Y_) += terrain_->GetDerivativeOfNormalizedBasisWrt(basis, Y_, x, y);

  return jac;
}

int
WheelsNonHolonomicConstraint::GetCol(int node, int dim) const
{
  int col = node*k6D + dim;
  int n_opt_variables = ee_wheels_motion_->GetNodeVariablesCount();

  if (col >= n_opt_variables)
    col = (node-1)*k6D + dim;

  return col;
}

void
WheelsNonHolonomicConstraint::UpdateJacobianAtInstance(double t, int k, std::string var_set, Jacobian& jac) const
{
  int row = k*n_constraints_per_node_;

  if (var_set == id::EEMotionNodes(ee_)) {
	Vector3d ee_pos_w = ee_wheels_motion_->GetPoint(t).p();
	Vector3d ee_vel_w = ee_wheels_motion_->GetPoint(t).v();
	Vector3d ee_acc_w = ee_wheels_motion_->GetPoint(t).a();

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

	Jacobian jac1 = ee_vel_w.x() * (Dex * ee_wheels_motion_->GetJacobianWrtNodes(t, kPos)).sparseView().row(Y);
	Jacobian jac2 = ee_wheels_motion_->GetJacobianWrtNodes(t, kVel).row(X) * ex.y();
	Jacobian jac3 = ee_vel_w.y() * (Dey * ee_wheels_motion_->GetJacobianWrtNodes(t, kPos)).sparseView().row(Y);
	Jacobian jac4 = ee_wheels_motion_->GetJacobianWrtNodes(t, kVel).row(Y) * ey.y();
	Jacobian jac5 = ee_vel_w.z() * (Dez * ee_wheels_motion_->GetJacobianWrtNodes(t, kPos)).sparseView().row(Y);
	Jacobian jac6 = ee_wheels_motion_->GetJacobianWrtNodes(t, kVel).row(Z) * ez.y();
	jac.row(row++) = (jac1 + jac2 + jac3 + jac4 + jac5 + jac6);
//	  jac.coeffRef(row++,col) = dR.row(Y) * ee_acc_w;
  }

  if (var_set == id::base_ang_nodes) {
	Vector3d ee_pos_w = ee_wheels_motion_->GetPoint(t).p();
	Vector3d ee_vel_w = ee_wheels_motion_->GetPoint(t).v();
	Vector3d ee_acc_w = ee_wheels_motion_->GetPoint(t).a();

	Vector3d n = terrain_->GetNormalizedBasis(HeightMap::Normal, ee_pos_w.x(), ee_pos_w.y());

	Eigen::Matrix3d w_R_b = base_angular_.GetRotationMatrixBaseToWorld(t);
	Jacobian Dheading = base_angular_.DerivOfRotVecMult(t,Vector3d::UnitX(),false);
	Jacobian Dey =   SkewSymmetricMatrix(n) * Dheading;
	Jacobian Dex = - SkewSymmetricMatrix(n) * Dey;

	jac.row(row++) = ee_vel_w.x() * Dex.row(Y) + ee_vel_w.y() * Dey.row(Y);
//    jac.row(row++) = base_angular_.DerivOfRotVecMult(t,ee_acc_w,true).row(Y);
  }
}

} /* namespace towr */



