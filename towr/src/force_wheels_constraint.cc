/*
 * force_wheels_constraint.cc
 *
 *  Created on: Apr 1, 2019
 *      Author: vivian
 */

#include <iostream>
#include <towr/constraints/force_wheels_constraint.h>
#include <towr/variables/variable_names.h>
#include <towr/terrain/height_map.h>

namespace towr {

ForceWheelsConstraint::ForceWheelsConstraint (const HeightMap::Ptr& terrain, double ft_max,
        									  double T, double dt, const EE& ee,
											  const SplineHolderDrive& spline_holder)
	:TimeDiscretizationConstraint(T, dt, "force-wheels-" + std::to_string(ee))
{
  terrain_ = terrain;
  ft_max_  = ft_max;
  mu_      = terrain_->GetFrictionCoeff();
  ee_      = ee;

  ee_wheels_motion_ = spline_holder.ee_wheels_motion_.at(ee_);
  ee_wheels_force_  = spline_holder.ee_wheels_force_.at(ee_);

  n_constraints_per_node_ = 1 + 1 + 2*k2D; // positive normal force + torque limit + 4 friction pyramid constraint

  SetRows(GetNumberOfNodes()*n_constraints_per_node_);
}

void
ForceWheelsConstraint::UpdateConstraintAtInstance (double t, int k, VectorXd& g) const
{
  int row = k*n_constraints_per_node_;
  Vector3d f = ee_wheels_force_->GetPoint(t).p();
  Vector3d p = ee_wheels_motion_->GetPoint(t).p();

  // basis
  Vector3d n  = terrain_->GetNormalizedBasis(HeightMap::Normal, p.x(), p.y());
  Vector3d t1 = terrain_->GetNormalizedBasis(HeightMap::Tangent1, p.x(), p.y());
  Vector3d t2 = terrain_->GetNormalizedBasis(HeightMap::Tangent2, p.x(), p.y());

  // unilateral force
  g(row++) = f.transpose() * n; // >0 (unilateral forces)

  // torque limit
  g(row++) = f.transpose() * t1; // traction limits

  // frictional pyramid
  g(row++) = f.transpose() * (t1 - mu_*n); // t1 < mu*n
  g(row++) = f.transpose() * (t1 + mu_*n); // t1 > -mu*n
  g(row++) = f.transpose() * (t2 - mu_*n); // t2 < mu*n
  g(row++) = f.transpose() * (t2 + mu_*n); // t2 > -mu*n

//  g.middleRows(row,3) = n;

}

void
ForceWheelsConstraint::UpdateBoundsAtInstance (double t, int k, VecBound& bounds) const
{
  int row = k*n_constraints_per_node_;

  bounds.at(row++) = ifopt::BoundGreaterZero; // unilateral forces
  bounds.at(row++) = ifopt::Bounds(-ft_max_, ft_max_); // traction force
  bounds.at(row++) = ifopt::BoundSmallerZero; // f_t1 <  mu*n
  bounds.at(row++) = ifopt::BoundGreaterZero; // f_t1 > -mu*n
  bounds.at(row++) = ifopt::BoundSmallerZero; // f_t2 <  mu*n
  bounds.at(row++) = ifopt::BoundGreaterZero; // f_t2 > -mu*n

//  bounds.at(row++) = ifopt::NoBound;
//  bounds.at(row++) = ifopt::NoBound;
//  bounds.at(row++) = ifopt::NoBound;
}

Eigen::Matrix3d
ForceWheelsConstraint::GetJacobianTerrainBasis(HeightMap::Direction basis, double x, double y) const
{
  Eigen::Matrix3d jac = Eigen::Matrix3d::Ones() * 1e-10;

  jac.col(X_) += terrain_->GetDerivativeOfNormalizedBasisWrt(basis, X_, x, y);
  jac.col(Y_) += terrain_->GetDerivativeOfNormalizedBasisWrt(basis, Y_, x, y);

  return jac;
}

//// Using this approach gives a problem with the sparsity structure of the Jacobian when building TOWR in Degug mode.
//// FIXED: if any of the basis vector have zeros, the number of nonZero elements of the Jacobian is not compatible to the
//// one defined in the original sparsity structure, which is not allowed by IFOPT, that's why the 1e-10 added to the derivatives.
//void
//ForceWheelsConstraint::UpdateJacobianAtInstance(double t, int k, std::string var_set, Jacobian& jac) const
//{
//  if (var_set == id::EEWheelsMotionNodes(ee_)) {
//	Vector3d f = ee_wheels_force_->GetPoint(t).p();
//	Vector3d p = ee_wheels_motion_->GetPoint(t).p();
//
//	Jacobian jac_n  = (f.transpose()*GetJacobianTerrainBasis(HeightMap::Normal, p.x(), p.y())).sparseView();
//	Jacobian jac_t1 = (f.transpose()*GetJacobianTerrainBasis(HeightMap::Tangent1, p.x(), p.y())).sparseView();
//	Jacobian jac_t2 = (f.transpose()*GetJacobianTerrainBasis(HeightMap::Tangent2, p.x(), p.y())).sparseView();
//
//	int row = k*n_constraints_per_node_;
//
////	jac.middleRows(row, 3) = GetJacobianTerrainBasis(HeightMap::Normal, p.x(), p.y()).sparseView() *
////							 ee_wheels_motion_->GetJacobianWrtNodes(t, kPos);
//
//	// unilateral force
//	jac.row(row++) = jac_n * ee_wheels_motion_->GetJacobianWrtNodes(t, kPos);
//
//	// traction force
//	jac.row(row++) = jac_t1 * ee_wheels_motion_->GetJacobianWrtNodes(t, kPos);
//
//	// friction force tangent 1 derivative
//	jac.row(row++) = (jac_t1-mu_*jac_n) * ee_wheels_motion_->GetJacobianWrtNodes(t, kPos);
//	jac.row(row++) = (jac_t1+mu_*jac_n) * ee_wheels_motion_->GetJacobianWrtNodes(t, kPos);
//
//	// friction force tangent 2 derivative
//	jac.row(row++) = (jac_t2-mu_*jac_n) * ee_wheels_motion_->GetJacobianWrtNodes(t, kPos);
//	jac.row(row++) = (jac_t2+mu_*jac_n) * ee_wheels_motion_->GetJacobianWrtNodes(t, kPos);
//
//  }
//
//  if (var_set == id::EEWheelsForceNodes(ee_)) {
//	Vector3d p  = ee_wheels_motion_->GetPoint(t).p();
//
//	Vector3d n  = terrain_->GetNormalizedBasis(HeightMap::Normal, p.x(), p.y()) + Vector3d::Constant(1e-10);
//	Vector3d t1 = terrain_->GetNormalizedBasis(HeightMap::Tangent1, p.x(), p.y()) + Vector3d::Constant(1e-10);
//	Vector3d t2 = terrain_->GetNormalizedBasis(HeightMap::Tangent2, p.x(), p.y()) + Vector3d::Constant(1e-10);
//
//	Jacobian jac_n  = n.transpose().sparseView();
//	Jacobian jac_t1 = t1.transpose().sparseView();
//	Jacobian jac_t2 = t2.transpose().sparseView();
//
//	int row = k*n_constraints_per_node_;
//
//	// unilateral force
//	jac.row(row++) = jac_n * ee_wheels_force_->GetJacobianWrtNodes(t, kPos);
//
//	// traction force
//	jac.row(row++) = jac_t1 * ee_wheels_force_->GetJacobianWrtNodes(t, kPos);
//
//	// friction force tangent 1 derivative
//	jac.row(row++) = (jac_t1-mu_*jac_n) * ee_wheels_force_->GetJacobianWrtNodes(t, kPos);
//	jac.row(row++) = (jac_t1+mu_*jac_n) * ee_wheels_force_->GetJacobianWrtNodes(t, kPos);
//
//	// friction force tangent 2 derivative
//	jac.row(row++) = (jac_t2-mu_*jac_n) * ee_wheels_force_->GetJacobianWrtNodes(t, kPos);
//	jac.row(row++) = (jac_t2+mu_*jac_n) * ee_wheels_force_->GetJacobianWrtNodes(t, kPos);
//
//  }
//}


int
ForceWheelsConstraint::GetCol(int node, int dim) const
{
  int col = node*k6D + dim;
  int n_opt_variables = ee_wheels_motion_->GetNodeVariablesCount();

  if (col >= n_opt_variables)
    col = (node-1)*k6D + dim;

  return col;
}

void
ForceWheelsConstraint::UpdateJacobianAtInstance(double t, int k, std::string var_set, Jacobian& jac) const
{
  if (var_set == id::EEWheelsMotionNodes(ee_)) {
	Vector3d f = ee_wheels_force_->GetPoint(t).p();
	Vector3d p = ee_wheels_motion_->GetPoint(t).p();
	int row_reset = k*n_constraints_per_node_;
	int row = row_reset;
	for (auto dim : {X_,Y_}) {
	  Vector3d dn  = terrain_->GetDerivativeOfNormalizedBasisWrt(HeightMap::Normal, dim, p.x(), p.y());
	  Vector3d dt1 = terrain_->GetDerivativeOfNormalizedBasisWrt(HeightMap::Tangent1, dim, p.x(), p.y());
	  Vector3d dt2 = terrain_->GetDerivativeOfNormalizedBasisWrt(HeightMap::Tangent2, dim, p.x(), p.y());

	  int col = GetCol(k, dim);

	  // unilateral force
	  jac.coeffRef(row++, col) = f.transpose()*dn;

	  // traction force
	  jac.coeffRef(row++, col) = f.transpose()*dt1;

	  // friction force tangent 1 derivative
	  jac.coeffRef(row++, col) = f.transpose()*(dt1-mu_*dn);
	  jac.coeffRef(row++, col) = f.transpose()*(dt1+mu_*dn);

	  // friction force tangent 2 derivative
	  jac.coeffRef(row++, col) = f.transpose()*(dt2-mu_*dn);
	  jac.coeffRef(row++, col) = f.transpose()*(dt2+mu_*dn);

	  row = row_reset;
	}

//	jac.middleRows(row, 3) = GetJacobianTerrainBasis(HeightMap::Normal, p.x(), p.y()).sparseView() *
//							 ee_wheels_motion_->GetJacobianWrtNodes(t, kPos);
  }

  if (var_set == id::EEWheelsForceNodes(ee_)) {
	Vector3d p  = ee_wheels_motion_->GetPoint(t).p();
	Vector3d n  = terrain_->GetNormalizedBasis(HeightMap::Normal, p.x(), p.y());
	Vector3d t1 = terrain_->GetNormalizedBasis(HeightMap::Tangent1, p.x(), p.y());
	Vector3d t2 = terrain_->GetNormalizedBasis(HeightMap::Tangent2, p.x(), p.y());

	int row_reset = k*n_constraints_per_node_;
	int row = row_reset;
	for (auto dim : {X,Y,Z}) {
	  int col = GetCol(k, dim);
	  jac.coeffRef(row++, col) = n(dim);              // normal force
	  jac.coeffRef(row++, col) = t1(dim);             // traction force
	  jac.coeffRef(row++, col) = t1(dim)-mu_*n(dim);  // f_t1 <  mu*n
	  jac.coeffRef(row++, col) = t1(dim)+mu_*n(dim);  // f_t1 > -mu*n
	  jac.coeffRef(row++, col) = t2(dim)-mu_*n(dim);  // f_t2 <  mu*n
	  jac.coeffRef(row++, col) = t2(dim)+mu_*n(dim);  // f_t2 > -mu*n
	  row = row_reset;
	}
  }
}

} /* namespace towr */




