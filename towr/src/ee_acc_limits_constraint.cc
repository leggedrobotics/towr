/*
 * wheels_acc_limits_constraint.cc
 *
 *  Created on: Apr 2, 2019
 *      Author: vivian
 */

#include <iostream>
#include <towr/constraints/ee_acc_limits_constraint.h>
#include <towr/variables/variable_names.h>

namespace towr {

EEAccLimitsConstraint::EEAccLimitsConstraint (Vector3d acc_max,
        									  const EE& ee,
											  const SplineHolder& ee_motion_holder)
		: ConstraintSet(kSpecifyLater, "acc-limits-" + id::EEMotionNodes(ee))
{
  acc_max_ = acc_max;
  ee_      = ee;

  ee_motion_ = ee_motion_holder.ee_motion_.at(ee_);
  n_polys_ = ee_motion_->GetPolynomialCount();

  int n_nodes = n_polys_ + 1;
  T_ = ee_motion_->GetPolyDurations();

  SetRows(n_nodes*k3D);  // 3 constraints per node (ax, ay, az)
}

Eigen::VectorXd
EEAccLimitsConstraint::GetValues () const
{
  VectorXd g(GetRows());

  int row = 0;

  // first node
  VectorXd acc = ee_motion_->GetPoint(0, 0.0).a();
  for (auto dim : {X, Y, Z})
  	g(row++) = acc(dim);

  // rest of nodes
  for (int j=0; j<n_polys_; ++j) {
    VectorXd acc = ee_motion_->GetPoint(j, T_.at(j)).a();

    for (auto dim : {X, Y, Z})
    	g(row++) = acc(dim);
  }

  return g;
}

void
EEAccLimitsConstraint::FillJacobianBlock (std::string var_set, Jacobian& jac) const
{
  if (var_set == id::EEMotionNodes(ee_)) {
	int row = 0;

    // first node
	Jacobian acc = ee_motion_->GetJacobianWrtNodes(0, 0.0, kAcc);
	for (auto dim : {X, Y, Z})
	  jac.row(row++) = acc.row(dim);

	// rest of nodes
    for (int j=0; j<n_polys_; ++j) {
      Jacobian acc = ee_motion_->GetJacobianWrtNodes(j, T_.at(j), kAcc);

      for (auto dim : {X, Y, Z})
      	jac.row(row++) = acc.row(dim);
    }
  }
}

EEAccLimitsConstraint::VecBound
EEAccLimitsConstraint::GetBounds () const
{
  VecBound bounds(GetRows());

  int i = 0;

  for (int j=0; j<(n_polys_+1); ++j) {
    for (auto dim : {X, Y, Z})
    	bounds.at(i++) = ifopt::Bounds(-acc_max_(dim), acc_max_(dim));
  }

  return bounds;

}

//int
//EEAccLimitsConstraint::GetRow (int node, int dim) const
//{
//  return node*k3D + dim;
//}
//
//void
//EEAccLimitsConstraint::UpdateConstraintAtInstance (double t, int k, VectorXd& g) const
//{
//  g.middleRows(GetRow(k, X), k3D) = ee_motion_->GetPoint(t).a();
//}
//
//void
//EEAccLimitsConstraint::UpdateBoundsAtInstance (double t, int k, VecBound& bounds) const
//{
//  for (int dim=0; dim<k3D; ++dim) {
//	bounds.at(GetRow(k,dim)) = ifopt::Bounds(-acc_max_(dim), acc_max_(dim));
//  }
//}
//
//void
//EEAccLimitsConstraint::UpdateJacobianAtInstance(double t, int k, std::string var_set, Jacobian& jac) const
//{
//  if (var_set == id::EEMotionNodes(ee_)) {
//	jac.middleRows(GetRow(k,X), k3D) = ee_motion_->GetJacobianWrtNodes(t,kAcc);
//  }
//}

} /* namespace towr */


