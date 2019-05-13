/*
 * wheels_acc_limits_constraint.cc
 *
 *  Created on: Apr 2, 2019
 *      Author: vivian
 */

#include <iostream>
#include <towr/constraints/wheels_acc_limits_constraint.h>
#include <towr/variables/variable_names.h>

namespace towr {

WheelsAccLimitsConstraint::WheelsAccLimitsConstraint (Vector3d acc_max,
        									  	  	  double T, double dt, const EE& ee,
													  const SplineHolderDrive& spline_holder)
	:TimeDiscretizationConstraint(T, dt, "wheels-acc-limits-" + std::to_string(ee))
{
  acc_max_ = acc_max;
  ee_      = ee;

  ee_wheels_motion_ = spline_holder.ee_wheels_motion_.at(ee_);

  SetRows(GetNumberOfNodes()*k3D);  // 3 constraints per node (ax, ay, az)
}

int
WheelsAccLimitsConstraint::GetRow (int node, int dim) const
{
  return node*k3D + dim;
}

void
WheelsAccLimitsConstraint::UpdateConstraintAtInstance (double t, int k, VectorXd& g) const
{
  g.middleRows(GetRow(k, X), k3D) = ee_wheels_motion_->GetPoint(t).a();
}

void
WheelsAccLimitsConstraint::UpdateBoundsAtInstance (double t, int k, VecBound& bounds) const
{
  for (int dim=0; dim<k3D; ++dim) {
	bounds.at(GetRow(k,dim)) = ifopt::Bounds(-acc_max_(dim), acc_max_(dim));
  }
}

void
WheelsAccLimitsConstraint::UpdateJacobianAtInstance(double t, int k, std::string var_set, Jacobian& jac) const
{
  if (var_set == id::EEWheelsMotionNodes(ee_)) {
	jac.middleRows(GetRow(k,X), k3D) = ee_wheels_motion_->GetJacobianWrtNodes(t,kAcc);
  }
}

} /* namespace towr */


