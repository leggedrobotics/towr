/*
 * wheels_non_holonomic_constraint.cc
 *
 *  Created on: Apr 3, 2019
 *      Author: vivian
 */
#include <iostream>
#include <towr/constraints/wheels_non_holonomic_constraint.h>
#include <towr/variables/variable_names.h>

namespace towr {

WheelsNonHolonomicConstraint::WheelsNonHolonomicConstraint (double T, double dt, const EE& ee,
													  const SplineHolderDrive& spline_holder)
	:TimeDiscretizationConstraint(T, dt, "wheels-nhc-" + std::to_string(ee))
{
  ee_      = ee;
  base_linear_  = spline_holder.base_linear_;
  base_angular_ = EulerConverter(spline_holder.base_angular_);
  ee_wheels_motion_ = spline_holder.ee_wheels_motion_.at(ee_);

  SetRows(GetNumberOfNodes());  // one constraint per node
}

void
WheelsNonHolonomicConstraint::UpdateConstraintAtInstance (double t, int k, VectorXd& g) const
{
  // This constraint would only work if the initial position of the base is not rotated in relation to the fixed frame.
  // TODO: use the base frame to define the non-holonomic constraint
  g(k) = ee_wheels_motion_->GetPoint(t).v()(Y);
}

void
WheelsNonHolonomicConstraint::UpdateBoundsAtInstance (double t, int k, VecBound& bounds) const
{
  bounds.at(k) = ifopt::BoundZero;
}

void
WheelsNonHolonomicConstraint::UpdateJacobianAtInstance(double t, int k, std::string var_set, Jacobian& jac) const
{
  if (var_set == id::EEWheelsMotionNodes(ee_)) {
	jac.middleRows(k,1) = ee_wheels_motion_->GetJacobianWrtNodes(t,kVel).row(Y);
  }
}

} /* namespace towr */



