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
  ee_ = ee;
  base_linear_  = spline_holder.base_linear_;
  base_angular_ = EulerConverter(spline_holder.base_angular_);
  ee_wheels_motion_ = spline_holder.ee_wheels_motion_.at(ee_);

  n_constraints_per_node_ = 2;  // lateral velocity and acceleration

  SetRows(GetNumberOfNodes()*n_constraints_per_node_);
}

void
WheelsNonHolonomicConstraint::UpdateConstraintAtInstance (double t, int k, VectorXd& g) const
{
  // This constraint would only work if the initial position of the base is not rotated in relation to the fixed frame.
  // TODO: use the base frame to define the non-holonomic constraint

  int row = k*n_constraints_per_node_;

  Eigen::Vector3d ee_vel_w = ee_wheels_motion_->GetPoint(t).v();
  Eigen::Vector3d ee_acc_w = ee_wheels_motion_->GetPoint(t).a();

  Eigen::Matrix3d b_R_w = base_angular_.GetRotationMatrixBaseToWorld(t).transpose();

  Eigen::Vector3d ee_vel_b = b_R_w * ee_vel_w;
  Eigen::Vector3d ee_acc_b = b_R_w * ee_acc_w;

  g(row++) = ee_vel_b(Y);
  g(row++) = ee_acc_b(Y);
}

void
WheelsNonHolonomicConstraint::UpdateBoundsAtInstance (double t, int k, VecBound& bounds) const
{
  int row = k*n_constraints_per_node_;
  bounds.at(row++) = ifopt::BoundZero;
  bounds.at(row++) = ifopt::BoundZero;
}

void
WheelsNonHolonomicConstraint::UpdateJacobianAtInstance(double t, int k, std::string var_set, Jacobian& jac) const
{
  int row = k*n_constraints_per_node_;
  EulerConverter::MatrixSXd b_R_w = base_angular_.GetRotationMatrixBaseToWorld(t).transpose();

  if (var_set == id::EEWheelsMotionNodes(ee_)) {
	jac.row(row++) = (b_R_w * ee_wheels_motion_->GetJacobianWrtNodes(t,kVel)).row(Y);
	jac.row(row++) = (b_R_w * ee_wheels_motion_->GetJacobianWrtNodes(t,kAcc)).row(Y);
  }

  if (var_set == id::base_ang_nodes) {
	Eigen::Vector3d ee_vel_w = ee_wheels_motion_->GetPoint(t).v();
	Eigen::Vector3d ee_acc_w = ee_wheels_motion_->GetPoint(t).a();

    jac.row(row++) = base_angular_.DerivOfRotVecMult(t,ee_vel_w,true).row(Y);
    jac.row(row++) = base_angular_.DerivOfRotVecMult(t,ee_acc_w,true).row(Y);
  }
}

} /* namespace towr */



