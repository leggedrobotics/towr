/*
 * stability_constraint.cc
 *
 *  Created on: Apr 12, 2019
 *      Author: vivian
 */

#include <iostream>

#include <towr/constraints/stability_constraint.h>

#include <towr/variables/variable_names.h>
#include <towr/variables/cartesian_dimensions.h>

namespace towr {

StabilityConstraint::StabilityConstraint (const DynamicModel::Ptr& m,
                                      	  double T, double dt,
										  const SplineHolderDrive& spline_holder)
    :TimeDiscretizationConstraint(T, dt, "stability")
{
  model_ = m;

  // link with up-to-date spline variables
  base_linear_  = spline_holder.base_linear_;
  base_angular_ = EulerConverter(spline_holder.base_angular_);
  ee_forces_    = spline_holder.ee_wheels_force_;
  ee_motion_    = spline_holder.ee_wheels_motion_;

  SetRows(GetNumberOfNodes());  // one constraint per node
}

void
StabilityConstraint::UpdateConstraintAtInstance(double t, int k, VectorXd& g) const
{
  UpdateModel(t);
  g(k) = model_->GetStabilityMeasure().second;// - 0.1745; // 10 deg minimum
}

void
StabilityConstraint::UpdateBoundsAtInstance(double t, int k, VecBound& bounds) const
{
  bounds.at(k) = ifopt::BoundGreaterZero;
}

void
StabilityConstraint::UpdateJacobianAtInstance(double t, int k, std::string var_set,
                                            Jacobian& jac) const
{
  UpdateModel(t);

  // sensitivity of dynamic constraint w.r.t base variables.
  if (var_set == id::base_lin_nodes) {
    Jacobian jac_base_lin_pos = base_linear_->GetJacobianWrtNodes(t,kPos);
    Jacobian jac_base_lin_acc = base_linear_->GetJacobianWrtNodes(t,kAcc);

    jac.row(k)  = (model_->GetStabilityDerivativeWrtBaseLinPos().transpose() * jac_base_lin_pos).sparseView();
    jac.row(k) += (model_->GetStabilityDerivativeWrtBaseLinAcc().transpose() * jac_base_lin_acc).sparseView();
  }

  // sensitivity of dynamic constraint w.r.t. endeffector variables
  for (int ee=0; ee<model_->GetEECount(); ++ee) {
    if (var_set == id::EEWheelsMotionNodes(ee)) {
      Jacobian jac_ee_pos = ee_motion_.at(ee)->GetJacobianWrtNodes(t,kPos);
      jac.row(k) += (model_->GetStabilityDerivativeWrtEEPos().col(ee).transpose() * jac_ee_pos).sparseView();
    }
  }

}

void
StabilityConstraint::UpdateModel (double t) const
{
  auto com = base_linear_->GetPoint(t);

  Eigen::Matrix3d w_R_b = base_angular_.GetRotationMatrixBaseToWorld(t);
  Eigen::Vector3d omega = base_angular_.GetAngularVelocityInWorld(t);
  Eigen::Vector3d omega_dot = base_angular_.GetAngularAccelerationInWorld(t);

  int n_ee = model_->GetEECount();
  std::vector<Eigen::Vector3d> ee_pos;
  std::vector<Eigen::Vector3d> ee_force;
  for (int ee=0; ee<n_ee; ++ee) {
    ee_force.push_back(ee_forces_.at(ee)->GetPoint(t).p());
    ee_pos.push_back(ee_motion_.at(ee)->GetPoint(t).p());
  }

  model_->SetCurrent(com.p(), com.a(), w_R_b, omega, omega_dot, ee_force, ee_pos);
}

} /* namespace towr */



