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
  decision_ = spline_holder.ee_decision_.at(ee_);

  terrain_ = terrain;

  n_constraints_per_node_ = 1;  // lateral velocity and acceleration

    T_ = ee_wheels_motion_->GetPolyDurations();
  SetRows(GetNumberOfNodes()*n_constraints_per_node_);
}

void
WheelsNonHolonomicConstraint::UpdateConstraintAtInstance (double t, int k, VectorXd& g) const
{
  Vector3d ee_d = decision_->GetPoint(t).p();
    EulerConverter::MatrixSXd w_C_b = base_angular_.GetRotationMatrixBaseToWorld(t).transpose();
    Vector3d v_wrt_b = w_C_b*ee_wheels_motion_->GetPoint(t).v();
    int row = k*n_constraints_per_node_;
    g(row++) = ee_d.x()*v_wrt_b(1);
}

void WheelsNonHolonomicConstraint::UpdateBoundsAtInstance(double t, int k, VecBound& bounds) const {
  int row = k * n_constraints_per_node_;
    bounds.at(row++) = ifopt::Bounds(-0.01, 0.01);
}

void
WheelsNonHolonomicConstraint::UpdateJacobianAtInstance(double t, int k, std::string var_set, Jacobian& jac) const
{
  EulerConverter::MatrixSXd w_C_b = base_angular_.GetRotationMatrixBaseToWorld(t).transpose();
    Vector3d ee_d = decision_->GetPoint(t).p();

    Vector3d v_wrt_b = w_C_b * ee_wheels_motion_->GetPoint(t).v();
    int row = k*n_constraints_per_node_;
    Vector3d v_W = ee_wheels_motion_->GetPoint(t).v();

    if (var_set == id::base_ang_nodes) {
        jac.row(row++) = base_angular_.DerivOfRotVecMult(t, ee_d.x()*v_W, true).row(Y);
    }

    if (var_set == id::EEMotionNodes(ee_)) {
        jac.row(row++) = (ee_d.x()*w_C_b*ee_wheels_motion_->GetJacobianWrtNodes(t, kVel)).row(Y);
    }

    if (var_set == id::EEDecision(ee_)) {
      //this values are not optimized over generally, but the jacobian has been put here still if one were to change that in the future
      jac.row(row++) = (v_wrt_b(1)*decision_->GetJacobianWrtNodes(t, kPos)).row(X);
    }

    if (var_set == id::EESchedule(ee_)) {
      jac.row(row++) =
          (ee_d.x() * w_C_b *
           ee_wheels_motion_->GetJacobianOfVelWrtDurations(t))
              .row(Y) +
          (v_wrt_b(1) * decision_->GetJacobianOfVelWrtDurations(t)).row(X);
    }
}

} /* namespace towr */



