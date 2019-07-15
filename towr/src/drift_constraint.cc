/*
 * drift_constraint.cc
 *
 *  Created on: Apr 25, 2019
 *      Author: dominic landolf
 */

#include <towr/constraints/drift_constraint.h>
#include <towr/variables/variable_names.h>
#include <towr/variables/spline_holder.h>
#include <towr/variables/cartesian_dimensions.h>
#include <Eigen/Dense>

using namespace Eigen;

namespace towr {

DriftConstraint::DriftConstraint (const KinematicModel::Ptr& model,
                                                  double T, double dt,
                                                  const EE& ee,
                                                  const SplineHolder& spline_holder)
    :TimeDiscretizationConstraint(T, dt, "drift_motion-" + std::to_string(ee))
{
  base_linear_  = spline_holder.base_linear_;
  base_angular_ = EulerConverter(spline_holder.base_angular_);
  base_angular_for_angles_ = spline_holder.base_angular_;

  ee_ = ee;
  ee_motion_    = spline_holder.ee_motion_.at(ee_);
  ee_force_ 	= spline_holder.ee_force_.at(ee_);

  n_constraints_per_node_ = 1;

  SetRows(GetNumberOfNodes()*n_constraints_per_node_);
}

int
DriftConstraint::GetRow (int node, int dim) const
{
  return node*n_constraints_per_node_ + dim;
}

//used for the Jacobians
DriftConstraint::Vector3d
DriftConstraint::NormDerivative (const Vector3d& v) const
{
  double p = pow(v.x(),2) + pow(v.y(),2) + pow(v.z(),2);

  if (p < 0.00000001)
	  return Vector3d(0,0,0);
  else
	  return Vector3d(v.x()/sqrt(p), v.y()/sqrt(p), v.z()/sqrt(p));
}

void
DriftConstraint::UpdateConstraintAtInstance (double t, int k, VectorXd& g) const
{
  EulerConverter::MatrixSXd b_C_w = base_angular_.GetRotationMatrixBaseToWorld(t);
  EulerConverter::MatrixSXd w_C_b = base_angular_.GetRotationMatrixBaseToWorld(t).transpose();

  //introduce a constant spinning velocity of the wheels during drift
  double des_v_x_wheel = params_.vx_wheel_during_drift_;
  double pitch_angle = base_angular_for_angles_->GetPoint(t).p()(Y);
  Vector3d v_x_wheel(-des_v_x_wheel, 0.0, 0.0); //in base frame
  Vector3d v_x_wheel_world_with_z = b_C_w*v_x_wheel;
  Vector3d v_x_wheel_world(v_x_wheel_world_with_z(X), v_x_wheel_world_with_z(Y), 0.0);

  Vector3d f_tang_world(ee_force_->GetPoint(t).p()(X), ee_force_->GetPoint(t).p()(Y), 0.0);
  Vector3d f_deriv_tang_world(ee_force_->GetPoint(t).v()(X), ee_force_->GetPoint(t).v()(Y), 0.0);

  //total velocity of the wheels
  Vector3d v_total = v_x_wheel_world + ee_motion_->GetPoint(t).v();

  double v_tot_abs_value = v_total.norm();

  double f_tang_abs_value = f_tang_world.norm();
  double f_deriv_tang_abs_value = f_deriv_tang_world.norm();

  int row = k*n_constraints_per_node_;

  //constrain the tangential force and the total velocity to point in opposite directions during drift.
  g(row++) = f_tang_world.transpose()*v_total + f_tang_abs_value*v_tot_abs_value;

}

void
DriftConstraint::UpdateBoundsAtInstance (double t, int k, VecBound& bounds) const
{
	//only have bounds during drift phase
	int row = k*n_constraints_per_node_;

	if (params_.just_drive_){
		bounds.at(row++) = ifopt::NoBound;
	}
	else {
		if (ee_ == 0 or ee_ == 1){
				bounds.at(row++) = ifopt::NoBound;
		}
		else if (ee_ == 2 or ee_ == 3){
			if ((t <= params_.phase_duration_drive_1) or (t > (params_.phase_duration_drive_1+params_.phase_duration_drift))){
				bounds.at(row++) = ifopt::NoBound;
			}
			//this is the drift phase
			else {
				bounds.at(row++) = ifopt::BoundZero;
			}
		}
	}
}

void
DriftConstraint::UpdateJacobianAtInstance (double t, int k,
                                                   std::string var_set,
                                                   Jacobian& jac) const
{
  EulerConverter::MatrixSXd w_C_b = base_angular_.GetRotationMatrixBaseToWorld(t).transpose();
  EulerConverter::MatrixSXd b_C_w = base_angular_.GetRotationMatrixBaseToWorld(t);

  int row = k*n_constraints_per_node_;

    Vector3d f_tang_world(ee_force_->GetPoint(t).p()(X), ee_force_->GetPoint(t).p()(Y), 0.0);
    f_tang_world = f_tang_world + Vector3d::Constant(1e-10);
    double pitch_angle = base_angular_for_angles_->GetPoint(t).p()(Y);

    double des_v_x_wheel = params_.vx_wheel_during_drift_;
    Vector3d v_x_wheel(-des_v_x_wheel, 0.0, 0.0); //in base frame
    Vector3d v_x_wheel_world_with_z = b_C_w*v_x_wheel;
    Vector3d v_x_wheel_world(v_x_wheel_world_with_z(X), v_x_wheel_world_with_z(Y), 0.0);
    Vector3d v_x_wheel_base_with_y = w_C_b*v_x_wheel_world;
    Vector3d v_x_wheel_base(v_x_wheel_base_with_y(X), 0.0, 0.0);

    Vector3d v_total = v_x_wheel_world + ee_motion_->GetPoint(t).v() + Vector3d::Constant(1e-10);
    double v_tot_abs_value =  v_total.norm();
    double f_tang_abs_value = f_tang_world.norm();

  Vector3d Dg = NormDerivative(v_total) + Vector3d::Constant(1e-10);
  Vector3d Dg_f = NormDerivative(f_tang_world) + Vector3d::Constant(1e-10);

  if (var_set == id::base_ang_nodes) {
	  Vector3d derivativeVec_base = (f_tang_world + f_tang_abs_value * Dg);
	  Jacobian derivativeJac_base = derivativeVec_base.transpose().sparseView();

	  jac.row(row++) = derivativeJac_base * base_angular_.DerivOfRotVecMult(t, v_x_wheel, false, false);
 }

  if (var_set == id::EEMotionNodes(ee_)) {
	  Vector3d derivativeVec_motion = (f_tang_world + f_tang_abs_value * Dg);
	  Jacobian derivativeJac_motion = derivativeVec_motion.transpose().sparseView();

	  jac.row(row++) = derivativeJac_motion * ee_motion_->GetJacobianWrtNodes(t, kVel, false);
  }

  if (var_set == id::EEForceNodes(ee_)) {
      Vector3d derivativeVec_force = (v_total + v_tot_abs_value * Dg_f);
      Jacobian derivativeJac_force = derivativeVec_force.transpose().sparseView();

	  jac.row(row++) = derivativeJac_force * ee_force_->GetJacobianWrtNodes(t, kPos, false);
  }
}

} /* namespace xpp */
