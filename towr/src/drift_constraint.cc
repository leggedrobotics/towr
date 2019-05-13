/*
 * drift_constraint.cc
 *
 *  Created on: Apr 25, 2019
 *      Author: dominic
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

DriftConstraint::Vector3d
DriftConstraint::NormDerivative (const Vector3d& v) const
{
  double p = pow(v.x(),2) + pow(v.y(),2) + pow(v.z(),2);

  if (p < 0.001)
	  return Vector3d(0,0,0);
  else
	  return Vector3d(v.x()/sqrt(p), v.y()/sqrt(p), v.z()/sqrt(p));
}

void
DriftConstraint::UpdateConstraintAtInstance (double t, int k, VectorXd& g) const
{
  EulerConverter::MatrixSXd b_C_w = base_angular_.GetRotationMatrixBaseToWorld(t);
  EulerConverter::MatrixSXd w_C_b = base_angular_.GetRotationMatrixBaseToWorld(t).transpose();

  //introduce a "fictional" constant spinning velocity of the wheels during drift
  //in the negative x-direction of the base projected onto the terrain
  double des_v_x_wheel = params_.vx_wheel_during_drift_;
  double pitch_angle = base_angular_for_angles_->GetPoint(t).p()(Y);
  Vector3d v_x_wheel(-des_v_x_wheel, 0.0, 0.0); //in base frame!!
  Vector3d v_x_wheel_world_with_z = b_C_w*v_x_wheel;
  Vector3d v_x_wheel_world(v_x_wheel_world_with_z(X), v_x_wheel_world_with_z(Y), 0.0);

  Vector3d f_tang_world(ee_force_->GetPoint(t).p()(X), ee_force_->GetPoint(t).p()(Y), 0.0);

  //total "fictional" velocity of the wheels
  Vector3d v_total = v_x_wheel_world + ee_motion_->GetPoint(t).v();

  double v_tot_abs_value =  v_total.norm();

  double f_tang_abs_value = f_tang_world.norm();

  int row = k*n_constraints_per_node_;

  //constrain the friction force to point in the opposite direction of the fictional total velocity
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
				  if ((t <= params_.ee_phase_durations_[0][0]) or (t > (params_.ee_phase_durations_[0][0]+params_.ee_phase_durations_[0][1]))){
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
    Vector3d v_x_wheel(-des_v_x_wheel, 0.0, 0.0); //in base frame!!!
    Vector3d v_x_wheel_world_with_z = b_C_w*v_x_wheel;
    Vector3d v_x_wheel_world(v_x_wheel_world_with_z(X), v_x_wheel_world_with_z(Y), 0.0);
    Vector3d v_x_wheel_base_with_y = w_C_b*v_x_wheel_world;
    Vector3d v_x_wheel_base(v_x_wheel_base_with_y(X), 0.0, 0.0);

    Vector3d v_total = v_x_wheel_world + ee_motion_->GetPoint(t).v() + Vector3d::Constant(1e-10);

    double v_tot_abs_value =  v_total.norm();

    double f_tang_abs_value = f_tang_world.norm();

  Vector3d Dg = NormDerivative(v_total) + Vector3d::Constant(1e-10);
  Jacobian g_deriv = Dg.transpose().sparseView();
  Vector3d Dg_f = NormDerivative(f_tang_world) + Vector3d::Constant(1e-10);
  Jacobian g_deriv_f = Dg_f.transpose().sparseView();

  if (var_set == id::base_ang_nodes) {

	  Vector3d tt_base = (f_tang_world + f_tang_abs_value * Dg);
	  Jacobian ttt_base = tt_base.transpose().sparseView();
//	  Vector3d d_vtot(0.0, 0.0, 0.0);
//	  d_vtot = d_vtot + Vector3d::Constant(1e-10);
//	  Jacobian d_vtot_jac = b_C_w*d_vtot; //sparseView()?

	  //only v_x_wheel is rotated (or v_x_wheel_base?)
	  jac.row(row++) = ttt_base * base_angular_.DerivOfRotVecMult(t, v_x_wheel, false, false);

	  // TODO: add the cosine
//	  jac.row(row++) = ttt_base * (base_angular_.DerivOfRotVecMult(t, v_x_wheel, false, false) + d_vtot_jac);
  }

  if (var_set == id::EEMotionNodes(ee_)) {
	  Vector3d tt_motion = (f_tang_world + f_tang_abs_value * Dg);
	  Jacobian ttt_motion = tt_motion.transpose().sparseView();

	  jac.row(row++) = ttt_motion * ee_motion_->GetJacobianWrtNodes(t, kVel, false);

//	  jac.row(row++) = f_tang_world * ee_motion_->GetJacobianWrtNodes(t, kVel, false) + f_tang_abs_value * g_deriv * ee_motion_->GetJacobianWrtNodes(t, kVel, false);
  }

  if (var_set == id::EEForceNodes(ee_)) {
      Vector3d tt = (v_total + v_tot_abs_value * Dg_f);
      Jacobian ttt = tt.transpose().sparseView();
	  jac.row(row++) = ttt * ee_force_->GetJacobianWrtNodes(t, kPos, false);

  }
}

} /* namespace xpp */
