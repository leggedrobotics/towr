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

  mu_   = 0.5;

  base_linear_  = spline_holder.base_linear_;
  base_angular_ = EulerConverter(spline_holder.base_angular_);
  base_angular_for_angles_ = spline_holder.base_angular_;

  ee_ = ee;
  ee_motion_    = spline_holder.ee_motion_.at(ee_);
  ee_force_ 	= spline_holder.ee_force_.at(ee_);


  n_constraints_per_node_ = 2;

  SetRows(GetNumberOfNodes()*n_constraints_per_node_);

//	  T_ = ee_motion_->GetPolyDurations();
}

int
DriftConstraint::GetRow (int node, int dim) const
{
  return node*n_constraints_per_node_ + dim;
}

void
DriftConstraint::UpdateConstraintAtInstance (double t, int k, VectorXd& g) const
{
//  EulerConverter::MatrixSXd b_R_w = base_angular_.GetRotationMatrixBaseToWorld(t).transpose();
  EulerConverter::MatrixSXd b_C_w = base_angular_.GetRotationMatrixBaseToWorld(t);
  EulerConverter::MatrixSXd w_C_b = base_angular_.GetRotationMatrixBaseToWorld(t).transpose();

  double des_v_x_wheel = params_.vx_wheel_during_drift_;
  double pitch_angle = base_angular_for_angles_->GetPoint(t).p()(Y);

  Vector3d v_x_wheel(-des_v_x_wheel/cos(pitch_angle), 0.0, 0.0); //in base frame!!!
  Vector3d v_x_wheel_world_with_z = b_C_w*v_x_wheel;
  Vector3d v_x_wheel_world(v_x_wheel_world_with_z(X), v_x_wheel_world_with_z(Y), 0.0);
  Vector3d f_tang_world(ee_force_->GetPoint(t).p()(X), ee_force_->GetPoint(t).p()(Y), 0.0);

  Vector3d v_total = v_x_wheel_world + ee_motion_->GetPoint(t).v();

  double v_tot_abs_value =  v_total.squaredNorm();
  Vector3d direction_vtot(0.0,0.0,0.0);
  if (v_tot_abs_value > 0.001)
  	  direction_vtot = v_total/v_tot_abs_value;

  double f_tang_abs_value = f_tang_world.squaredNorm();
  Vector3d direction_f_tang(0.0,0.0,0.0);
  if (f_tang_abs_value > 0.001)
	  direction_f_tang = f_tang_world/f_tang_abs_value;

//  Vector3d n(0,0,1);
//  Vector3d t1(1,0,0);
//  Vector3d t2(0,1,0);
//  double f_tang = sqrt(pow(f.transpose()*t1, 2) + pow(f.transpose()*t2, 2));

  int row = k*n_constraints_per_node_;

  g(row++) = direction_vtot(X) + direction_f_tang(X);		// set x and y direction (negative) equal
  g(row++) = direction_vtot(Y) + direction_f_tang(Y);

}

void
DriftConstraint::UpdateBoundsAtInstance (double t, int k, VecBound& bounds) const
{
	int row = k*n_constraints_per_node_;

	if (params_.just_drive_){
		bounds.at(row++) = ifopt::NoBound;
		bounds.at(row++) = ifopt::NoBound;
	}
	else {
		if (ee_ == 0 or ee_ == 1){
				bounds.at(row++) = ifopt::NoBound;
				bounds.at(row++) = ifopt::NoBound;
			  }
			  else if (ee_ == 2 or ee_ == 3){
				  if ((t <= params_.ee_phase_durations_[0][0]) or (t > (params_.ee_phase_durations_[0][0]+params_.ee_phase_durations_[0][1]))){
					  bounds.at(row++) = ifopt::NoBound;
					  bounds.at(row++) = ifopt::NoBound;
				  }
				  //this is the drift phase
				  else {
					  bounds.at(row++) = ifopt::BoundZero;
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

  Vector3d v_wrt_b = w_C_b*ee_motion_->GetPoint(t).v();
    Vector3d f_wrt_b = w_C_b*ee_force_->GetPoint(t).p();
    Vector3d f_tang_world(ee_force_->GetPoint(t).p()(X), ee_force_->GetPoint(t).p()(Y), 0.0);
    Vector3d f_tang_wrt_b(f_wrt_b(X), f_wrt_b(Y), 0.0);
    double pitch_angle = base_angular_for_angles_->GetPoint(t).p()(Y);

    double des_v_x_wheel = params_.vx_wheel_during_drift_;
    Vector3d v_x_wheel(-des_v_x_wheel/cos(pitch_angle), 0.0, 0.0); //in base frame!!! //caution: this is also in world z-direction, need to adjust that...
    Vector3d v_x_wheel_world_with_z = b_C_w*v_x_wheel;
    Vector3d v_x_wheel_world(v_x_wheel_world_with_z(X), v_x_wheel_world_with_z(Y), 0.0);
    Vector3d v_x_wheel_base_with_y = w_C_b*v_x_wheel_world;
    Vector3d v_x_wheel_base(v_x_wheel_base_with_y(X), 0.0, 0.0);
//    Vector3d v_total = v_x_wheel + v_wrt_b;

    Vector3d v_total = v_x_wheel_world + ee_motion_->GetPoint(t).v();

    double v_tot_abs_value =  v_total.squaredNorm();
    Vector3d direction_vtot(0.0,0.0,0.0);
    if (v_tot_abs_value > 0.001)
    	  direction_vtot = v_total/v_tot_abs_value;

    double f_tang_abs_value = f_tang_world.squaredNorm();

    Vector3d direction_f_tang(0.0,0.0,0.0);
    if (f_tang_abs_value > 0.001)
  	  direction_f_tang = f_tang_world/f_tang_abs_value;

    Vector3d Dg = direction_vtot + Vector3d::Constant(1e-10); //auch norm von force beachten, JA!?
    Vector3d Df = direction_f_tang + Vector3d::Constant(1e-10);
    Jacobian g_deriv = Dg.transpose().sparseView();
    Jacobian f_deriv = Df.transpose().sparseView();

//    Vector3d x_dir(1,0,0);
//    Vector3d y_dir(0,1,0);
//    Jacobian x_dir_jac = x_dir.transpose().sparseView();
//    Jacobian y_dir_jac = y_dir.transpose().sparseView();

  if (var_set == id::base_ang_nodes) {

//v_tot_abs_value ist auch rotiert...!?

//	  Vector3d v_W(0.0,0.0,0.0);
//	  Vector3d f_W(0.0,0.0,0.0);
//	  if (v_tot_abs_value != 0)
//		  v_W = (ee_motion_->GetPoint(t).v() + v_x_wheel_world)/v_tot_abs_value;
////		  v_W = (ee_motion_->GetPoint(t).v() + v_x_wheel)/v_tot_abs_value;
//	  if (f_tang_abs_value != 0)
//		  f_W = f_tang_world/f_tang_abs_value;
//
//	  Vector3d v_f_W = v_W + f_W;

	  Vector3d v_wheel_norm = v_x_wheel_base;
	  if (v_tot_abs_value > 0.001)
		  v_wheel_norm = v_x_wheel_base/v_tot_abs_value;


//    	jac.row(row++) = base_angular_.DerivOfRotVecMult(t, v_W, true, false).row(X) + base_angular_.DerivOfRotVecMult(t, f_W, true, false).row(X);
//    	jac.row(row++) = base_angular_.DerivOfRotVecMult(t, v_W, true, false).row(Y) + base_angular_.DerivOfRotVecMult(t, f_W, true, false).row(Y);

	  //forces nicht beachtet...!! (ist aber auch unabhaenging von base angular...?)
	  if (v_tot_abs_value > 0.001) {
		  jac.row(row++) = (base_angular_.DerivOfRotVecMult(t, v_x_wheel_base, false, false).row(X)*v_tot_abs_value - v_total(X)*g_deriv)/(v_tot_abs_value*v_tot_abs_value);
		  jac.row(row++) = (base_angular_.DerivOfRotVecMult(t, v_x_wheel_base, false, false).row(Y)*v_tot_abs_value - v_total(Y)*g_deriv)/(v_tot_abs_value*v_tot_abs_value);

//		  jac.row(row++) = (base_angular_.DerivOfRotVecMult(t, v_x_wheel_base, false, false).row(X)*v_tot_abs_value - v_total(X)*(v_total(X)*base_angular_.DerivOfRotVecMult(t, v_x_wheel_base, false, false).row(X))/v_tot_abs_value)/(v_tot_abs_value*v_tot_abs_value);
//		  jac.row(row++) = (base_angular_.DerivOfRotVecMult(t, v_x_wheel_base, false, false).row(Y)*v_tot_abs_value - v_total(Y)*(v_total(Y)*base_angular_.DerivOfRotVecMult(t, v_x_wheel_base, false, false).row(Y))/v_tot_abs_value)/(v_tot_abs_value*v_tot_abs_value);

//		  jac.row(row++) = base_angular_.DerivOfRotVecMult(t, v_wheel_norm, false, false).row(X);
//		  jac.row(row++) = base_angular_.DerivOfRotVecMult(t, v_wheel_norm, false, false).row(Y);
	  }
	  else {
		  jac.row(row++) = base_angular_.DerivOfRotVecMult(t, v_x_wheel_base, false, false).row(X);
		  jac.row(row++) = base_angular_.DerivOfRotVecMult(t, v_x_wheel_base, false, false).row(Y);
	  }
  }

  if (var_set == id::EEMotionNodes(ee_)) {

	  if (v_tot_abs_value > 0.001){
//		  jac.row(row++) = (ee_motion_->GetJacobianWrtNodes(t, kVel, false).row(X)*v_tot_abs_value - v_total(X)*(v_total(X)*ee_motion_->GetJacobianWrtNodes(t, kVel, false).row(X))/v_tot_abs_value)/(v_tot_abs_value*v_tot_abs_value);
//		  jac.row(row++) = (ee_motion_->GetJacobianWrtNodes(t, kVel, false).row(Y)*v_tot_abs_value - v_total(Y)*(v_total(Y)*ee_motion_->GetJacobianWrtNodes(t, kVel, false).row(Y))/v_tot_abs_value)/(v_tot_abs_value*v_tot_abs_value);

		  jac.row(row++) = (ee_motion_->GetJacobianWrtNodes(t, kVel, false).row(X)*v_tot_abs_value - v_total(X)*g_deriv)/(v_tot_abs_value*v_tot_abs_value);
		  jac.row(row++) = (ee_motion_->GetJacobianWrtNodes(t, kVel, false).row(Y)*v_tot_abs_value - v_total(Y)*g_deriv)/(v_tot_abs_value*v_tot_abs_value);

	  }
	  else {
		  jac.row(row++) = ee_motion_->GetJacobianWrtNodes(t, kVel, false).row(X);
		  jac.row(row++) = ee_motion_->GetJacobianWrtNodes(t, kVel, false).row(Y);
	  }
	  }

  if (var_set == id::EEForceNodes(ee_)) {

	  if (f_tang_abs_value > 0.001){
//		  jac.row(row++) = (ee_force_->GetJacobianWrtNodes(t, kPos, false).row(X)*f_tang_abs_value - f_tang_world(X)*f_deriv)/(f_tang_abs_value*f_tang_abs_value);
//		  jac.row(row++) = (ee_force_->GetJacobianWrtNodes(t, kPos, false).row(Y)*f_tang_abs_value - f_tang_world(Y)*f_deriv)/(f_tang_abs_value*f_tang_abs_value);

		  jac.row(row++) = (ee_force_->GetJacobianWrtNodes(t, kPos, false).row(X)*f_tang_abs_value - f_tang_world(X)*f_tang_world(X)*ee_force_->GetJacobianWrtNodes(t, kPos, false).row(X)/f_tang_abs_value)/(f_tang_abs_value*f_tang_abs_value);
		  jac.row(row++) = (ee_force_->GetJacobianWrtNodes(t, kPos, false).row(Y)*f_tang_abs_value - f_tang_world(Y)*f_tang_world(Y)*ee_force_->GetJacobianWrtNodes(t, kPos, false).row(Y)/f_tang_abs_value)/(f_tang_abs_value*f_tang_abs_value);

//		  jac.row(row++) = (x_dir_jac*ee_force_->GetJacobianWrtNodes(t, kPos, false)*f_tang_abs_value - f_tang_world(X)*f_tang_world(X)*x_dir_jac/f_tang_abs_value)/(f_tang_abs_value*f_tang_abs_value);
//		  jac.row(row++) = (y_dir_jac*ee_force_->GetJacobianWrtNodes(t, kPos, false)*f_tang_abs_value - f_tang_world(Y)*f_tang_world(Y)*y_dir_jac/f_tang_abs_value)/(f_tang_abs_value*f_tang_abs_value);

//		  jac.row(row++) = ee_force_->GetJacobianWrtNodes(t, kPos, false).row(X)/(f_tang_abs_value*f_tang_abs_value);
//		  jac.row(row++) = ee_force_->GetJacobianWrtNodes(t, kPos, false).row(Y)/(f_tang_abs_value*f_tang_abs_value);
	  }
	  else {
		  jac.row(row++) = ee_force_->GetJacobianWrtNodes(t, kPos, false).row(X);
		  jac.row(row++) = ee_force_->GetJacobianWrtNodes(t, kPos, false).row(Y);
	  }
  	  }
}

} /* namespace xpp */
