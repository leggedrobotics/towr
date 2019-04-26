/*
 * drive_constraint.cc
 *
 *  Created on: Apr 25, 2019
 *      Author: dominic
 */

#include <towr/constraints/drive_constraint.h>
#include <towr/variables/variable_names.h>
#include <towr/variables/spline_holder.h>
#include <towr/variables/cartesian_dimensions.h>

namespace towr {

DriveConstraint::DriveConstraint (const KinematicModel::Ptr& model,
                                                  double T, double dt,
                                                  const EE& ee,
                                                  const SplineHolder& spline_holder)
    :TimeDiscretizationConstraint(T, dt, "drive_motion-" + std::to_string(ee))
{
  base_linear_  = spline_holder.base_linear_;
  base_angular_ = EulerConverter(spline_holder.base_angular_);
  ee_motion_    = spline_holder.ee_motion_.at(ee);

  ee_ = ee;
  n_constraints_per_node_ = 2; // 2: y-vel. and y-acc

	  SetRows(GetNumberOfNodes()*n_constraints_per_node_);

	  T_ = ee_motion_->GetPolyDurations();
}

int
DriveConstraint::GetRow (int node, int dim) const
{
  return node*n_constraints_per_node_ + dim;
}

void
DriveConstraint::UpdateConstraintAtInstance (double t, int k, VectorXd& g) const
{
//  EulerConverter::MatrixSXd b_R_w = base_angular_.GetRotationMatrixBaseToWorld(t).transpose();
  //EulerConverter::MatrixSXd b_C_w = base_angular_.GetRotationMatrixBaseToWorld(t);
  EulerConverter::MatrixSXd w_C_b = base_angular_.GetRotationMatrixBaseToWorld(t).transpose();

  int poly_id = ee_motion_->GetSegmentID(t, T_);


  Vector3d v_wrt_b = w_C_b*ee_motion_->GetPoint(t).v();
  Vector3d a_wrt_b = w_C_b*ee_motion_->GetPoint(poly_id, T_.at(poly_id)).a(); //k is not poly id!
  Vector3d v_b_y = {0, v_wrt_b(1), 0};
  Vector3d a_b_y = {0, a_wrt_b(1), 0};
//  Vector3d v_b_x = {v_wrt_b(0), 0, 0};
//  Vector3d v_x = b_C_w*v_b_x;
//  Vector3d v_y = b_C_w*v_b_y;

  int row = k*n_constraints_per_node_;

  g(row++) = v_b_y(Y);
  g(row++) = a_b_y(Y);

}

void
DriveConstraint::UpdateBoundsAtInstance (double t, int k, VecBound& bounds) const
{
	int row = k*n_constraints_per_node_;

	if (ee_ == 0 or ee_ == 1){
		bounds.at(row++) = ifopt::BoundZero;
		bounds.at(row++) = ifopt::BoundZero;
	  }
	  else if (ee_ == 2 or ee_ == 3){
		  if ((t <= params_.ee_phase_durations_[0][0]) or (t > (params_.ee_phase_durations_[0][0]+params_.ee_phase_durations_[0][1]))){
			  bounds.at(row++) = ifopt::BoundZero;
			  bounds.at(row++) = ifopt::BoundZero;
		  }
		  else {
			  bounds.at(row++) = ifopt::NoBound;
		  	  bounds.at(row++) = ifopt::NoBound;
		  }
		  }
}

void
DriveConstraint::UpdateJacobianAtInstance (double t, int k,
                                                   std::string var_set,
                                                   Jacobian& jac) const
{
  EulerConverter::MatrixSXd b_R_w = base_angular_.GetRotationMatrixBaseToWorld(t).transpose();
//  EulerConverter::MatrixSXd w_C_b = base_angular_.GetRotationMatrixBaseToWorld(t).transpose();

  int row = k*n_constraints_per_node_;

//  if (var_set == id::base_lin_nodes) {
//
//		  jac.row(row++) = -1*(b_R_w*base_linear_->GetJacobianWrtNodes(t, kVel, false)).row(Y); //kPos bec only pos of base is of interest!
//		  jac.row(row++) = -1*(b_R_w*base_linear_->GetJacobianWrtNodes(t, kAcc, false)).row(Y);
//
//  }

  if (var_set == id::base_ang_nodes) {

//    Vector3d base_W   = base_linear_->GetPoint(t).p();
//    Vector3d ee_pos_W = ee_motion_->GetPoint(t).p();
//    Vector3d r_W = ee_pos_W - base_W;
//
//    Vector3d v_base_W   = base_linear_->GetPoint(t).v();
//    Vector3d ee_v_W = ee_motion_->GetPoint(t).v();
////    Vector3d ee_v_W = ee_motion_->GetPoint(poly_id, T_.at(poly_id)).v();
//    Vector3d v_W = ee_v_W - v_base_W;
//
    int poly_id = ee_motion_->GetSegmentID(t, T_);
//    Vector3d a_base_W   = base_linear_->GetPoint(t).a();
//    Vector3d ee_a_W = ee_motion_->GetPoint(poly_id, T_.at(poly_id)).a();
//    Vector3d a_W = ee_a_W - a_base_W;

    Vector3d v_W = ee_motion_->GetPoint(t).v();
    Vector3d a_W = ee_motion_->GetPoint(poly_id, T_.at(poly_id)).a(); //k is not poly id!

    	jac.row(row++) = base_angular_.DerivOfRotVecMult(t, v_W, true, false).row(Y);
    	jac.row(row++) = base_angular_.DerivOfRotVecMult(t, a_W, true, false).row(Y);

//    	jac.row(row++) = base_angular_.DerivOfRotVecMult(t, r_W, true, false).row(Y); //do I need v_W?
//    	jac.row(row++) = base_angular_.DerivOfRotVecMult(t, r_W, true, false).row(Y); //do I need a_W?

  }

  if (var_set == id::EEMotionNodes(ee_)) {
	  int poly_id = ee_motion_->GetSegmentID(t, T_);

		  jac.row(row++) = (b_R_w*ee_motion_->GetJacobianWrtNodes(t, kVel, false)).row(Y);
		  jac.row(row++) = (b_R_w*ee_motion_->GetJacobianWrtNodes(poly_id, T_.at(poly_id), kAcc, false)).row(Y);

	  }

//  if (var_set == id::EESchedule(ee_)) {
//
//		  jac.row(row++) = b_R_w*ee_motion_->GetJacobianOfPosWrtDurations(t);
//		  jac.row(row++) = b_R_w*ee_motion_->GetJacobianOfPosWrtDurations(t);
//
//  }
}

} /* namespace xpp */
