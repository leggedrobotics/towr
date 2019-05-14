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

  // set 3 if y-vel. and y-acc constraint, 2 if only y-vel constraint (and pos x vel)
  n_constraints_per_node_ = 2; // 3: if y-vel. and y-acc constraint

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

//  Vector3d v_y = b_C_w*v_b_y;

  int row = k*n_constraints_per_node_;

  g(row++) = v_wrt_b(X);
  g(row++) = v_wrt_b(Y);
  if (n_constraints_per_node_ == 3)
	  g(row++) = a_b_y(Y);
}

void
DriveConstraint::UpdateBoundsAtInstance (double t, int k, VecBound& bounds) const
{
	int row = k*n_constraints_per_node_;

	if (params_.just_drive_){
		bounds.at(row++) = ifopt::BoundGreaterZero;
		bounds.at(row++) = ifopt::BoundZero;
		if (n_constraints_per_node_ == 3)
			bounds.at(row++) = ifopt::BoundZero;
	}
	else {
		if (ee_ == 0 or ee_ == 1){
				bounds.at(row++) = ifopt::BoundGreaterZero;
				bounds.at(row++) = ifopt::BoundZero;
				if (n_constraints_per_node_ == 3)
					bounds.at(row++) = ifopt::BoundZero;
			  }
			  else if (ee_ == 2 or ee_ == 3){
				  if ((t <= params_.ee_phase_durations_[0][0]) or (t > (params_.ee_phase_durations_[0][0]+params_.ee_phase_durations_[0][1]))){
					  bounds.at(row++) = ifopt::BoundGreaterZero;
					  bounds.at(row++) = ifopt::BoundZero;
					  if (n_constraints_per_node_ == 3)
						  bounds.at(row++) = ifopt::BoundZero;
				  }
				  else {
					  bounds.at(row++) = ifopt::BoundGreaterZero; //during drift also just in positive x dir
					  bounds.at(row++) = ifopt::NoBound;
					  if (n_constraints_per_node_ == 3)
						  bounds.at(row++) = ifopt::NoBound;
				  }
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

  if (var_set == id::base_ang_nodes) {

    int poly_id = ee_motion_->GetSegmentID(t, T_);

    Vector3d v_W = ee_motion_->GetPoint(t).v();
    Vector3d a_W = ee_motion_->GetPoint(poly_id, T_.at(poly_id)).a(); //k is not poly id!

    	jac.row(row++) = base_angular_.DerivOfRotVecMult(t, v_W, true, false).row(X);
    	jac.row(row++) = base_angular_.DerivOfRotVecMult(t, v_W, true, false).row(Y);
    	if (n_constraints_per_node_ == 3)
    		jac.row(row++) = base_angular_.DerivOfRotVecMult(t, a_W, true, false).row(Y);

  }

  if (var_set == id::EEMotionNodes(ee_)) {
	  int poly_id = ee_motion_->GetSegmentID(t, T_);
	  	  jac.row(row++) = (b_R_w*ee_motion_->GetJacobianWrtNodes(t, kVel, false)).row(X);
		  jac.row(row++) = (b_R_w*ee_motion_->GetJacobianWrtNodes(t, kVel, false)).row(Y);
		  if (n_constraints_per_node_ == 3)
			  jac.row(row++) = (b_R_w*ee_motion_->GetJacobianWrtNodes(poly_id, T_.at(poly_id), kAcc, false)).row(Y);

	  }
}

} /* namespace xpp */
