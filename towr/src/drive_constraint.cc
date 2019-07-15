/*
 * drive_constraint.cc
 *
 *  Created on: Apr 25, 2019
 *      Author: dominic landolf
 */

#include <towr/constraints/drive_constraint.h>
#include <towr/variables/variable_names.h>
#include <towr/variables/spline_holder.h>
#include <towr/variables/cartesian_dimensions.h>
#include <iostream>

using namespace std;

namespace towr {

// constrains the lateral velocity during drive phases to be 0
// constrains the heading velocity to be positive to avoid oscillations
// optional: constrains the lateral Acceleration to lie inside bounds to avoid infeasible motions.

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

  // set is_lateralAccBounds true in params, if y-Acc wants to be constrained
  if (!params_.is_lateralAccBounds)
	  n_constraints_per_node_ = 2;
  else
	  n_constraints_per_node_ = 3;

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

  EulerConverter::MatrixSXd w_C_b = base_angular_.GetRotationMatrixBaseToWorld(t).transpose();

  int poly_id = ee_motion_->GetSegmentID(t, T_);


  Vector3d v_wrt_b = w_C_b*ee_motion_->GetPoint(t).v();
  Vector3d a_wrt_b = w_C_b*ee_motion_->GetPoint(poly_id, T_.at(poly_id)).a();

  int row = k*n_constraints_per_node_;

  //only positive wheel velocity in x-direction (heading direction) to avoid oscillations of EE
  g(row++) = v_wrt_b(X);

  //no velocity in y-direction (lateral direction) during drive phases
  g(row++) = v_wrt_b(Y);

  //optional constraint: bounds for lateral acceleration
  if (params_.is_lateralAccBounds)
	  g(row++) = a_wrt_b(Y);
}

void
DriveConstraint::UpdateBoundsAtInstance (double t, int k, VecBound& bounds) const
{
	int row = k*n_constraints_per_node_;

	if (params_.just_drive_){
		bounds.at(row++) = ifopt::NoBound;
		bounds.at(row++) = ifopt::BoundZero;
		if (params_.is_lateralAccBounds)
			bounds.at(row++) = ifopt::NoBound;
	}
	else {

		if (ee_ == 0 or ee_ == 1){
				bounds.at(row++) = ifopt::BoundGreaterZero;
				bounds.at(row++) = ifopt::BoundZero;
				if (params_.is_lateralAccBounds)
					bounds.at(row++) = ifopt::NoBound;
		}
		else if (ee_ == 2 or ee_ == 3){
				if ((t <= params_.phase_duration_drive_1) or (t > (params_.phase_duration_drive_1+params_.phase_duration_drift))){
					  bounds.at(row++) = ifopt::BoundGreaterZero;
					  bounds.at(row++) = ifopt::BoundZero;
					  if (params_.is_lateralAccBounds)
						  bounds.at(row++) = ifopt::NoBound;
				 }
				 else {
					  bounds.at(row++) = ifopt::BoundGreaterZero;
					  bounds.at(row++) = ifopt::NoBound;
					  if (params_.is_lateralAccBounds)
						  bounds.at(row++) = ifopt::Bounds(-params_.lateral_AccBound, params_.lateral_AccBound);
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

  int row = k*n_constraints_per_node_;

  if (var_set == id::base_ang_nodes) {

    int poly_id = ee_motion_->GetSegmentID(t, T_);

    Vector3d v_W = ee_motion_->GetPoint(t).v();
    Vector3d a_W = ee_motion_->GetPoint(poly_id, T_.at(poly_id)).a();

    jac.row(row++) = base_angular_.DerivOfRotVecMult(t, v_W, true, false).row(X);
    jac.row(row++) = base_angular_.DerivOfRotVecMult(t, v_W, true, false).row(Y);
    if (params_.is_lateralAccBounds)
    	jac.row(row++) = base_angular_.DerivOfRotVecMult(t, a_W, true, false).row(Y);
  }

  if (var_set == id::EEMotionNodes(ee_)) {
	  int poly_id = ee_motion_->GetSegmentID(t, T_);
	  jac.row(row++) = (b_R_w*ee_motion_->GetJacobianWrtNodes(t, kVel, false)).row(X);
	  jac.row(row++) = (b_R_w*ee_motion_->GetJacobianWrtNodes(t, kVel, false)).row(Y);
	  if (params_.is_lateralAccBounds)
		  jac.row(row++) = (b_R_w*ee_motion_->GetJacobianWrtNodes(poly_id, T_.at(poly_id), kAcc, false)).row(Y);
  }
}

} /* namespace xpp */
