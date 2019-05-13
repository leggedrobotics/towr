/*
 * wheels_motion_constraint.cc
 *
 *  Created on: Apr 5, 2019
 *      Author: vivian
 */
#include <iostream>
#include <towr/constraints/wheels_motion_constraint.h>
#include <towr/variables/variable_names.h>
#include <towr/variables/cartesian_dimensions.h>
#include <towr/variables/spline_holder.h>

namespace towr {


WheelsMotionConstraint::WheelsMotionConstraint (const KinematicModel::Ptr& model,
												double T, double dt, const EE& ee,
                                                const SplineHolderDrive& spline_holder)
    :TimeDiscretizationConstraint(T, dt, "wheels-motion-" + std::to_string(ee))
{
  base_linear_  = spline_holder.base_linear_;
  base_angular_ = EulerConverter(spline_holder.base_angular_);
  ee_      		= ee;

  ee_wheels_motion_ = spline_holder.ee_wheels_motion_.at(ee_);

  max_leg_extension_ = 0.50;
  min_leg_extension_ = 0.35;

  nominal_pos_hip_B_ 		  = model->GetHipPositionInBaseFrame().at(ee);
  max_deviation_from_nominal_ = model->GetMaximumDeviationFromNominal();
  nominal_ee_pos_B_           = model->GetNominalStanceInBase().at(ee);

  n_constraints_per_node_ = 3;  // 3 constraints per node: leg extension + (x,y) displacement

  SetRows(GetNumberOfNodes()*n_constraints_per_node_);
}

void
WheelsMotionConstraint::UpdateConstraintAtInstance (double t, int k,
                                                    VectorXd& g) const
{
  Vector3d base_W  = base_linear_->GetPoint(t).p();
  Vector3d pos_ee_W = ee_wheels_motion_->GetPoint(t).p();
  EulerConverter::MatrixSXd w_R_b = base_angular_.GetRotationMatrixBaseToWorld(t);
  Vector3d pos_hip_W = w_R_b*(nominal_pos_hip_B_) + base_W;

  Vector3d vector_base_to_ee_W = pos_ee_W - base_W;
  Vector3d vector_base_to_ee_B = w_R_b.transpose()*(vector_base_to_ee_W);

  int row = k*n_constraints_per_node_;
  g(row++) = (pos_ee_W - pos_hip_W).norm();
  g(row++) = vector_base_to_ee_B(X);
  g(row++) = vector_base_to_ee_B(Y);
}

void
WheelsMotionConstraint::UpdateBoundsAtInstance (double t, int k, VecBound& bounds) const
{
  int row = k*n_constraints_per_node_;

  bounds.at(row++) = ifopt::Bounds(min_leg_extension_,  max_leg_extension_);

  bounds.at(row++) = ifopt::Bounds(nominal_ee_pos_B_(X) - max_deviation_from_nominal_(X), nominal_ee_pos_B_(X) + max_deviation_from_nominal_(X));
  bounds.at(row++) = ifopt::Bounds(nominal_ee_pos_B_(Y) - max_deviation_from_nominal_(Y), nominal_ee_pos_B_(Y) + max_deviation_from_nominal_(Y));
}

WheelsMotionConstraint::Vector3d
WheelsMotionConstraint::NormDerivative (const Vector3d& v) const
{
  double q = pow(v.x(),2) + pow(v.y(),2) + pow(v.z(),2);

  // TODO: protection against division by 0
  return Vector3d(v.x()/sqrt(q), v.y()/sqrt(q), v.z()/sqrt(q));
}

int
WheelsMotionConstraint::GetCol(int node, int dim) const
{
  int col = node*k6D+dim;
  int n_opt_variables = ee_wheels_motion_->GetNodeVariablesCount();

  if (col >= n_opt_variables)
    col = (node-1)*k6D+dim;

  return col;
}

void
WheelsMotionConstraint::UpdateJacobianAtInstance (double t, int k,
                                                  std::string var_set,
                                                  Jacobian& jac) const
{
  Vector3d base_W  = base_linear_->GetPoint(t).p();
  Vector3d pos_ee_W = ee_wheels_motion_->GetPoint(t).p();
  EulerConverter::MatrixSXd w_R_b = base_angular_.GetRotationMatrixBaseToWorld(t);
  EulerConverter::MatrixSXd b_R_w = w_R_b.transpose();
  Vector3d pos_hip_W = w_R_b*(nominal_pos_hip_B_) + base_W;

  Vector3d Dg = NormDerivative(pos_ee_W - pos_hip_W) + Vector3d::Constant(1e-10);
  Jacobian g_deriv = Dg.transpose().sparseView();

  if (var_set == id::base_lin_nodes) {
    int row = k*n_constraints_per_node_;
    jac.row(row++) = -g_deriv * base_linear_->GetJacobianWrtNodes(t, kPos);
    jac.row(row++) = -1*(b_R_w*base_linear_->GetJacobianWrtNodes(t, kPos)).row(X);
    jac.row(row++) = -1*(b_R_w*base_linear_->GetJacobianWrtNodes(t, kPos)).row(Y);
  }

  if (var_set == id::base_ang_nodes) {
    int row = k*n_constraints_per_node_;
    jac.row(row++) = -g_deriv * base_angular_.DerivOfRotVecMult(t,nominal_pos_hip_B_, false);
    Vector3d r_W = pos_ee_W - base_W;
    jac.row(row++) = base_angular_.DerivOfRotVecMult(t,r_W, true).row(X);
    jac.row(row++) = base_angular_.DerivOfRotVecMult(t,r_W, true).row(Y);
  }

  if (var_set == id::EEWheelsMotionNodes(ee_)) {
	int row = k*n_constraints_per_node_;
	jac.row(row++) = g_deriv * ee_wheels_motion_->GetJacobianWrtNodes(t, kPos);
	jac.row(row++) = (b_R_w*ee_wheels_motion_->GetJacobianWrtNodes(t,kPos)).row(X);
	jac.row(row++) = (b_R_w*ee_wheels_motion_->GetJacobianWrtNodes(t,kPos)).row(Y);
  }
}

} /* namespace towr */




