/*
 * base_acc_limits_constraint.cc
 *
 *  Created on: Apr 3, 2019
 *      Author: vivian
 */
#include <iostream>
#include <towr/constraints/base_acc_limits_constraint.h>
#include <towr/variables/variable_names.h>

namespace towr {

BaseAccLimitsConstraint::BaseAccLimitsConstraint (const std::vector<Vector3d>& acc_max,
	  	  	   	   	   	   	   	   	   	   	   	  double T, double dt,
												  const SplineHolder& spline_holder)
	: TimeDiscretizationConstraint(T, dt, "base-acc-limits")
{
  acc_max_lin_ = acc_max.at(0);
  acc_max_ang_ = acc_max.at(1);

  base_linear_  = spline_holder.base_linear_;
  base_angular_ = spline_holder.base_angular_;

  SetRows(GetNumberOfNodes()*k6D);  // 6 constraints per node (ax, ay, az, roll, pitch, yaw)
}

int
BaseAccLimitsConstraint::GetRow (int node, int dim) const
{
  return node*k6D + dim;
}

void
BaseAccLimitsConstraint::UpdateConstraintAtInstance (double t, int k, VectorXd& g) const
{
  g.middleRows(GetRow(k, LX), k3D) = base_linear_->GetPoint(t).a();
  g.middleRows(GetRow(k, AX), k3D) = base_angular_->GetPoint(t).a();
}

void
BaseAccLimitsConstraint::UpdateBoundsAtInstance (double t, int k, VecBound& bounds) const
{
  for (int dim : {LX, LY, LZ}) {
	bounds.at(GetRow(k,dim)) = ifopt::Bounds(-acc_max_lin_(dim-k3D), acc_max_lin_(dim-k3D));
  }

  for (int dim : {AX, AY, AZ}) {
	bounds.at(GetRow(k,dim)) = ifopt::Bounds(-acc_max_ang_(dim), acc_max_ang_(dim));
  }
}

void
BaseAccLimitsConstraint::UpdateJacobianAtInstance(double t, int k, std::string var_set, Jacobian& jac) const
{
  if (var_set == id::base_ang_nodes)
	jac.middleRows(GetRow(k,AX), k3D) = base_angular_->GetJacobianWrtNodes(t, kAcc);

  if (var_set == id::base_lin_nodes)
	jac.middleRows(GetRow(k,LX), k3D) = base_linear_->GetJacobianWrtNodes(t, kAcc);
}

} /* namespace towr */



