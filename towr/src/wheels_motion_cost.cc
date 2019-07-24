/*
 * wheels_motion_cost.cc
 *
 *  Created on: Jul 20, 2019
 *      Author: vivian
 */

#include <towr/costs/wheels_motion_cost.h>
#include <towr/variables/variable_names.h>
#include <towr/variables/state.h>
#include <towr/models/endeffector_mappings.h>

#include <cmath>
#include <iostream>

namespace towr {

WheelsMotionCost::WheelsMotionCost (const double dt, double weight)
    : CostTerm("wheels_motion_cost")
{
  weight_ = weight;
  dt_ = dt;
}

void
WheelsMotionCost::InitVariableDependedQuantities (const VariablesPtr& x)
{
  LF_motion_nodes_  = x->GetComponent<NodesVariablesEEWheelsMotion>(id::EEWheelsMotionNodes(LF));
  RF_motion_nodes_  = x->GetComponent<NodesVariablesEEWheelsMotion>(id::EEWheelsMotionNodes(RF));
}

double
WheelsMotionCost::GetCost () const
{
  double cost = 0.0;
  std::vector<Node> lf_motion = LF_motion_nodes_->GetNodes();
  std::vector<Node> rf_motion = RF_motion_nodes_->GetNodes();
  for (int i = 0; i < lf_motion.size(); ++i ) {
    double lf_z = lf_motion[i].at(kPos).z();
    double rf_z = rf_motion[i].at(kPos).z();
    if (i == 0 || i == lf_motion.size()-1) {
      cost += pow(lf_z - rf_z, 2);
    }
    else {
      cost += 2*pow(lf_z - rf_z, 2);
    }
  }

  return weight_*cost*(dt_/2.0);
}

int
WheelsMotionCost::GetCol(int node, int dim) const
{
  int col = node*k6D + dim;

  return col;
}

void
WheelsMotionCost::FillJacobianBlock (std::string var_set, Jacobian& jac) const
{
  std::vector<Node> lf_motion = LF_motion_nodes_->GetNodes();
  std::vector<Node> rf_motion = RF_motion_nodes_->GetNodes();

  if (var_set == id::EEWheelsMotionNodes(LF)) {
    for (int i = 0; i < lf_motion.size(); ++i ) {
	  double lf_z = lf_motion[i].at(kPos).z();
	  double rf_z = rf_motion[i].at(kPos).z();
	  if (i == 0 || i == lf_motion.size()-1) {
	    jac.coeffRef(0, GetCol(i, Z)) = 2*(lf_z - rf_z)*weight_*dt_/2.0;
	  }
	  else {
		jac.coeffRef(0, GetCol(i, Z)) = 2*(lf_z - rf_z)*weight_*dt_;
	  }
    }
  }

  if (var_set == id::EEWheelsMotionNodes(RF)) {
    for (int i = 0; i < rf_motion.size(); ++i ) {
  	  double lf_z = lf_motion[i].at(kPos).z();
  	  double rf_z = rf_motion[i].at(kPos).z();
	  if (i == 0 || i == rf_motion.size()-1) {
	    jac.coeffRef(0, GetCol(i, Z)) = -2*(lf_z - rf_z)*weight_*dt_/2.0;
	  }
	  else {
		jac.coeffRef(0, GetCol(i, Z)) = -2*(lf_z - rf_z)*weight_*dt_;
	  }
    }
  }
}

} /* namespace towr */





