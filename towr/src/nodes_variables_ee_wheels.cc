/*
 * nodes_variables_ee_wheels.cc
 *
 *  Created on: Mar 29, 2019
 *      Author: vivian
 */

#include <iostream>
#include <towr/variables/nodes_variables_ee_wheels.h>
#include <towr/variables/cartesian_dimensions.h>

namespace towr {

// EE Motion
NodesVariablesEEWheelsMotion::NodesVariablesEEWheelsMotion (int n_nodes, int n_dim, std::string variable_id, const HeightMap::Ptr& terrain)
    : NodesVariablesAll(n_nodes, n_dim, variable_id)
{
  terrain_ = terrain;
}

void
NodesVariablesEEWheelsMotion::SetByLinearInterpolation(const VectorXd& initial_val,
                                         	 	 	   const VectorXd& final_val,
													   double t_total)
{
  VectorXd dp = final_val-initial_val;
  VectorXd average_velocity = dp / t_total;
  int num_nodes = nodes_.size();

  for (int idx=0; idx<GetRows(); ++idx) {
    for (auto nvi : GetNodeValuesInfo(idx)) {

      if (nvi.deriv_ == kPos) {
        VectorXd pos = initial_val + nvi.id_/static_cast<double>(num_nodes-1)*dp;
        pos(Z) = terrain_->GetHeight(pos(X),pos(Y)); // enforces terrain contact!
        nodes_.at(nvi.id_).at(kPos)(nvi.dim_) = pos(nvi.dim_);
      }

      if (nvi.deriv_ == kVel) {
        nodes_.at(nvi.id_).at(kVel)(nvi.dim_) = average_velocity(nvi.dim_);
      }
    }
  }
}

// EE force
NodesVariablesEEWheelsForce::NodesVariablesEEWheelsForce (int n_nodes, int n_dim, std::string variable_id)
    : NodesVariablesAll(n_nodes, n_dim, variable_id) {}

void
NodesVariablesEEWheelsForce::SetAllAtConstantValue (const VectorXd& force)
{
  for (int idx=0; idx<GetRows(); ++idx) {
    for (auto nvi : GetNodeValuesInfo(idx)) {
      if (nvi.deriv_ == kPos) {
	    nodes_.at(nvi.id_).at(kPos)(nvi.dim_) = force(nvi.dim_);
      }

//      if (nvi.deriv_ == kVel) {
//        nodes_.at(nvi.id_).at(kVel)(nvi.dim_) = 0.0;
//      }
    }
  }
}


} /* namespace towr */
