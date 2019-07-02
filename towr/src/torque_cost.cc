/*
 * torque_cost.cc
 *
 *  Created on: Jul 1, 2019
 *      Author: vivian
 */

#include <towr/costs/torque_cost.h>
#include <towr/variables/variable_names.h>
#include <towr/variables/state.h>

#include <cmath>
#include <iostream>

namespace towr {

TorqueCost::TorqueCost (const HeightMap::Ptr& terrain, double weight, const uint& ee)
    : CostTerm("wheels_torque_cost_" + std::to_string(ee))
{
  terrain_ = terrain;
  weight_ = weight;
  ee_ = ee;
}

void
TorqueCost::InitVariableDependedQuantities (const VariablesPtr& x)
{
  ee_force_nodes_  = x->GetComponent<NodesVariablesEEWheelsForce>(id::EEWheelsForceNodes(ee_));
  ee_motion_nodes_ = x->GetComponent<NodesVariablesEEWheelsMotion>(id::EEWheelsMotionNodes(ee_));
}

double
TorqueCost::GetCost () const
{
  double cost = 0.0;
  std::vector<Node> f = ee_force_nodes_->GetNodes();
  std::vector<Node> p = ee_motion_nodes_->GetNodes();
  for (int i = 0; i < f.size(); ++i ) {
    Vector3d ee_force = f[i].at(kPos);
    Vector3d ee_pos   = p[i].at(kPos);
    Vector3d t1 = terrain_->GetNormalizedBasis(HeightMap::Tangent1, ee_pos.x(), ee_pos.y());
    double ft = ee_force.transpose() * t1;
    cost += weight_*std::pow(ft,2);
  }

  return cost;
}

int
TorqueCost::GetCol(int node, int dim) const
{
  int col = node*k6D + dim;

  return col;
}

void
TorqueCost::FillJacobianBlock (std::string var_set, Jacobian& jac) const
{
  std::vector<Node> f = ee_force_nodes_->GetNodes();
  std::vector<Node> p = ee_motion_nodes_->GetNodes();

  if (var_set == id::EEWheelsForceNodes(ee_)) {
    for (int i = 0; i < f.size(); ++i ) {
	  Vector3d ee_force = f[i].at(kPos);
	  Vector3d ee_pos   = p[i].at(kPos);
	  Vector3d t1 = terrain_->GetNormalizedBasis(HeightMap::Tangent1, ee_pos.x(), ee_pos.y());
	  double ft = ee_force.transpose() * t1;
	  for (auto dim : {X, Y, Z})
		jac.coeffRef(0, GetCol(i, dim)) = 2*weight_*ft*t1(dim);
    }
  }

  if (var_set == id::EEWheelsMotionNodes(ee_)) {
    for (int i = 0; i < f.size(); ++i ) {
	  Vector3d ee_force = f[i].at(kPos);
	  Vector3d ee_pos   = p[i].at(kPos);
	  Vector3d t1 = terrain_->GetNormalizedBasis(HeightMap::Tangent1, ee_pos.x(), ee_pos.y());
	  double ft = ee_force.transpose() * t1;
	  for (auto dim : {X_, Y_}) {
		Vector3d dt1 = terrain_->GetDerivativeOfNormalizedBasisWrt(HeightMap::Tangent1, dim, ee_pos.x(), ee_pos.y());
		double aux = ee_force.transpose() * dt1;
		jac.coeffRef(0, GetCol(i, dim)) = 2*weight_*ft*aux;
	  }
    }
  }
}

} /* namespace towr */


