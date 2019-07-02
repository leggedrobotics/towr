/*
 * torque_cost.h
 *
 *  Created on: Jul 1, 2019
 *      Author: vivian
 */

#ifndef TOWR_COSTS_TORQUE_COST_H_
#define TOWR_COSTS_TORQUE_COST_H_

#include <memory>
#include <string>

#include <ifopt/cost_term.h>

#include <towr/variables/nodes_variables_ee_wheels.h>
#include <towr/terrain/height_map.h>

namespace towr {

class TorqueCost : public ifopt::CostTerm {
public:
  using Vector3d = Eigen::Vector3d;

  TorqueCost (const HeightMap::Ptr& terrain, double weight, const uint& ee);
  virtual ~TorqueCost () = default;

  void InitVariableDependedQuantities(const VariablesPtr& x) override;

  double GetCost () const override;

private:
  std::shared_ptr<NodesVariablesEEWheelsForce> ee_force_nodes_;
  std::shared_ptr<NodesVariablesEEWheelsMotion> ee_motion_nodes_;
  HeightMap::Ptr terrain_;
  uint ee_;
  double weight_;

  int GetCol(int node, int dim) const;

  void FillJacobianBlock(std::string var_set, Jacobian&) const override;
};

} /* namespace towr */

#endif /* TOWR_COSTS_TORQUE_COST_H_ */
