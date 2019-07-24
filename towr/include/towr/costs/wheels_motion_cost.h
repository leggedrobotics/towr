/*
 * wheels_motion_cost.h
 *
 *  Created on: Jul 20, 2019
 *      Author: vivian
 */

#ifndef TOWR_TOWR_INCLUDE_TOWR_COSTS_WHEELS_MOTION_COST_H_
#define TOWR_TOWR_INCLUDE_TOWR_COSTS_WHEELS_MOTION_COST_H_

#include <memory>
#include <string>

#include <ifopt/cost_term.h>

#include <towr/variables/nodes_variables_ee_wheels.h>
#include <towr/terrain/height_map.h>

namespace towr {

class WheelsMotionCost : public ifopt::CostTerm {
public:
  using Vector3d = Eigen::Vector3d;

  WheelsMotionCost (const double dt, double weight);
  virtual ~WheelsMotionCost () = default;

  void InitVariableDependedQuantities(const VariablesPtr& x) override;

  double GetCost () const override;

private:
  std::shared_ptr<NodesVariablesEEWheelsMotion> LF_motion_nodes_;
  std::shared_ptr<NodesVariablesEEWheelsMotion> RF_motion_nodes_;
  double weight_;
  double dt_;

  int GetCol(int node, int dim) const;

  void FillJacobianBlock(std::string var_set, Jacobian&) const override;
};

} /* namespace towr */


#endif /* TOWR_TOWR_INCLUDE_TOWR_COSTS_WHEELS_MOTION_COST_H_ */
