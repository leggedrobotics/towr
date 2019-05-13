/*
 * nodes_variables_ee_wheels.h
 *
 *  Created on: Mar 29, 2019
 *      Author: vivian
 */

#ifndef TOWR_INCLUDE_TOWR_VARIABLES_NODES_VARIABLES_EE_WHEELS_H_
#define TOWR_INCLUDE_TOWR_VARIABLES_NODES_VARIABLES_EE_WHEELS_H_

#include "nodes_variables_all.h"
#include <towr/terrain/height_map.h>
#include <towr/models/dynamic_model.h>

namespace towr {

class NodesVariablesEEWheelsMotion : public NodesVariablesAll {
public:
  /**
   * @param n_nodes  Number of nodes to construct the spline.
   * @param n_dim    Number of dimensions of each node.
   * @param variable_id  Name of this variables set in the optimization.
   */
  NodesVariablesEEWheelsMotion (int n_nodes, int n_dim, std::string variable_id, const HeightMap::Ptr& terrain);
  virtual ~NodesVariablesEEWheelsMotion () = default;

  // Initialize the variables as a linear interpolation between the initial and final value.
  // The difference between this method and the base method is that enforces contact with the terrain.
  void SetByLinearInterpolation(const VectorXd& initial_val,
                                const VectorXd& final_val,
                                double t_total);

private:
  HeightMap::Ptr terrain_;

};

class NodesVariablesEEWheelsForce : public NodesVariablesAll {
public:
  /**
   * @param n_nodes  Number of nodes to construct the spline.
   * @param n_dim    Number of dimensions of each node.
   * @param variable_id  Name of this variables set in the optimization.
   */
  NodesVariablesEEWheelsForce (int n_nodes, int n_dim, std::string variable_id);
  virtual ~NodesVariablesEEWheelsForce () = default;

  void SetAllAtConstantValue (const VectorXd& force);

};

}  // namespace towr

#endif /* TOWR_INCLUDE_TOWR_VARIABLES_NODES_VARIABLES_EE_WHEELS_H_ */
