/*
 * spline_holder_drive.h
 *
 *  Created on: Mar 28, 2019
 *      Author: vivian
 */

#ifndef TOWR_TOWR_INCLUDE_TOWR_VARIABLES_SPLINE_HOLDER_DRIVE_H_
#define TOWR_TOWR_INCLUDE_TOWR_VARIABLES_SPLINE_HOLDER_DRIVE_H_

#include "node_spline.h"
#include "spline_holder.h"
#include "nodes_variables.h"

namespace towr {

/**
 * SplineHolder for driving motions only. The nodes representing the wheels motion and forces are no longer phase based.
 * Instead, they are optimized in fixed intervals similar to the base variables.
 */
struct SplineHolderDrive
{
  /**
   * @brief Fully construct all splines.
   * @param base_lin_nodes  		 The nodes describing the base linear motion.
   * @param base_ang_nodes  		 The nodes describing the base angular motion.
   * @param base_poly_durations 	 The durations of each base polynomial.
   * @param ee_wheels_motion_nodes 	 The nodes describing the wheels motions.
   * @param ee_wheels_force_nodes	 The nodes describing the wheels forces.
   * @param ee_wheels_poly_durations The durations of each wheel's polynomial.
   */
  SplineHolderDrive (NodesVariables::Ptr base_lin_nodes,
                	 NodesVariables::Ptr base_ang_nodes,
					 const std::vector<double>& base_poly_durations,
					 std::vector<NodesVariables::Ptr> ee_wheels_motion_nodes,
					 std::vector<NodesVariables::Ptr> ee_wheels_force_nodes,
					 const std::vector<double>& ee_wheels_poly_durations);

  SplineHolderDrive () = default;

  NodeSpline::Ptr base_linear_;
  NodeSpline::Ptr base_angular_;
  std::vector<NodeSpline::Ptr> ee_wheels_motion_;
  std::vector<NodeSpline::Ptr> ee_wheels_force_;
};

} /* namespace towr */

#endif /* TOWR_TOWR_INCLUDE_TOWR_VARIABLES_SPLINE_HOLDER_DRIVE_H_ */
