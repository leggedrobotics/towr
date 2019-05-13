/*
 * spline_holder_drive.cc
 *
 *  Created on: Mar 28, 2019
 *      Author: vivian
 */

#include <towr/variables/spline_holder_drive.h>

namespace towr{

SplineHolderDrive::SplineHolderDrive (NodesVariables::Ptr base_lin_nodes,
									  NodesVariables::Ptr base_ang_nodes,
									  const std::vector<double>& base_poly_durations,
									  std::vector<NodesVariables::Ptr> ee_wheels_motion_nodes,
									  std::vector<NodesVariables::Ptr> ee_wheels_force_nodes,
									  const std::vector<double>& ee_wheels_poly_durations)
{
  base_linear_  = std::make_shared<NodeSpline>(base_lin_nodes.get(), base_poly_durations);
  base_angular_ = std::make_shared<NodeSpline>(base_ang_nodes.get(), base_poly_durations);

  for (uint ee=0; ee<ee_wheels_motion_nodes.size(); ++ee) {
      ee_wheels_motion_.push_back(std::make_shared<NodeSpline>(ee_wheels_motion_nodes.at(ee).get(), ee_wheels_poly_durations));
      ee_wheels_force_.push_back (std::make_shared<NodeSpline>(ee_wheels_force_nodes.at(ee).get(), ee_wheels_poly_durations));
  }
}

} /* namespace towr */



