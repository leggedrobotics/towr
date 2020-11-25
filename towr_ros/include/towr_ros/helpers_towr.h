/*
 * helpers_towr.h
 *
 *  Created on: Mar 28, 2019
 *      Author: vivian
 */
#include <iostream>

#include <ros/ros.h>
#include <rosbag/bag.h>
#include <rosbag/view.h>
#include <rosbag/message_instance.h>
#include <boost/foreach.hpp>

#include <xpp_states/robot_state_cartesian.h>
#include <xpp_states/endeffector_mappings.h>
#include <xpp_states/convert.h>

#include <xpp_msgs/topic_names.h>
#include <xpp_msgs/RobotStateCartesianTrajectory.h>

#include "towr/variables/spline_holder.h"
#include "towr/variables/euler_converter.h"
#include "towr/terrain/height_map.h"

#include <towr_ros/towr_xpp_ee_map.h>

namespace towr {

std::vector<xpp::RobotStateCartesian> GetTrajectoryFromSolution (const SplineHolder& solution);

// method to check if the bag was created correctly
void getDataFromBag (std::string bagname);

void ExtractGeometryMessagesFromTrajectoryBag (const std::string bag_file, const SplineHolder& solution);

}

