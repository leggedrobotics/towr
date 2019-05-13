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

#include "towr/variables/euler_converter.h"
#include "towr/variables/spline_holder_drive.h"
#include "towr/terrain/height_map.h"

#include <towr_ros/towr_xpp_ee_map.h>

namespace towr {

std::vector<xpp::RobotStateCartesian> GetTrajectoryFromSolution (const SplineHolder& solution);
std::vector<xpp::RobotStateCartesian> GetDrivingTrajectoryFromSolution (const SplineHolderDrive& solution);

void SaveTrajectoryInRosbag (const SplineHolder& solution, const std::string &bag_file);
void SaveDrivingTrajectoryInRosbag (const SplineHolderDrive& solution, const std::string &bag_file);
void SaveDrivingTrajectoryStatesInRosbag (const SplineHolderDrive& solution, const std::string &bag_file);

void SaveTerrainNormalsInFile (const SplineHolderDrive& solution, int terrain_id, const std::string &bag_file);

void ExtractGeometryMessagesFromTrajectoryBag (const std::string bag_file);

}

