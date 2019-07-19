/*
 * send_trajectory_from_bag.cc
 *
 *  Created on: Jul 17, 2019
 *      Author: vivian
 */

#include <iostream>
#include <string>
#include <vector>

#include <ros/ros.h>
#include <std_msgs/Bool.h>

#include <rosbag/bag.h>
#include <rosbag/view.h>
#include <rosbag/message_instance.h>

#include <anymal_wheels_ctrl_track_msgs/TerrainMap.h>
#include <xpp_msgs/RobotStateCartesian.h>
#include <xpp_msgs/RobotStateCartesianTrajectory.h>
#include <xpp_msgs/topic_names.h>

int main(int argc, char *argv[])
{
  if (argc==1) {
	std::cerr << "Error: Please enter path to bag file\n";
	return 0;
  }

  std::string bag_file = argv[1];

  ros::init(argc, argv, "towr_trajectory_sender");

  ros::NodeHandle n;

  ros::Publisher trajectory_pub;
  ros::Publisher terrain_map_pub;
  ros::Publisher new_trajectory_pub;

  trajectory_pub = n.advertise<xpp_msgs::RobotStateCartesianTrajectory>(xpp_msgs::robot_trajectory_desired, 1);
  new_trajectory_pub = n.advertise<std_msgs::Bool>("/xpp/new_trajectory", 1);
  terrain_map_pub = n.advertise<anymal_wheels_ctrl_track_msgs::TerrainMap>("/towr/terrain_map", 1);

  // Get trajectory and terrain map from rosbag ----------------------
  rosbag::Bag bag_r;
  bag_r.open(bag_file, rosbag::bagmode::Read);
  std::cout << "Reading from bag " + bag_r.getFileName() << std::endl;

  std::vector<std::string> topics;
  topics.push_back(xpp_msgs::robot_state_desired);
  topics.push_back(xpp_msgs::terrain_info);

  rosbag::View view(bag_r, rosbag::TopicQuery(topics));
//  std::cout << "Size: " << view.size() << std::endl;
  if (view.size() == 0) {
    std::cerr << "Error: At least one of the topics doesn't exist\n";
    return 0;
  }

  xpp_msgs::RobotStateCartesianTrajectory traj_ros;
  anymal_wheels_ctrl_track_msgs::TerrainMap map;

  for (rosbag::MessageInstance const m: view)
  {
	double t = m.getTime().toSec();
//	std::cout << t << std::endl;
	xpp_msgs::RobotStateCartesianPtr state_ros = m.instantiate<xpp_msgs::RobotStateCartesian>();
	if (state_ros != NULL) {
	  traj_ros.points.push_back(*state_ros);
	}

	xpp_msgs::TerrainInfoPtr terrain_msg = m.instantiate<xpp_msgs::TerrainInfo>();
	if (terrain_msg != NULL) {
	  map.terrain.push_back(*terrain_msg);
	}
  }

  bag_r.close();
  // -------------------------------------------------------------------

  std_msgs::Bool new_traj;
  new_traj.data = true;
  new_trajectory_pub.publish(new_traj);

  ros::Rate loop_rate(10);

  int i = 0;
  while (ros::ok() && i < 2)
  {
//	std::cout << i << std::endl;
	trajectory_pub.publish(traj_ros);
	terrain_map_pub.publish(map);

	ros::spinOnce();
	loop_rate.sleep();
	i++;
  }

  ros::spin();

  return 1;

}


