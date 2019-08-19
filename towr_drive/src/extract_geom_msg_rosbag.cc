/*
 * extract_geom_msg_rosbag.cc
 *
 *  Created on: Aug 11, 2019
 *      Author: vivian
 */

#include <iostream>
#include <string>
#include <vector>

#include <rosbag/bag.h>
#include <rosbag/view.h>
#include <rosbag/message_instance.h>
#include <boost/foreach.hpp>

#include <xpp_msgs/RobotStateCartesian.h>

int main(int argc, char *argv[])
{
	if (argc==1) {
		std::cerr << "Error: Please enter path to bag file\n";
		return 0;
	}

	std::string bag_file = argv[1];

	rosbag::Bag bag_r;
	bag_r.open(bag_file, rosbag::bagmode::Read);
	std::cout << "Reading from bag " << bag_r.getFileName() << std::endl;

	// select which iterations (message topics) to be included in bag file
	std::string topic = "/xpp/state_des";
	rosbag::View view(bag_r, rosbag::TopicQuery(topic));
	if (view.size() == 0) {
	  std::cerr << "Error: Topic " << topic << " doesn't exist\n";
	}

	// write the message with modified timestamp into new bag file
	rosbag::Bag bag_w;
	std::string path = ros::package::getPath("towr_drive") + "/bags/anymal_wheels_matlab_extracted.bag";
	bag_w.open(path, rosbag::bagmode::Write);

	BOOST_FOREACH(rosbag::MessageInstance const m, view)
	{
	  ros::Time t = m.getTime();
	  auto state_msg = m.instantiate<xpp_msgs::RobotStateCartesian>();
	  bag_w.write("base_pose", t, state_msg->base.pose);
	  bag_w.write("base_vel_lin", t, state_msg->base.twist.linear);
	  bag_w.write("base_vel_ang", t, state_msg->base.twist.angular);
	  bag_w.write("base_acc_lin", t, state_msg->base.accel.linear);
	  bag_w.write("base_acc_ang", t, state_msg->base.accel.angular);

	  int n_feet = state_msg->ee_motion.size();

	  for (int i=0; i<n_feet; ++i) {
		bag_w.write("foot_pos_"+std::to_string(i), t, state_msg->ee_motion.at(i).pos);
		bag_w.write("foot_vel_"+std::to_string(i), t, state_msg->ee_motion.at(i).vel);
		bag_w.write("foot_acc_"+std::to_string(i), t, state_msg->ee_motion.at(i).acc);
		bag_w.write("foot_force_"+std::to_string(i), t, state_msg->ee_forces.at(i));
	  }
	}

	bag_r.close();
	std::cout << "Successfully created bag  " << bag_w.getFileName() << std::endl;
	bag_w.close();
}
