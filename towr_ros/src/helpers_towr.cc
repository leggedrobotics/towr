/*
 * helpers_towr.cc
 *
 *  Created on: Mar 28, 2019
 *      Author: vivian
 */

#include "towr_ros/helpers_towr.h"

#include <xpp_msgs/TerrainInfo.h>

#include <std_msgs/UInt8.h>
#include <std_msgs/Float64MultiArray.h>

namespace towr {

std::vector<xpp::RobotStateCartesian> GetTrajectoryFromSolution (const SplineHolder& solution)
{
  std::vector<xpp::RobotStateCartesian> trajectory;
  double t = 0.0;
  double T = solution.base_linear_->GetTotalTime();

  EulerConverter base_angular(solution.base_angular_);

  while (t<=T+1e-5) {
    int n_ee = solution.ee_motion_.size();
    xpp::RobotStateCartesian state(n_ee);

    state.base_.lin = ToXpp(solution.base_linear_->GetPoint(t));

    state.base_.ang.q  = base_angular.GetQuaternionBaseToWorld(t);
    state.base_.ang.w  = base_angular.GetAngularVelocityInWorld(t);
    state.base_.ang.wd = base_angular.GetAngularAccelerationInWorld(t);

    for (int ee_towr=0; ee_towr<n_ee; ++ee_towr) {
      int ee_xpp = ToXppEndeffector(n_ee, ee_towr).first;

      state.ee_contact_.at(ee_xpp) = solution.phase_durations_.at(ee_towr)->IsContactPhase(t);
      state.ee_motion_.at(ee_xpp)  = ToXpp(solution.ee_motion_.at(ee_towr)->GetPoint(t));
      state.ee_forces_.at(ee_xpp) = solution.ee_force_.at(ee_towr)->GetPoint(t).p();
    }

    state.t_global_ = t;
    trajectory.push_back(state);
    t += 0.0025;
  }

  return trajectory;
}

void getDataFromBag (std::string bagname)
{
  std::vector<xpp_msgs::RobotStateCartesian> desiredStates_;
  std::map<double,std::vector<Eigen::Vector3d>> terrainNormalMap_;

  rosbag::Bag bag;
  bag.open(bagname, rosbag::bagmode::Read);

  std::vector<std::string> topics;
  topics.push_back(xpp_msgs::robot_state_desired);
  topics.push_back(xpp_msgs::terrain_info);

  rosbag::View view(bag, rosbag::TopicQuery(topics));

  double t = 0.0;
  for (rosbag::MessageInstance const m: view)
  {
	t = m.getTime().toSec();
	xpp_msgs::RobotStateCartesianPtr traj_msg = m.instantiate<xpp_msgs::RobotStateCartesian>();
	if (traj_msg != NULL) {
	  desiredStates_.push_back(*traj_msg);
	}

    std::cout << std::fixed;
    std::cout << std::setprecision(5);
	std::cout << t << std::endl;

	xpp_msgs::TerrainInfoPtr terrain_msg = m.instantiate<xpp_msgs::TerrainInfo>();
	if (terrain_msg != NULL) {
	  auto n = xpp::Convert::ToXpp<geometry_msgs::Vector3>(terrain_msg->surface_normals).ToImpl();
	  std::vector<Eigen::Vector3d> normals({n.begin(), n.end()});
	  terrainNormalMap_.insert( std::pair<double,std::vector<Eigen::Vector3d>>(t,normals) );
	}
  }
  bag.close();

  std::cout << "States: " << desiredStates_.size() << std::endl;
  std::cout << "Normals: " << terrainNormalMap_.size() << std::endl;
  for (auto const& x : terrainNormalMap_)
  {
      std::cout << x.first  // string (key)
                << ':'
                << x.second.at(0).transpose() << ", "// string's value
				<< x.second.at(1).transpose() << ", "
				<< x.second.at(2).transpose() << ", "
				<< x.second.at(3).transpose()
                << std::endl ;
  }
}

void
ExtractGeometryMessagesFromTrajectoryBag (const std::string bag_file, const SplineHolder& solution)
{
	rosbag::Bag bag_r;
	bag_r.open(bag_file, rosbag::bagmode::Read);
//	std::cout << "Reading from bag " << bag_r.getFileName() << std::endl;

	// select which iterations (message topics) to be included in bag file
	std::string topic = "/xpp/state_des";
	rosbag::View view(bag_r, rosbag::TopicQuery(topic));
	if (view.size() == 0) {
	  std::cerr << "Error: Topic " << topic << " doesn't exist\n";
	}

	// write the message with modified timestamp into new bag file
	rosbag::Bag bag_w;
	std::string path = ros::package::getPath("towr_ros") + "/bags/anymal_wheels_matlab.bag";
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

		std_msgs::UInt8 c;
		c.data = state_msg->ee_contact.at(i);
		bag_w.write("foot_contact_"+std::to_string(i), t, c);
	  }
	}

	int n_ee = solution.ee_motion_.size();
	auto timestamp = ::ros::Time(1e-6);

	for (int ee=0; ee<n_ee; ++ee) {
	  std::vector<double> motion_poly_durations = solution.ee_motion_.at(ee)->GetPolyDurations();
	  std::vector<double> force_poly_durations = solution.ee_force_.at(ee)->GetPolyDurations();

	  std_msgs::Float64MultiArray poly_motion;
	  for (auto p : motion_poly_durations)
		  poly_motion.data.push_back(p);

	  std_msgs::Float64MultiArray poly_force;
	  for (auto p : force_poly_durations)
		  poly_force.data.push_back(p);

	  bag_w.write("poly_dur_pos_"+std::to_string(ee), timestamp, poly_motion);
	  bag_w.write("poly_dur_force_"+std::to_string(ee), timestamp, poly_force);
	}

	bag_r.close();
	std::cout << "Successfully created bag  " << bag_w.getFileName() << std::endl;
	bag_w.close();
}

void
AddPolyDurationsToBag (const std::string bag_file, const SplineHolder& solution) {
	rosbag::Bag bag;
	bag.open(bag_file, rosbag::bagmode::Write);

	int n_ee = solution.ee_motion_.size();
	auto timestamp = ::ros::Time(1e-6);

	for (int ee=0; ee<n_ee; ++ee) {
	  std::vector<double> motion_poly_durations = solution.ee_motion_.at(ee)->GetPolyDurations();
	  std::vector<double> force_poly_durations = solution.ee_force_.at(ee)->GetPolyDurations();

	  std_msgs::Float64MultiArray poly_motion;
	  for (auto p : motion_poly_durations)
		  poly_motion.data.push_back(p);

	  std_msgs::Float64MultiArray poly_force;
	  for (auto p : force_poly_durations)
		  poly_force.data.push_back(p);

	  bag.write("poly_dur_pos_"+std::to_string(ee), timestamp, poly_motion);
	  bag.write("poly_dur_force_"+std::to_string(ee), timestamp, poly_force);
	}

}

} // namespace towr


