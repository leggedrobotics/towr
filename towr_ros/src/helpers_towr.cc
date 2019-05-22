/*
 * helpers_towr.cc
 *
 *  Created on: Mar 28, 2019
 *      Author: vivian
 */

#include "towr_ros/helpers_towr.h"

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
      state.ee_forces_ .at(ee_xpp) = solution.ee_force_.at(ee_towr)->GetPoint(t).p();
    }

    state.t_global_ = t;
    trajectory.push_back(state);
    t += 0.0025;
  }

  return trajectory;
}

std::vector<xpp::RobotStateCartesian> GetDrivingTrajectoryFromSolution (const SplineHolderDrive& solution)
{
  std::vector<xpp::RobotStateCartesian> trajectory;
  double t = 0.0;
  double T = solution.base_linear_->GetTotalTime();

  EulerConverter base_angular(solution.base_angular_);

  while (t<=T+1e-5) {
    int n_ee = solution.ee_wheels_motion_.size();
    xpp::RobotStateCartesian state(n_ee);

    state.base_.lin = ToXpp(solution.base_linear_->GetPoint(t));

    state.base_.ang.q  = base_angular.GetQuaternionBaseToWorld(t);
    state.base_.ang.w  = base_angular.GetAngularVelocityInWorld(t);
    state.base_.ang.wd = base_angular.GetAngularAccelerationInWorld(t);

    for (int ee_towr=0; ee_towr<n_ee; ++ee_towr) {
      int ee_xpp = ToXppEndeffector(n_ee, ee_towr).first;

      state.ee_contact_.at(ee_xpp) = true;
      state.ee_motion_.at(ee_xpp)  = ToXpp(solution.ee_wheels_motion_.at(ee_towr)->GetPoint(t));
      state.ee_forces_ .at(ee_xpp) = solution.ee_wheels_force_.at(ee_towr)->GetPoint(t).p();
    }

    state.t_global_ = t;
    trajectory.push_back(state);
    t += 0.0025;
  }

  return trajectory;
}

void SaveTrajectoryInRosbag (const SplineHolder& solution, const std::string &bag_file)
{
  rosbag::Bag bag;
  bag.open(bag_file, rosbag::bagmode::Write);

  std::vector<xpp::RobotStateCartesian> traj = GetTrajectoryFromSolution(solution);
  xpp_msgs::RobotStateCartesianTrajectory msg = xpp::Convert::ToRos(traj);

  ::ros::Time t0(1e-6);
  bag.write(xpp_msgs::robot_trajectory_desired, t0, msg);

  std::cout << "Successfully created bag  " << bag.getFileName() << std::endl;
  bag.close();

}

void SaveDrivingTrajectoryInRosbag (const SplineHolderDrive& solution, const std::string &bag_file)
{
  rosbag::Bag bag;
  bag.open(bag_file, rosbag::bagmode::Write);

  std::vector<xpp::RobotStateCartesian> traj = GetDrivingTrajectoryFromSolution(solution);
  xpp_msgs::RobotStateCartesianTrajectory msg = xpp::Convert::ToRos(traj);

  ::ros::Time t0(1e-6);
  bag.write(xpp_msgs::robot_trajectory_desired, t0, msg);

  std::cout << "Successfully created bag  " << bag.getFileName() << std::endl;
  bag.close();
}

void
SaveDrivingTrajectoryStatesInRosbag (const SplineHolderDrive& solution, const std::string &bag_file)
{
  std::string topic = xpp_msgs::robot_state_desired;
  std::vector<xpp::RobotStateCartesian> traj = GetDrivingTrajectoryFromSolution(solution);

  rosbag::Bag bag;
  bag.open(bag_file, rosbag::bagmode::Write);

  for (const auto state : traj) {
	auto timestamp = ::ros::Time(state.t_global_ + 1e-6); // t=0.0 throws ROS exception

	xpp_msgs::RobotStateCartesian msg;
	msg = xpp::Convert::ToRos(state);
	bag.write(topic, timestamp, msg);
  }

  std::cout << "Successfully created bag  " << bag.getFileName() << std::endl;
  bag.close();
}

void SaveTerrainNormalsInFile (const SplineHolderDrive& solution, int terrain, const std::string &filename)
{
  double t = 0.0;
  double T = solution.base_linear_->GetTotalTime();

  auto terrain_id = static_cast<HeightMap::TerrainID>(terrain);
  auto terrain_ = HeightMap::MakeTerrain(terrain_id);

  std::ofstream file;
  file.open(filename);
  while (t<=T+1e-5) {
	int n_ee = solution.ee_wheels_motion_.size();
	file << t << " ";

	for (int ee=0; ee<n_ee; ++ee) {
	  Eigen::Vector3d ee_pos = solution.ee_wheels_motion_.at(ee)->GetPoint(t).p();
	  Eigen::Vector3d n = terrain_->GetNormalizedBasis(HeightMap::Normal,ee_pos(X),ee_pos(Y));
	  file << n.transpose() << " ";
	}

	file << "\n";
	t += 0.0025;
  }

  file.close();
  std::cout << "Successfully created file " + filename << std::endl;
}

void SaveTerrainNormalsInFile (const SplineHolder& solution, int terrain, const std::string &filename)
{
  double t = 0.0;
  double T = solution.base_linear_->GetTotalTime();

  auto terrain_id = static_cast<HeightMap::TerrainID>(terrain);
  auto terrain_ = HeightMap::MakeTerrain(terrain_id);

  std::ofstream file;
  file.open(filename);
  while (t<=T+1e-5) {
	int n_ee = solution.ee_motion_.size();
	file << t << " ";

	for (int ee=0; ee<n_ee; ++ee) {
	  Eigen::Vector3d ee_pos = solution.ee_motion_.at(ee)->GetPoint(t).p();
	  Eigen::Vector3d n = terrain_->GetNormalizedBasis(HeightMap::Normal,ee_pos(X),ee_pos(Y));
	  file << n.transpose() << " ";
	}

	file << "\n";
	t += 0.0025;
  }

  file.close();
  std::cout << "Successfully created file " + filename << std::endl;
}

void
ExtractGeometryMessagesFromTrajectoryBag (const std::string bag_file)
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
	std::string path = ros::package::getPath("towr_test") + "/bags/anymal_wheels_matlab.bag";
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

} // namespace towr


