/*
 * tower_drive_ros_interface.cc
 *
 *  Created on: Jun 11, 2019
 *      Author: vivian
 */

#include <towr_drive/towr_drive_ros_interface.h>

#include <std_msgs/Float64.h>
#include <std_msgs/Int32.h>
#include <std_msgs/Bool.h>

#include <xpp_states/convert.h>
#include <xpp_msgs/topic_names.h>
#include <xpp_msgs/TerrainInfo.h>

#include <towr/terrain/height_map.h>
#include <towr/variables/euler_converter.h>
#include <towr_ros/topic_names.h>
#include <towr_ros/towr_xpp_ee_map.h>
#include <towr_ros/helpers_towr.h>

namespace towr {


TowrDriveRosInterface::TowrDriveRosInterface () : current_state_robot_(4)
{
  ::ros::NodeHandle n;

  trajectory_pub_ = n.advertise<xpp_msgs::RobotStateCartesianTrajectory>(xpp_msgs::robot_trajectory_desired, 1);

  new_trajectory_pub_ = n.advertise<std_msgs::Bool>("/xpp/new_trajectory", 1);

  terrain_map_pub_ = n.advertise<anymal_wheels_ctrl_track_msgs::TerrainMap>("/towr/terrain_map", 1);

  towr_command_pub_ = n.advertise<towr_ros::TowrCommand>(towr_msgs::user_command, 1);

  terrain_ref_height_pub_ = n.advertise<std_msgs::Float64>("/towr/terrain_ref_height", 1);

  plan_service_ = n.advertiseService("towr_drive/plan_motion", &TowrDriveRosInterface::planServiceCallback, this);
  replay_service_ = n.advertiseService("towr_drive/replay_motion", &TowrDriveRosInterface::replayServiceCallback, this);

  initial_state_pub_  = n.advertise<xpp_msgs::RobotStateCartesian>(xpp_msgs::robot_state_desired, 1);

  robot_parameters_pub_  = n.advertise<xpp_msgs::RobotParameters>(xpp_msgs::robot_parameters, 1);

  current_state_sub_ = n.subscribe(xpp_msgs::robot_state_current, 1, &TowrDriveRosInterface::currentStateCallback, this);

  solver_ = std::make_shared<ifopt::IpoptSolver>();

  visualization_dt_ = 0.01;
  trajectory_dt_ = 0.0025;
}

bool
TowrDriveRosInterface::planServiceCallback(std_srvs::Trigger::Request  &req,
		  	  	  	       	   	   	   	   std_srvs::Trigger::Response &res)
{
  // robot model
  formulation_.model_ = RobotModel(RobotModel::AnymalWheels);
  auto robot_params_msg = BuildRobotParametersMsg(formulation_.model_);
  robot_parameters_pub_.publish(robot_params_msg);

  // set initial and goal position for TOWR + towr command message
  TowrCommandMsg msg = BuildTowrCommandMsg(current_state_robot_);
  towr_command_pub_.publish(msg);

  int n_ee = formulation_.model_.kinematic_model_->GetNumberOfEndeffectors();
  formulation_.params_drive_ = GetTowrDriveParameters(n_ee, msg);
  formulation_.final_base_ = GetGoalState(msg);

  // solver parameters
  SetIpoptParameters(msg);

  // visualization
  PublishInitialState();

  // Defaults to /home/user/.ros/
  std::string bag_file = "towr_trajectory.bag";
  if (msg.optimize || msg.play_initialization) {
    nlp_ = ifopt::Problem();
    for (auto c : formulation_.GetVariableSets(solution))
      nlp_.AddVariableSet(c);
    for (auto c : formulation_.GetConstraints(solution))
      nlp_.AddConstraintSet(c);
    for (auto c : formulation_.GetCosts())
      nlp_.AddCostSet(c);

    solver_->Solve(nlp_);
    SaveOptimizationAsRosbag(bag_file, robot_params_msg, msg, false);
  }

  nlp_.PrintCurrent();

  // playback using terminal commands
  if (msg.replay_trajectory || msg.play_initialization || msg.optimize) {
    int success = system(("rosbag play --topics "
        + xpp_msgs::robot_state_desired + " "
        + xpp_msgs::terrain_info
        + " -r " + std::to_string(msg.replay_speed)
        + " --quiet " + bag_file).c_str());
  }

  if (msg.plot_trajectory) {
    int success = system(("killall rqt_bag; rqt_bag " + bag_file + "&").c_str());
  }

  // publish trajectory and terrain map for the controller
  XppVec trajectory = GetTrajectory(trajectory_dt_);
  std_msgs::Bool new_traj;
  new_traj.data = true;
  xpp_msgs::RobotStateCartesianTrajectory xpp_msg = xpp::Convert::ToRos(trajectory);
  new_trajectory_pub_.publish(new_traj);
  trajectory_pub_.publish(xpp_msg);
  anymal_wheels_ctrl_track_msgs::TerrainMap terrain_msg = GetTerrainMap(trajectory);
  terrain_map_pub_.publish(terrain_msg);

  // save bags for controller and matlab
  ExtractGeometryMessagesFromTrajectoryBag(bag_file);
  bag_file = ros::package::getPath("anymal_wheels_ctrl_track_ros") + "/data/anymal_wheels_traj.bag";
  SaveDrivingMotionTerrainInRosbag (solution, msg.terrain, bag_file);

  res.success = true;
  res.message = "optimization done";

  return true;
}

bool
TowrDriveRosInterface::replayServiceCallback(std_srvs::Trigger::Request  &req,
		  	  	  	       	   	   	   	     std_srvs::Trigger::Response &res)
{
  std::string bag_file = "towr_trajectory.bag";
  int success = system(("rosbag play --topics "
		  	  	  	  + xpp_msgs::robot_state_desired + " "
					  + xpp_msgs::terrain_info
					  + " --quiet " + bag_file).c_str());

  res.success = true;
  res.message = "replay done";

  return true;
}

void
TowrDriveRosInterface::currentStateCallback (const xpp_msgs::RobotStateCartesianConstPtr& msg)
{
  xpp::RobotStateCartesian xpp = xpp::Convert::ToXpp(*msg);

  current_state_robot_ = xpp;
}

BaseState
TowrDriveRosInterface::GetGoalState(const TowrCommandMsg& msg) const
{
  BaseState goal;
  goal.lin.at(kPos) = xpp::Convert::ToXpp(msg.goal_lin.pos);
  goal.lin.at(kVel) = xpp::Convert::ToXpp(msg.goal_lin.vel);
  goal.ang.at(kPos) = xpp::Convert::ToXpp(msg.goal_ang.pos);
  goal.ang.at(kVel) = xpp::Convert::ToXpp(msg.goal_ang.vel);

  return goal;
}

void
TowrDriveRosInterface::PublishInitialState()
{
  int n_ee = formulation_.initial_ee_W_.size();
  xpp::RobotStateCartesian xpp(n_ee);
  xpp.base_.lin.p_ = formulation_.initial_base_.lin.p();
  xpp.base_.ang.q  = EulerConverter::GetQuaternionBaseToWorld(formulation_.initial_base_.ang.p());

  for (int ee_towr=0; ee_towr<n_ee; ++ee_towr) {
    int ee_xpp = ToXppEndeffector(n_ee, ee_towr).first;
    xpp.ee_contact_.at(ee_xpp)   = true;
    xpp.ee_motion_.at(ee_xpp).p_ = formulation_.initial_ee_W_.at(ee_towr);
    xpp.ee_forces_.at(ee_xpp).setZero(); // zero for visualization
  }

  initial_state_pub_.publish(xpp::Convert::ToRos(xpp));
}

std::vector<TowrDriveRosInterface::XppVec>
TowrDriveRosInterface::GetIntermediateSolutions ()
{
  std::vector<XppVec> trajectories;

  for (int iter=0; iter<nlp_.GetIterationCount(); ++iter) {
    nlp_.SetOptVariables(iter);
    trajectories.push_back(GetTrajectory(visualization_dt_));
  }

  return trajectories;
}

TowrDriveRosInterface::XppVec
TowrDriveRosInterface::GetTrajectory (const double dt) const
{
  XppVec trajectory;
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
    t += dt;
  }

  return trajectory;
}

xpp_msgs::RobotParameters
TowrDriveRosInterface::BuildRobotParametersMsg(const RobotModel& model) const
{
  xpp_msgs::RobotParameters params_msg;
  auto max_dev_xyz = model.kinematic_model_->GetMaximumDeviationFromNominal();
  params_msg.ee_max_dev = xpp::Convert::ToRos<geometry_msgs::Vector3>(max_dev_xyz);

  auto nominal_B = model.kinematic_model_->GetNominalStanceInBase();
  int n_ee = nominal_B.size();
  for (int ee_towr=0; ee_towr<n_ee; ++ee_towr) {
    Vector3d pos = nominal_B.at(ee_towr);
    params_msg.nominal_ee_pos.push_back(xpp::Convert::ToRos<geometry_msgs::Point>(pos));
    params_msg.ee_names.push_back(ToXppEndeffector(n_ee, ee_towr).second);
  }

  params_msg.base_mass = model.dynamic_model_->m();

  return params_msg;
}

void
TowrDriveRosInterface::SaveOptimizationAsRosbag (const std::string& bag_name,
                                   const xpp_msgs::RobotParameters& robot_params,
                                   const TowrCommandMsg user_command_msg,
                                   bool include_iterations)
{
  rosbag::Bag bag;
  bag.open(bag_name, rosbag::bagmode::Write);
  ::ros::Time t0(1e-6); // t=0.0 throws ROS exception

  // save the a-priori fixed optimization variables
  bag.write(xpp_msgs::robot_parameters, t0, robot_params);
  bag.write(towr_msgs::user_command+"_saved", t0, user_command_msg);

  // save the trajectory of each iteration
  if (include_iterations) {
    auto trajectories = GetIntermediateSolutions();
    int n_iterations = trajectories.size();
    for (int i=0; i<n_iterations; ++i)
      SaveTrajectoryInRosbag(bag, trajectories.at(i), towr_msgs::nlp_iterations_name + std::to_string(i));

    // save number of iterations the optimizer took
    std_msgs::Int32 m;
    m.data = n_iterations;
    bag.write(towr_msgs::nlp_iterations_count, t0, m);
  }

  // save the final trajectory
  auto final_trajectory = GetTrajectory(visualization_dt_);
  SaveTrajectoryInRosbag(bag, final_trajectory, xpp_msgs::robot_state_desired);

  bag.close();
}

void
TowrDriveRosInterface::SaveTrajectoryInRosbag (rosbag::Bag& bag,
                                 const XppVec& traj,
                                 const std::string& topic) const
{
  for (const auto state : traj) {
    auto timestamp = ::ros::Time(state.t_global_ + 1e-6); // t=0.0 throws ROS exception

    xpp_msgs::RobotStateCartesian msg;
    msg = xpp::Convert::ToRos(state);
    bag.write(topic, timestamp, msg);

    xpp_msgs::TerrainInfo terrain_msg;
    for (auto ee : state.ee_motion_.ToImpl()) {
      Vector3d n = formulation_.terrain_->GetNormalizedBasis(HeightMap::Normal, ee.p_.x(), ee.p_.y());
      terrain_msg.surface_normals.push_back(xpp::Convert::ToRos<geometry_msgs::Vector3>(n));
      terrain_msg.friction_coeff = formulation_.terrain_->GetFrictionCoeff();
    }

    bag.write(xpp_msgs::terrain_info, timestamp, terrain_msg);
  }
}

anymal_wheels_ctrl_track_msgs::TerrainMap
TowrDriveRosInterface::GetTerrainMap (const XppVec& traj) const
{
  anymal_wheels_ctrl_track_msgs::TerrainMap map;

  for (const auto state : traj) {
    auto timestamp = ::ros::Time(state.t_global_ + 1e-6); // t=0.0 throws ROS exception

    xpp_msgs::TerrainInfo terrain_msg;
    for (auto ee : state.ee_motion_.ToImpl()) {
      Vector3d n = formulation_.terrain_->GetNormalizedBasis(HeightMap::Normal, ee.p_.x(), ee.p_.y());
      terrain_msg.surface_normals.push_back(xpp::Convert::ToRos<geometry_msgs::Vector3>(n));
      terrain_msg.friction_coeff = formulation_.terrain_->GetFrictionCoeff();
    }

    map.terrain.push_back(terrain_msg);

  }

  return map;
}

} /* namespace towr */



