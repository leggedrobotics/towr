/*
 * towr_drive_ros_interface.h
 *
 *  Created on: Jun 11, 2019
 *      Author: vivian
 */

#ifndef TOWR_DRIVE_INCLUDE_TOWR_DRIVE_TOWR_DRIVE_ROS_INTERFACE_H_
#define TOWR_DRIVE_INCLUDE_TOWR_DRIVE_TOWR_DRIVE_ROS_INTERFACE_H_

#include <string>

#include <ros/ros.h>
#include <rosbag/bag.h>

#include <std_srvs/Empty.h>
#include <std_srvs/Trigger.h>

#include <xpp_states/robot_state_cartesian.h>
#include <xpp_msgs/RobotStateCartesian.h>
#include <xpp_msgs/RobotParameters.h>
#include <towr_ros/TowrCommand.h>
#include <terrain_msgs/TerrainMap.h>

#include <towr/nlp_formulation.h>
#include <towr/nlp_formulation_drive.h>
#include <ifopt/ipopt_solver.h>

namespace towr {

class TowrDriveRosInterface {
public:
  using XppVec         = std::vector<xpp::RobotStateCartesian>;
  using TowrCommandMsg = towr_ros::TowrCommand;
  using Vector3d       = Eigen::Vector3d;

protected:
  TowrDriveRosInterface ();
  virtual ~TowrDriveRosInterface () = default;

  virtual ParametersDrive GetTowrDriveParameters(int n_ee, const TowrCommandMsg& msg) const = 0;

  virtual void SetIpoptParameters(const TowrCommandMsg& msg) = 0;

  virtual TowrCommandMsg BuildTowrCommandMsg (const xpp::RobotStateCartesian& init_state) = 0;

  NlpFormulationDrive formulation_;         ///< the default formulation, can be adapted

  ifopt::IpoptSolver::Ptr solver_; ///< NLP solver, could also use SNOPT.

private:
  SplineHolderDrive solution; ///< the solution splines linked to the opt-variables.
  ifopt::Problem nlp_;        ///< the actual nonlinear program to be solved.
  double visualization_dt_;   ///< duration between two rviz visualization states.
  double trajectory_dt_;      ///< interval between the states in the trajectory.

  xpp::RobotStateCartesian current_state_robot_;

  ::ros::Publisher trajectory_pub_;
  ::ros::Publisher terrain_map_pub_;
  ::ros::Publisher new_trajectory_pub_;
  ::ros::Publisher initial_state_pub_;
  ::ros::Publisher robot_parameters_pub_;
  ::ros::ServiceServer plan_service_;
  ::ros::ServiceServer replay_service_;
  ::ros::Publisher  towr_command_pub_;
  ::ros::Subscriber current_state_sub_;

  terrain_msgs::TerrainMap GetTerrainMap (const XppVec& traj) const;
  bool planServiceCallback(std_srvs::Trigger::Request  &req,
	   	   	   	   	   	   std_srvs::Trigger::Response &res);
  bool replayServiceCallback(std_srvs::Trigger::Request  &req,
	   	   	   	   	   	     std_srvs::Trigger::Response &res);
  void currentStateCallback (const xpp_msgs::RobotStateCartesianConstPtr& msg);
  XppVec GetTrajectory(const double dt) const;
  virtual BaseState GetGoalState(const TowrCommandMsg& msg) const;
  void PublishInitialState();
  std::vector<XppVec>GetIntermediateSolutions();
  xpp_msgs::RobotParameters BuildRobotParametersMsg(const RobotModel& model) const;
  void SaveOptimizationAsRosbag(const std::string& bag_name,
                                const xpp_msgs::RobotParameters& robot_params,
                                const TowrCommandMsg user_command_msg,
                                bool include_iterations=false);
  void SaveTrajectoryInRosbag (rosbag::Bag&,
                               const std::vector<xpp::RobotStateCartesian>& traj,
                               const std::string& topic) const;
};

} /* namespace towr */

#endif /* TOWR_DRIVE_INCLUDE_TOWR_DRIVE_TOWR_DRIVE_ROS_INTERFACE_H_ */
