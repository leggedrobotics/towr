/*
 * anymal_wheels_test.cc
 *
 *  Created on: Apr 2, 2019
 *      Author: vivian
 */

#include <cmath>
#include <iostream>
#include <fstream>

#include <ros/ros.h>
#include <ros/package.h>

#include "yaml_tools/yaml_tools.hpp"

#include <towr/terrain/examples/height_map_examples.h>
#include <towr/initialization/gait_generator.h>
#include <towr/nlp_formulation_drive.h>
#include <ifopt/ipopt_solver.h>

#include <xpp_msgs/TerrainInfo.h>

#include "towr_ros/helpers_towr.h"

using namespace towr;

/**
* @brief Sets the parameters required to formulate the TOWR problem.
*/
bool SetTowrParameters(NlpFormulationDrive *formulation, const std::string& filename, int& terrain_id)
{
  yaml_tools::YamlNode basenode = yaml_tools::YamlNode::fromFile(filename);
  Eigen::Vector3d initial_position_, final_position_;

  std::string terrain = basenode["terrain"].as<std::string>();
  auto towr_terrain_id = terrain_ids.find(terrain)->second;
  formulation->terrain_ = HeightMap::MakeTerrain(towr_terrain_id);
  terrain_id = (int) towr_terrain_id;

  if (basenode.isNull())
	throw std::runtime_error("CONFIGURATION LOADING FAILED");

  initial_position_.x() = basenode[terrain]["initial_pose"]["x"].as<double>();
  initial_position_.y() = basenode[terrain]["initial_pose"]["y"].as<double>();
  initial_position_.z() = basenode[terrain]["initial_pose"]["z"].as<double>();

  final_position_.x() = basenode[terrain]["final_pose"]["x"].as<double>();
  final_position_.y() = basenode[terrain]["final_pose"]["y"].as<double>();
  final_position_.z() = basenode[terrain]["final_pose"]["z"].as<double>();

	// set the initial position of the robot
  formulation->initial_base_.lin.at(kPos) = initial_position_;
  formulation->initial_base_.ang.at(kPos) = Eigen::Vector3d(0.0, 0.0, 0.3491);

  // define the desired goal state of the robot
  formulation->final_base_.lin.at(kPos) = final_position_;

  // Number of end-effectors
  int n_ee = 4;
  formulation->params_drive_.n_ee_ = n_ee;

  // initial feet position
  auto nominal_stance_B = formulation->model_.kinematic_model_->GetNominalStanceInBase();
  formulation->initial_ee_W_.resize(n_ee);
  for (int ee = 0; ee < n_ee; ++ee) {
	  formulation->initial_ee_W_.at(ee) = nominal_stance_B.at(ee) + formulation->initial_base_.lin.at(kPos);
	  Eigen::Vector3d ee_pos = formulation->initial_ee_W_.at(ee);
	  formulation->initial_ee_W_.at(ee)(Z) = formulation->terrain_->GetHeight(ee_pos.x(), ee_pos.y());
  }

  bool shift_initialization = basenode["shift_initialization"].as<bool>();
  // hack to initialize the left foot ahead of the right
  if (shift_initialization)
  {
	Eigen::Vector3d pos_phase (0.15, 0.0, 0.0);
	formulation->initial_ee_W_.at(LF) += pos_phase;
	formulation->initial_ee_W_.at(LH) += pos_phase;
	formulation->initial_ee_W_.at(RF) -= pos_phase;
	formulation->initial_ee_W_.at(RH) -= pos_phase;
  }

  // Time duration of the motion
  double total_duration = basenode[terrain]["total_time"].as<double>();
  formulation->params_drive_.total_time_ = total_duration;

  bool use_wheels_motion_constraint = basenode[terrain]["use_wheels_motion_constraint"].as<bool>();
  if (use_wheels_motion_constraint)
	formulation->params_drive_.SetWheelsMotionConstraint();
  else
	formulation->params_drive_.SetEndeffectorRomConstraint();

  bool use_non_holonomic_constraint = basenode["use_non_holonomic_constraint"].as<bool>();
  if (use_non_holonomic_constraint)
	formulation->params_drive_.SetNonHolonomicConstraint();

  formulation->params_drive_.force_limit_in_x_direction_ = basenode[terrain]["traction_limit"].as<double>();

  return basenode["run_derivative_test"].as<bool>();

}

int main(int argc, char **argv)
{
  ros::init(argc, argv, "anymal_wheels_to");
  ros::NodeHandle nh;

  // NLP formulation
  NlpFormulationDrive formulation;

  // Robot Model
  formulation.model_ = RobotModel(RobotModel::AnymalWheels);

  // Load parameters from file
  const std::string config_file = ros::package::getPath("towr_drive") + "/config/parameters.yaml";
  int terrain_id = 0;
//  std::cout << "terrain_id (antes): " << terrain_id << std::endl;
  bool run_derivative_test = SetTowrParameters(&formulation, config_file, terrain_id);
//  std::cout << "terrain_id (depois): " << terrain_id << std::endl;

  // Initialize the nonlinear-programming problem with the variables,
  // constraints and costs.
  ifopt::Problem nlp;
  SplineHolderDrive solution;
  for (auto c : formulation.GetVariableSets(solution))
    nlp.AddVariableSet(c);
  for (auto c : formulation.GetConstraints(solution))
    nlp.AddConstraintSet(c);
  for (auto c : formulation.GetCosts())
    nlp.AddCostSet(c);

  auto solver = std::make_shared<ifopt::IpoptSolver>();
  solver->SetOption("linear_solver", "ma97"); // ma27, ma57, ma77, ma86, ma97
  solver->SetOption("jacobian_approximation", "exact"); // "finite difference-values"
  solver->SetOption("max_cpu_time", 180.0); // 3 min
  solver->SetOption("print_level", 5);

  // derivative test
  if (run_derivative_test) {
    solver->SetOption("max_iter", 0);
    solver->SetOption("derivative_test", "first-order");
    solver->SetOption("print_level", 4);
    solver->SetOption("derivative_test_tol", 1e-3);
    //solver->SetOption("derivative_test_perturbation", 1e-4);
    //solver->SetOption("derivative_test_print_all", "yes");
  }

  solver->Solve(nlp);

  // Print NLP information and some nodes solution
  nlp.PrintCurrent();
//  formulation.PrintSolution(solution, 0.5);

  // save trajectory (to send to the controller)
//  std::string bag_file = ros::package::getPath("towr_test") + "/bags/anymal_wheels_traj.bag";
  std::string bag_file = ros::package::getPath("anymal_wheels_ctrl_track_ros") + "/data/anymal_wheels_traj.bag";
  SaveDrivingMotionTerrainInRosbag (solution, terrain_id, bag_file);

  // save trajectory as collection of states (for plots)
  bag_file = ros::package::getPath("towr_drive") + "/bags/anymal_wheels_states.bag";
  SaveDrivingTrajectoryStatesInRosbag (solution, bag_file);

  // Create a new bag with geometry messages to plot in MATLAB
  ExtractGeometryMessagesFromTrajectoryBag(bag_file);

  // Create a map of terrain normals in a txt file (not necessary for the controller anymore!!)
//  std::string filename = ros::package::getPath("anymal_wheels_ctrl_track") + "/data/terrain_normals.txt";
//  SaveTerrainNormalsInFile(solution, terrain_id, filename);

  return 0;
}







