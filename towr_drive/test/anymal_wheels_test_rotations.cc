/*
 * anymal_wheels_test_rotations.cc
 *
 *  Created on: Apr 10, 2019
 *      Author: vivian
 */

#include <cmath>
#include <iostream>

#include <ros/ros.h>
#include <ros/package.h>

#include "yaml_tools/yaml_tools.hpp"

#include <towr/terrain/examples/height_map_examples.h>
#include <towr/initialization/gait_generator.h>
#include <towr/nlp_formulation_drive.h>
#include <ifopt/ipopt_solver.h>

#include "towr_ros/helpers_towr.h"

using namespace towr;

Eigen::Matrix3d GetRotationMatrixBaseToWorld (const Eigen::Vector3d xyz)
{
  double x = xyz(X);
  double y = xyz(Y);
  double z = xyz(Z);

  Eigen::Matrix3d M;
  //  http://docs.leggedrobotics.com/kindr/cheatsheet_latest.pdf (Euler ZYX)
  M << cos(y)*cos(z), cos(z)*sin(x)*sin(y) - cos(x)*sin(z), sin(x)*sin(z) + cos(x)*cos(z)*sin(y),
       cos(y)*sin(z), cos(x)*cos(z) + sin(x)*sin(y)*sin(z), cos(x)*sin(y)*sin(z) - cos(z)*sin(x),
             -sin(y),                        cos(y)*sin(x),                        cos(x)*cos(y);

  return M;
}

/**
* @brief Sets the parameters required to formulate the TOWR problem.
*/
bool SetTowrParameters(NlpFormulationDrive *formulation, const std::string& filename, const std::string terrain)
{
  yaml_tools::YamlNode basenode = yaml_tools::YamlNode::fromFile(filename);
  Eigen::Vector3d initial_position_, final_position_;

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
//	  std::cout << "Initial ee " << ee << " : " << formulation->initial_ee_W_.at(ee).transpose() << std::endl;
  }
//  formulation->initial_ee_W_ = nominal_stance_B;
//  std::for_each(formulation->initial_ee_W_.begin(), formulation->initial_ee_W_.end(),
//				[&](Eigen::Vector3d& p){ p.z() = 0.0; } // feet at 0 height
//  );

  // Time duration of the motion.
  double total_duration = basenode[terrain]["total_time"].as<double>();
  formulation->params_drive_.total_time_ = total_duration;

  return basenode["run_derivative_test"].as<bool>();

}

int main(int argc, char **argv)
{
  ros::init(argc, argv, "anymal_wheels_to");
  ros::NodeHandle nh;

  // NLP formulation
  NlpFormulationDrive formulation;

  // terrain
  HeightMap::TerrainID terrain_id = HeightMap::StepID;
  formulation.terrain_ = HeightMap::MakeTerrain(terrain_id);

  // Robot Model
  formulation.model_ = RobotModel(RobotModel::AnymalWheels);

  // Load parameters from file
  const std::string config_file = ros::package::getPath("towr_drive") + "/config/parameters.yaml";
  const std::string terrain_name = terrain_names.find(terrain_id)->second;
  bool run_derivative_test = SetTowrParameters(&formulation, config_file, terrain_name);

  auto ee_pos_W = formulation.model_.kinematic_model_->GetNominalStanceInBase();
  std::for_each(ee_pos_W.begin(), ee_pos_W.end(),
				[&](Eigen::Vector3d& p){ p.z() = 0.0; } // feet at 0 height
  );

  Eigen::Vector3d base_pos_W = formulation.initial_base_.lin.at(kPos);
  auto hips_pos_B = formulation.model_.kinematic_model_->GetHipPositionInBaseFrame();
  std::vector<Eigen::Vector3d> hips_pos_W(4);
  std::vector<double> legs(4, 0.0);

  std::cout << ">>> No rotation!!\n";
  std::cout << "Base pos: " << base_pos_W.transpose() << std::endl;
  for (int ee = 0; ee < 4; ee++) {
	  std::cout << "ee " << ee << ": " << ee_pos_W.at(ee).transpose() << std::endl;
	  hips_pos_W.at(ee) = hips_pos_B.at(ee) + base_pos_W;
	  std::cout << "hip " << ee << ": " << hips_pos_W.at(ee).transpose() << std::endl;
	  legs.at(ee) = (ee_pos_W.at(ee) - hips_pos_W.at(ee)).norm();
	  std::cout << "leg extension " << ee << ": " << legs.at(ee) << std::endl;
  }

  Eigen::Vector3d angles (0.0, 0.3, 0.0); // roll, pitch, yaw

  std::cout << "\n";
  std::cout << ">>> Rotation pitch = 17.18 deg!!\n";
  std::cout << "Base pos: " << base_pos_W.transpose() << std::endl;
  for (int ee = 0; ee < 4; ee++) {
	  std::cout << "ee " << ee << ": " << ee_pos_W.at(ee).transpose() << std::endl;
	  Eigen::Matrix3d w_R_b = GetRotationMatrixBaseToWorld(angles);//.transpose();
	  Eigen::Vector3d pos_hip_W = w_R_b*(hips_pos_B.at(ee)) + base_pos_W;
	  std::cout << "hip " << ee << ": " << pos_hip_W.transpose() << std::endl;
	  legs.at(ee) = (ee_pos_W.at(ee) - pos_hip_W).norm();
	  std::cout << "leg extension " << ee << ": " << legs.at(ee) << std::endl;
  }

  std::vector<double> v{3, 7, 4, 1, 5, 9};
  auto it = std::min_element(std::begin(v), std::end(v));
  auto min_idx = std::distance(std::begin(v), it);
  std::cout << "min element at: " << min_idx << ", value: " << *it << std::endl;

  Eigen::Vector3d x (1, 2, 3);
  Eigen::Vector3d y (4, 5, 6);
  Eigen::Vector3d z (7, 8, 9);
  Eigen::Matrix3d M;
  M << x, y, z;

  std::cout << M << std::endl;


  return 0;
}



