/*
 * towr_drive_app.cc
 *
 *  Created on: Jun 11, 2019
 *      Author: vivian
 */

#include <towr_drive/towr_drive_ros_interface.h>

#include <towr/initialization/gait_generator.h>

#include <xpp_states/convert.h>

#include "yaml_tools/yaml_tools.hpp"

namespace towr {

class TowrDriveRosApp : public TowrDriveRosInterface {
public:

  TowrDriveRosApp (std::string config_file)
  {
	config_file_ = config_file;

	initial_position_ = Eigen::Vector3d::Zero();

	goal_geom_.lin.p_.setZero();
	goal_geom_.ang.p_ << 0.0, 0.0, 0.0; // roll, pitch, yaw angle applied Z->Y'->X''

	shift_initialization_ = false;
  }

  static Eigen::Vector3d toEulerAngle(const Eigen::Quaterniond& q)
  {
	double roll, pitch, yaw;

  	// roll (x-axis rotation)
  	double sinr_cosp = +2.0 * (q.w() * q.x() + q.y() * q.z());
  	double cosr_cosp = +1.0 - 2.0 * (q.x() * q.x() + q.y() * q.y());
  	roll = atan2(sinr_cosp, cosr_cosp);

  	// pitch (y-axis rotation)
  	double sinp = +2.0 * (q.w() * q.y() - q.z() * q.x());
  	if (fabs(sinp) >= 1)
  		pitch = copysign(M_PI / 2, sinp); // use 90 degrees if out of range
  	else
  		pitch = asin(sinp);

  	// yaw (z-axis rotation)
  	double siny_cosp = +2.0 * (q.w() * q.z() + q.x() * q.y());
  	double cosy_cosp = +1.0 - 2.0 * (q.y() * q.y() + q.z() * q.z());
  	yaw = atan2(siny_cosp, cosy_cosp);

  	return Eigen::Vector3d(roll, pitch, yaw);
  }

  Eigen::Matrix3d GetRotationMatrixBaseToWorld (const Eigen::Vector3d& xyz)
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

  void SetTowrInitialStateFromController(const xpp::RobotStateCartesian& init_state)
  {
	formulation_.initial_base_.lin.at(kPos) = init_state.base_.lin.GetByIndex(xpp::kPos);
	formulation_.initial_base_.ang.at(kPos) = toEulerAngle(init_state.base_.ang.q);

	// initial feet position
	int n_ee = formulation_.model_.kinematic_model_->GetNumberOfEndeffectors();
	formulation_.initial_ee_W_.resize(n_ee);
	for (int ee = 0; ee < n_ee; ++ee) {
	  formulation_.initial_ee_W_.at(ee) = init_state.ee_motion_.at(ee).GetByIndex(xpp::kPos);
	}

  }

  void SetTowrInitialStateFromFile()
  {
	formulation_.initial_base_.lin.at(kPos) = initial_position_;
//	formulation_.initial_base_.ang.at(kPos) = Eigen::Vector3d(0.0, 0.0, 0.3491);

    // initial feet position
    int n_ee = formulation_.model_.kinematic_model_->GetNumberOfEndeffectors();
    auto nominal_stance_B = formulation_.model_.kinematic_model_->GetNominalStanceInBase();
    formulation_.initial_ee_W_.resize(n_ee);
    for (int ee = 0; ee < n_ee; ++ee) {
  	  formulation_.initial_ee_W_.at(ee) = GetRotationMatrixBaseToWorld(formulation_.initial_base_.ang.at(kPos)) * nominal_stance_B.at(ee)
  										  + formulation_.initial_base_.lin.at(kPos);
  	  Eigen::Vector3d ee_pos = formulation_.initial_ee_W_.at(ee);
  	  formulation_.initial_ee_W_.at(ee)(Z) = formulation_.terrain_->GetHeight(ee_pos.x(), ee_pos.y());
    }
  }


  ParametersDrive GetTowrDriveParameters(int n_ee, const TowrCommandMsg& msg) const override
  {
    ParametersDrive params_drive;

    yaml_tools::YamlNode basenode = yaml_tools::YamlNode::fromFile(config_file_);

    if (basenode.isNull())
  	throw std::runtime_error("CONFIGURATION LOADING FAILED");

    auto terrain_id = static_cast<HeightMap::TerrainID>(msg.terrain);
    std::string terrain = terrain_names.find(terrain_id)->second;

    // Number of end-effectors
    params_drive.n_ee_ = n_ee;

    // Time duration of the motion
    double total_duration = basenode[terrain]["total_time"].as<double>();
    params_drive.total_time_ = total_duration;

    bool use_wheels_motion_constraint = basenode[terrain]["use_wheels_motion_constraint"].as<bool>();
    if (use_wheels_motion_constraint)
    	params_drive.SetWheelsMotionConstraint();
    else
    	params_drive.SetEndeffectorRomConstraint();

    bool use_non_holonomic_constraint = basenode["use_non_holonomic_constraint"].as<bool>();
    if (use_non_holonomic_constraint)
  	  params_drive.SetNonHolonomicConstraint();

    params_drive.force_limit_in_x_direction_ = basenode[terrain]["traction_limit"].as<double>();

    bool constrain_initial_y_base = basenode["constrain_final_y_base"].as<bool>();
    if (constrain_initial_y_base)
    {
    	params_drive.bounds_initial_lin_pos_ = {X, Y, Z};
    }

    bool constrain_final_y_base = basenode["constrain_final_y_base"].as<bool>();
    bool constrain_final_z_base = basenode["constrain_final_z_base"].as<bool>();
    if (constrain_final_y_base)
    {
    	if (constrain_final_z_base)
    		params_drive.bounds_final_lin_pos_ = {X, Y, Z};
    	else
    		params_drive.bounds_final_lin_pos_ = {X, Y};
    }
    else
    {
    	if (constrain_final_z_base)
    		params_drive.bounds_final_lin_pos_ = {X, Z};
    	else
    		params_drive.bounds_final_lin_pos_ = {X};
    }

    bool constrain_final_base_orientation = basenode["constrain_final_base_orientation"].as<bool>();
    if (constrain_final_base_orientation)
    	params_drive.bounds_final_ang_pos_ = {X,Y,Z};
    else
    	params_drive.bounds_final_ang_pos_.clear();

    params_drive.constrain_final_ee_pos_ = basenode["constrain_final_ee_pos"].as<bool>();

    return params_drive;
  }

  towr_ros::TowrCommand BuildTowrCommandMsg (const xpp::RobotStateCartesian& init_state) override
  {
    towr_ros::TowrCommand msg;

    yaml_tools::YamlNode basenode = yaml_tools::YamlNode::fromFile(config_file_);

    if (basenode.isNull())
  	throw std::runtime_error("CONFIGURATION LOADING FAILED");
    Eigen::Vector3d final_position_;

    std::string terrain = basenode["terrain"].as<std::string>();
    auto towr_terrain_id = terrain_ids.find(terrain)->second;
    formulation_.terrain_ = HeightMap::MakeTerrain(towr_terrain_id);

    bool plot_trajectory = basenode["plot_trajectory"].as<bool>();
    bool play_initialization = basenode["play_initialization"].as<bool>();

    shift_initialization_ = basenode["shift_initialization"].as<bool>();

    // initial position now comes from the controller
    initial_position_.x() = basenode[terrain]["initial_pose"]["x"].as<double>();
    initial_position_.y() = basenode[terrain]["initial_pose"]["y"].as<double>();
    initial_position_.z() = basenode[terrain]["initial_pose"]["z"].as<double>();

    bool init_state_from_file = basenode["init_state_from_file"].as<bool>();
    if (init_state_from_file)
    	SetTowrInitialStateFromFile();
    else
    	SetTowrInitialStateFromController(init_state);

    // hack to initialize the left foot ahead of the right
    if (shift_initialization_)
    {
	  Eigen::Vector3d pos_phase (0.15, 0.0, 0.0);
	  formulation_.initial_ee_W_.at(LF) += pos_phase;
	  formulation_.initial_ee_W_.at(LH) += pos_phase;
	  formulation_.initial_ee_W_.at(RF) -= pos_phase;
	  formulation_.initial_ee_W_.at(RH) -= pos_phase;
    }

    final_position_.x() = basenode[terrain]["final_pose"]["x"].as<double>();
    final_position_.y() = basenode[terrain]["final_pose"]["y"].as<double>();
    final_position_.z() = basenode[terrain]["final_pose"]["z"].as<double>();
    goal_geom_.lin.p_ = final_position_;
    goal_geom_.lin.p_.x() += formulation_.initial_base_.lin.at(kPos).x();

    float total_duration = basenode[terrain]["total_time"].as<float>();

    msg.goal_lin                 = xpp::Convert::ToRos(goal_geom_.lin);
    msg.goal_ang                 = xpp::Convert::ToRos(goal_geom_.ang);
    msg.total_duration           = total_duration;
    msg.replay_trajectory        = true;
    msg.play_initialization      = play_initialization;
    msg.replay_speed             = 0.6; //1.0;
    msg.optimize                 = true;
    msg.terrain                  = (int) towr_terrain_id;
    msg.gait                     = towr::GaitGenerator::C0;   // it's not used for driving motions!!
    msg.robot                    = towr::RobotModel::AnymalWheels;
    msg.optimize_phase_durations = false;
    msg.plot_trajectory          = plot_trajectory;

    return msg;

  }

  /**
   * @brief Sets the paramters for IPOPT.
   */
  void SetIpoptParameters(const TowrCommandMsg& msg) override
  {
    solver_->SetOption("linear_solver", "ma97"); // ma27, ma57, ma77, ma86, ma97
	solver_->SetOption("jacobian_approximation", "exact"); // "finite difference-values"
	solver_->SetOption("max_cpu_time", 60.0); // 1 min
	solver_->SetOption("print_level", 5);

    if (msg.play_initialization)
      solver_->SetOption("max_iter", 0);
    else
      solver_->SetOption("max_iter", 3000);
  }

protected:
  std::string config_file_;
  xpp::State3dEuler goal_geom_;
  Eigen::Vector3d initial_position_;
  bool shift_initialization_;

};

} // namespace towr

int main(int argc, char *argv[])
{
  ros::init(argc, argv, "towr_drive_app");

  std::string config_file = ros::package::getPath("towr_drive") + "/config/parameters.yaml";
  towr::TowrDriveRosApp towr_drive_app (config_file);

  ros::spin();

  return 1;
}



