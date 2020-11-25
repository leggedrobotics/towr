/******************************************************************************
Copyright (c) 2018, Alexander W. Winkler. All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

* Neither the name of the copyright holder nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
******************************************************************************/

#include <towr/initialization/gait_generator.h>
#include <towr_ros/towr_ros_interface.h>

#include "yaml_tools/yaml_tools.hpp"


namespace towr {

/**
 * @brief An example application of using TOWR together with ROS.
 *
 * Build your own application with your own formulation using the building
 * blocks provided in TOWR and following the example below.
 */
class TowrRosApp : public TowrRosInterface {
public:
  using VecTimes = std::vector<double>;

  TowrRosApp (std::string config_file)
  {
	config_file_ = config_file;
  }

  void PrintPhaseDurations(std::vector<VecTimes> phase_durations) const
  {
    std::cout << "Initial Phase Durations: \n";
    for (int ee = 0; ee < phase_durations.size(); ++ee) {
      std::cout << "EE " << ee << " : ";
      for (int i = 0; i < phase_durations.at(ee).size(); i++)
        std::cout << phase_durations.at(ee).at(i) << " ";
      std::cout << std::endl;
    }
  }

  /**
   * @brief Sets the feet to nominal position on flat ground and base above.
   */
  void SetTowrInitialState() override
  {
    auto nominal_stance_B = formulation_.model_.kinematic_model_->GetNominalStanceInBase();

    double z_ground = 0.0;
    formulation_.initial_ee_W_ =  nominal_stance_B;
    std::for_each(formulation_.initial_ee_W_.begin(), formulation_.initial_ee_W_.end(),
                  [&](Vector3d& p){ p.z() = z_ground; } // feet at 0 height
    );

    formulation_.initial_base_.lin.at(kPos).z() = - nominal_stance_B.front().z() + z_ground;
  }

  /**
   * @brief Sets the parameters required to formulate the TOWR problem.
   */
  Parameters GetTowrParameters(int n_ee, const TowrCommandMsg& msg) const override
  {
    Parameters params;

    yaml_tools::YamlNode basenode = yaml_tools::YamlNode::fromFile(config_file_);
    if (basenode.isNull())
    	throw std::runtime_error("CONFIGURATION LOADING FAILED");

    auto terrain_id = static_cast<HeightMap::TerrainID>(msg.terrain);
    std::string terrain = terrain_names.find(terrain_id)->second;

    // Define gait
    auto gait_gen_ = GaitGenerator::MakeGaitGenerator(n_ee);
    auto id_gait   = static_cast<GaitGenerator::Combos>(msg.gait);
    gait_gen_->SetCombo(id_gait);
    for (int ee=0; ee<n_ee; ++ee) {
      params.ee_phase_durations_.push_back(gait_gen_->GetPhaseDurations(msg.total_duration, ee));
      params.ee_in_contact_at_start_.push_back(gait_gen_->IsInContactAtStart(ee));
    }

    int n_polys_in_swing_phase = basenode["n_polys_in_swing_phase"].as<int>();
    params.n_polynomials_per_swing_phase_ = n_polys_in_swing_phase;
    // These gaits (for floating steps) converge faster for 3 polynomials in swing phase
    // All the other gaits work well for 2 polynomials in swing phase
    // (this is optional, you can comment the if bellow to test other number of polynomials)
    // TODO: include this parameter in the user interface?
    if (msg.gait == 9 || msg.gait == 11)
    	params.n_polynomials_per_swing_phase_ = 3;
    else
    	params.n_polynomials_per_swing_phase_ = 2;

    // Print phase durations (optional)
    //PrintPhaseDurations(params.ee_phase_durations_);

    // increases optimization time, but sometimes helps find a solution for
    // more difficult terrain.
    if (msg.optimize_phase_durations)
      params.OptimizePhaseDurations();

    bool constrain_final_z_base = basenode["constrain_final_z_base"].as<bool>();
    if (constrain_final_z_base)
    {
    	params.bounds_final_lin_pos_ = {X, Y, Z};
    }

    // Get parameters specific to that terrain type
    bool use_non_holonomic_constraint = basenode[terrain]["use_non_holonomic_constraint"].as<bool>();
    if (use_non_holonomic_constraint)
  	  params.SetNonHolonomicConstraint();

    params.limit_base_angles_ = false;
    bool limit_base_angles = basenode[terrain]["limit_base_angles"].as<bool>();
    if (limit_base_angles)
  	  params.limit_base_angles_ = true;

    float min_distance_above_terrain = basenode["min_distance_above_terrain"].as<float>();
    params.min_distance_above_terrain_ = min_distance_above_terrain;

	double max_wheels_acc_z = basenode[terrain]["max_wheels_acc_z"].as<double>();
	params.max_wheels_acc_.at(Z) = max_wheels_acc_z;

//	    bool constrain_base_acc = basenode[terrain]["constrain_base_acc"].as<bool>();
//	    if (constrain_base_acc)
//	  	  params.SetBaseAccLimitsContraint();

    return params;
  }

  /**
   * @brief Sets the paramters for IPOPT.
   */
  void SetIpoptParameters(const TowrCommandMsg& msg) override
  {
	yaml_tools::YamlNode basenode = yaml_tools::YamlNode::fromFile(config_file_);

	if (basenode.isNull())
	throw std::runtime_error("CONFIGURATION LOADING FAILED");

	bool run_derivative_test = basenode["run_derivative_test"].as<bool>();
	solver_->SetOption("linear_solver", "ma57"); // ma27, ma57, ma77, ma86, ma97
	solver_->SetOption("jacobian_approximation", "exact"); // "finite difference-values"
	solver_->SetOption("max_cpu_time", 60.0); // 1 min
	solver_->SetOption("print_level", 5);

	if (msg.play_initialization)
	  solver_->SetOption("max_iter", 0);
	else
	  solver_->SetOption("max_iter", 3000);

    // derivative test
    if (run_derivative_test) {
      solver_->SetOption("max_iter", 0);
      solver_->SetOption("derivative_test", "first-order");
      solver_->SetOption("print_level", 4);
      solver_->SetOption("derivative_test_tol", 1e-3);
      //solver_->SetOption("derivative_test_perturbation", 1e-4);
      //solver_->SetOption("derivative_test_print_all", "yes");
    }

  }

protected:
  std::string config_file_;

};

} // namespace towr


int main(int argc, char *argv[])
{
  ros::init(argc, argv, "my_towr_ros_app");

  std::string config_file = ros::package::getPath("towr_ros") + "/config/parameters.yaml";
  towr::TowrRosApp towr_app (config_file);
  ros::spin();

  return 1;
}
