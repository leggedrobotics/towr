/*
 * anymal_example.cc
 *
 *  Created on: Mar 19, 2019
 *      Author: vivian
 */
#include <cmath>
#include <iostream>

#include <towr/terrain/examples/height_map_examples.h>
#include <towr/initialization/gait_generator.h>
#include <towr/nlp_formulation.h>
#include <ifopt/ipopt_solver.h>


using namespace towr;

void PrintPhaseDurations(NlpFormulation *formulation)
{
  auto phase_durations = formulation->params_.ee_phase_durations_;

  std::cout << "Initial Phase Durations: \n";
  for (int ee = 0; ee < phase_durations.size(); ++ee) {
    std::cout << "EE " << ee << " : ";
    for (int i = 0; i < phase_durations.at(ee).size(); i++)
      std::cout << phase_durations.at(ee).at(i) << " ";
    std::cout << std::endl;
  }

}

// A minimal example to use TOWR with ANYmal.
int main()
{
  NlpFormulation formulation;

  // terrain
  formulation.terrain_ = std::make_shared<Gap>(); //(0.0);

  // Kinematic limits and dynamic parameters of the Anymal Wheels
  formulation.model_ = RobotModel(RobotModel::AnymalWheels);
  auto nominal_stance_B = formulation.model_.kinematic_model_->GetNominalStanceInBase();
  formulation.initial_ee_W_ = nominal_stance_B;
  std::for_each(formulation.initial_ee_W_.begin(), formulation.initial_ee_W_.end(),
                [&](Eigen::Vector3d& p){ p.z() = 0.0; } // feet at 0 height
  );

  // set the initial position of the robot
  formulation.initial_base_.lin.at(kPos).z() = 0.48;

  // define the desired goal state of the robot
  formulation.final_base_.lin.at(towr::kPos) << 2.0, 0.0, 0.48;

  // Parameters that define the motion.
  double total_duration = 2.0;
  auto id_gait = towr::GaitGenerator::C1; //GapCross; //C0;
  auto gait_gen_ = GaitGenerator::MakeGaitGenerator(4,id_gait, 2.0, formulation.terrain_);
  gait_gen_->SetCombo(id_gait);
  for (int ee=0; ee<4; ++ee) {
	formulation.params_.ee_phase_durations_.push_back(gait_gen_->GetPhaseDurations(total_duration, ee));
	formulation.params_.ee_in_contact_at_start_.push_back(gait_gen_->IsInContactAtStart(ee));
  }

  PrintPhaseDurations (&formulation);

  //formulation.params_.OptimizePhaseDurations();  // if optimization over timing durations is necessary..

  // Initialize the nonlinear-programming problem with the variables,
  // constraints and costs.
  ifopt::Problem nlp;
  SplineHolder solution;
  for (auto c : formulation.GetVariableSets(solution))
    nlp.AddVariableSet(c);
  for (auto c : formulation.GetConstraints(solution))
    nlp.AddConstraintSet(c);
  for (auto c : formulation.GetCosts())
    nlp.AddCostSet(c);

  // You can add your own elements to the nlp as well, simply by calling:
  // nlp.AddVariablesSet(your_custom_variables);
  // nlp.AddConstraintSet(your_custom_constraints);

  // Choose ifopt solver (IPOPT or SNOPT), set some parameters and solve.
  // solver->SetOption("derivative_test", "first-order");
  auto solver = std::make_shared<ifopt::IpoptSolver>();
  solver->SetOption("linear_solver", "ma97");
  solver->SetOption("jacobian_approximation", "exact"); // "finite difference-values"
  solver->SetOption("max_cpu_time", 60.0);
  solver->SetOption("max_iter", 3000);
  solver->SetOption("print_level", 5);

  // derivative test
  solver->SetOption("max_iter", 0);
  solver->SetOption("derivative_test", "first-order");
  solver->SetOption("print_level", 4);
  solver->SetOption("derivative_test_tol", 1e-3);

  solver->Solve(nlp);

  nlp.PrintCurrent();

}




