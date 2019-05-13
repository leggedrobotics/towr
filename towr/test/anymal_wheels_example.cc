/*
 * anymal_wheels_example.cc
 *
 *  Created on: Mar 29, 2019
 *      Author: vivian
 */

#include <cmath>
#include <iostream>

#include <towr/terrain/examples/height_map_examples.h>
#include <towr/initialization/gait_generator.h>
#include <towr/nlp_formulation_drive.h>
#include <ifopt/ipopt_solver.h>

using namespace towr;

void SetTowrParameters(NlpFormulationDrive *formulation)
{
  // set the initial position of the robot
  formulation->initial_base_.lin.at(kPos) << 0.0, 0.0, 0.5;

  // define the desired goal state of the robot
  formulation->final_base_.lin.at(kPos) << 0.5, 0.0, 0.5;

  formulation->params_drive_.total_time_ = 1.0;
  formulation->params_drive_.n_ee_ = 4;
}

int main()
{
  // NLP formulation
  NlpFormulationDrive formulation;

  // terrain
  formulation.terrain_ = std::make_shared<FlatGround>(0.0);

  // Kinematic limits and dynamic parameters of the Anymal with Wheels
  formulation.model_ = RobotModel(RobotModel::AnymalWheels);
  auto nominal_stance_B = formulation.model_.kinematic_model_->GetNominalStanceInBase();
  formulation.initial_ee_W_ = nominal_stance_B;
  std::for_each(formulation.initial_ee_W_.begin(), formulation.initial_ee_W_.end(),
				[&](Eigen::Vector3d& p){ p.z() = 0.0; } // feet at 0 height
  );

  // Load TOWR parameters
  SetTowrParameters(&formulation);

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
  solver->SetOption("linear_solver", "ma57");
  solver->SetOption("jacobian_approximation", "exact"); // "finite difference-values"
  solver->SetOption("max_cpu_time", 300.0);
  solver->Solve(nlp);

  nlp.PrintCurrent();
  formulation.PrintSolution(solution, 0.1);
}


