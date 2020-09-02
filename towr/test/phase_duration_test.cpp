//
// Created by chaoni on 8/15/20.
//


#include <gtest/gtest.h>
#include <cmath>
#include <iostream>
#include <fstream>

#include <towr/terrain/examples/height_map_examples.h>
#include <towr/variables/spline_holder.h>
#include <towr/variables/phase_spline.h>

#include <towr/nlp_formulation.h>
#include <ifopt/ipopt_solver.h>
#include <towr/initialization/gait_generator.h>


using namespace towr;

int main(){
  int n_ee = 4;
//  auto gait_gen_ = GaitGenerator::MakeGaitGenerator(n_ee);
//  auto id_gait   = static_cast<GaitGenerator::Combos>(1); //msg.gait=1: trot
//  gait_gen_->SetCombo(id_gait);
//  Parameters params;
//  for (int ee=0; ee < n_ee; ++ee){
//    params.ee_phase_durations_.push_back((gait_gen_->GetPhaseDurations(2.4,ee)));
//  }
//  std::cout << "without smart initialization\n";
//  for (int ee=0; ee<n_ee; ++ee){
//    for ( auto a: params.ee_phase_durations_.at(ee)){
//      std::cout<<a<<" ,";
//    }
//    std::cout<<"\n";
//  }

  NlpFormulation formulation;

  formulation.terrain_ = std::make_shared<StepFlat>();
  formulation.model_ = RobotModel(RobotModel::AnymalWheels);
  auto nominal_stance_B = formulation.model_.kinematic_model_->GetNominalStanceInBase();
//  std::cout << nominal_stance_B.front().z();

  formulation.initial_ee_W_ = nominal_stance_B;
  double z_ground = 0.0;
  std::for_each(formulation.initial_ee_W_.begin(), formulation.initial_ee_W_.end(), [&](Eigen::Vector3d &p){ p.z() = z_ground; });
  formulation.initial_base_.lin.at(kPos).z() = -nominal_stance_B.front().z() + z_ground;
  formulation.final_base_.lin.at(kPos) << 1.4,0,0.05;

//  SplineHolder spline_holder;
//  auto c = formulation.GetVariableSets(spline_holder);

  Parameters params;
  // Instead of manually defining the initial durations for each foot and
  // step, for convenience we use a GaitGenerator with some predefined gaits
  // for a variety of robots (walk, trot, pace, ...).
  auto gait_gen_ = GaitGenerator::MakeGaitGenerator(n_ee);
  auto id_gait   = static_cast<GaitGenerator::Combos>(1); // 1:trot
  gait_gen_->SetCombo(id_gait);
  for (int ee=0; ee<n_ee; ++ee) {
    params.ee_phase_durations_.push_back(gait_gen_->GetPhaseDurations(2.4, ee)); //total duration: 2.4s
    params.ee_in_contact_at_start_.push_back(gait_gen_->IsInContactAtStart(ee));
  }

//  std::cout << "without smart initialization\n";
//  for (int ee=0; ee<n_ee; ++ee){
//    for ( auto a: params.ee_phase_durations_.at(ee)){
//      std::cout<<a<<" ,";
//    }
//    std::cout<<"\n";
//  }

  params.number_of_polys_per_phase_motion_.clear();
  params.number_of_polys_per_phase_force_.clear();
  params.number_of_polys_per_phase_decision_.clear();

  for (int ee=0; ee<n_ee; ++ee) {
    bool contact = gait_gen_->IsInContactAtStart(ee);
    std::vector<int> temp_motion;
    std::vector<int> temp_force;
    std::vector<int> temp_decision;
    for (auto const &value : params.ee_phase_durations_.at(ee)) {
      if (contact) {
        temp_motion.push_back(params.polynomials2_motion_per_stance_phase_);
        temp_force.push_back(params.polynomials2_force_per_stance_phase_);
        temp_decision.push_back(params.polynomials2_decision_per_stance_phase_);
        contact = false;
      } else {
        temp_motion.push_back(params.polynomials2_motion_per_swing_phase_);
        temp_force.push_back(params.polynomials2_force_per_swing_phase_);
        temp_decision.push_back(params.polynomials2_decision_per_swing_phase_);
        contact = true;
      }
    }
    params.number_of_polys_per_phase_motion_.push_back(temp_motion);
    params.number_of_polys_per_phase_force_.push_back(temp_force);
    params.number_of_polys_per_phase_decision_.push_back(temp_decision);
  }

  formulation.params_ = params;


  std::cout << "here?" << std::endl;
  SplineHolder solution;
  ifopt::Problem nlp;
  for (auto c : formulation.GetVariableSets(solution))
    nlp.AddVariableSet(c);
  for (auto c : formulation.GetConstraints(solution))
    nlp.AddConstraintSet(c);
  for (auto c : formulation.GetCosts())
    nlp.AddCostSet(c);

//  auto solver_ = std::make_shared<ifopt::IpoptSolver>();
//  solver_->SetOption("linear_solver", "ma57");
//  solver_->SetOption("jacobian_approximation", "exact");
//  solver_->SetOption("max_cpu_time", 8000.0);
//  solver_->SetOption("print_level", 5);
//  solver_->SetOption("max_iter", 0);
//  solver_->Solve(nlp);
//
//  nlp.PrintCurrent();

  // debug: read date intio csv file
  std::string name = "Initdebug";
  std::ofstream myfile("/home/chaoni/git/towr/towr_ros/python/"+ name+ ".csv");

  double t = 0;
  while (t < solution.base_linear_->GetTotalTime()+ 1e-5){
    myfile << solution.ee_motion_.at(0)->GetPoint(t).p() << "\n";
    myfile << solution.ee_decision_.at(0)->GetPoint(t).p() << "\n";
    myfile << solution.phase_durations_.at(0)->IsContactPhase(t) << "\n";

    myfile << solution.ee_motion_.at(1)->GetPoint(t).p() << "\n";
    myfile << solution.ee_decision_.at(1)->GetPoint(t).p() << "\n";
    myfile << solution.phase_durations_.at(1)->IsContactPhase(t) << "\n";

    myfile << solution.ee_motion_.at(2)->GetPoint(t).p() << "\n";
    myfile << solution.ee_decision_.at(2)->GetPoint(t).p() << "\n";
    myfile << solution.phase_durations_.at(2)->IsContactPhase(t) << "\n";

    myfile << solution.ee_motion_.at(3)->GetPoint(t).p() << "\n";
    myfile << solution.ee_decision_.at(3)->GetPoint(t).p() << "\n";
    myfile << solution.phase_durations_.at(3)->IsContactPhase(t) << "\n";
    t += 0.02;
  }
  myfile.close();

  std::string name_1 = "NodesValue";
  std::ofstream nodefile("/home/chaoni/git/towr/towr_ros/python/"+ name_1+ ".csv");

  auto ee_motion = formulation.MakeEndeffectorVariables();
  for (int ee=0; ee< n_ee; ++ee){
    auto foot = ee_motion.at(ee)->GetNodes();
    for (int i=0; i<foot.size(); ++i)
      nodefile << foot.at(i).p() << "\n";
  }

  nodefile.close();

  auto contact_schedule = formulation.MakeContactScheduleVariables();

  /** to change the nodes value manually*/
  for (int ee=0; ee<n_ee; ++ee){
    std::cout << "EE: " << ee << " " << ee_motion.at(ee)->nodes_.size() << "\n" << std::endl;
  }
   ee_motion.at(0)->nodes_.at(5).at(kPos).z() = 0.05;
   ee_motion.at(1)->nodes_.at(5).at(kPos).z() = 0.05;
   ee_motion.at(2)->nodes_.at(6).at(kPos).z() = 0.05;
   ee_motion.at(3)->nodes_.at(7).at(kPos).z() = 0.05;



   // print out manual node values
  std::string name_3 = "ManualNodesValue";
  std::ofstream manualnodefile("/home/chaoni/git/towr/towr_ros/python/"+ name_3+ ".csv");

  for (int ee=0; ee< n_ee; ++ee){
    auto foot = ee_motion.at(ee)->GetNodes();
    for (int i=0; i<foot.size(); ++i)
      manualnodefile << foot.at(i).p() << "\n";
  }
  manualnodefile.close();

  // print manual ee positions
  std::vector<NodeSpline::Ptr> ee_motion_manual;
  for (int ee=0; ee<n_ee; ++ee){
    ee_motion_manual.push_back(std::make_shared<PhaseSpline>(ee_motion.at(ee), contact_schedule.at(ee).get()));
  }

  std::string name_2 = "manualEEPos";
  std::ofstream manualEEfile("/home/chaoni/git/towr/towr_ros/python/"+ name_2+ ".csv");

  t = 0;
  while (t < solution.base_linear_->GetTotalTime()+ 1e-5){
    manualEEfile << ee_motion_manual.at(0)->GetPoint(t).p() << "\n";
    manualEEfile << ee_motion_manual.at(1)->GetPoint(t).p() << "\n";
    manualEEfile << ee_motion_manual.at(2)->GetPoint(t).p() << "\n";
    manualEEfile << ee_motion_manual.at(3)->GetPoint(t).p() << "\n";
    t += 0.02;
  }

  manualEEfile.close();



  //  std::cout << formulation.final_base_v_.lin.p() ;
}
