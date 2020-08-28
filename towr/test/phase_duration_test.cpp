//
// Created by chaoni on 8/15/20.
//


#include <gtest/gtest.h>
#include <cmath>
#include <iostream>

#include <towr/terrain/examples/height_map_examples.h>
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
  formulation.final_base_.lin.at(kPos) << 2.1,0,0.5;

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

  std::cout << "without smart initialization\n";
  for (int ee=0; ee<n_ee; ++ee){
    for ( auto a: params.ee_phase_durations_.at(ee)){
      std::cout<<a<<" ,";
    }
    std::cout<<"\n";
  }

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

  SplineHolder spline_holder;
  auto c = formulation.GetVariableSets(spline_holder);

//  std::cout << "with smart initialization\n";
//  for (int ee=0; ee<n_ee; ++ee){
//    for ( auto a: formulation.params_.ee_phase_durations_.at(ee)){
//      std::cout<<a<<" ,";
//    }
//    std::cout<<"\n";
//  }

//  std::cout << formulation.final_base_v_.lin.p() ;
}
