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

#include <towr/nlp_formulation.h>

namespace towr {

NlpFormulation::NlpFormulation ()
{
  using namespace std;
  cout << "\n";
  cout << "************************************************************\n";
  cout << " TOWR - Trajectory Optimization for Walking Robots (v1.4)\n";
  cout << "                \u00a9 Alexander W. Winkler\n";
  cout << "           https://github.com/ethz-adrl/towr\n";
  cout << "************************************************************";
  cout << "\n\n";
}

NlpFormulation::VariablePtrVec
NlpFormulation::GetVariableSets (SplineHolder& spline_holder)
{
  VariablePtrVec vars;

  auto base_motion = MakeBaseVariables();
  vars.insert(vars.end(), base_motion.begin(), base_motion.end());

  auto ee_motion = MakeEndeffectorVariables();
  vars.insert(vars.end(), ee_motion.begin(), ee_motion.end());

  auto ee_force = MakeForceVariables();
  vars.insert(vars.end(), ee_force.begin(), ee_force.end());

  auto contact_schedule = MakeContactScheduleVariables();
  // can also just be fixed timings that aren't optimized over, but still added
  // to spline_holder.
  if (params_.IsOptimizeTimings()) {
    vars.insert(vars.end(), contact_schedule.begin(), contact_schedule.end());
  }

  auto ee_decision = MakeDecisionVariables();

  // stores these readily constructed spline
  spline_holder = SplineHolder(base_motion.at(0), // linear
                               base_motion.at(1), // angular
                               params_.GetBasePolyDurations(),
                               ee_motion,
                               ee_force,
                               ee_decision,
                               contact_schedule,
                               params_.IsOptimizeTimings());
  return vars;
}

std::vector<NodesVariables::Ptr>
NlpFormulation::MakeBaseVariables () const
{
  std::vector<NodesVariables::Ptr> vars;

  int n_nodes = params_.GetBasePolyDurations().size() + 1;

  auto spline_lin = std::make_shared<NodesVariablesAll>(n_nodes, k3D, id::base_lin_nodes);

  double x = final_base_.lin.p().x();
  double y = final_base_.lin.p().y();
  double z = terrain_->GetHeight(x,y) - model_.kinematic_model_->GetNominalStanceInBase().front().z();
  Vector3d final_pos(x, y, z);
  double x2 = initial_base_.lin.p().x();
  double y2 = initial_base_.lin.p().y();
  double z2 = terrain_->GetHeight(x2,y2) - model_.kinematic_model_->GetNominalStanceInBase().front().z();
  Vector3d init_pos(x2, y2, z2);

  spline_lin->AdvancedInititialisationBase(
      init_pos, final_pos, params_.GetTotalTime(),
      params_.duration_base_polynomial_, final_base_v_.ang.p().z(),
      final_base_v_.lin.p().x(),
      final_base_v_.lin.p().y(), initial_base_.ang.p().z());


  Eigen::Vector3d euler(0.0, 0.0, final_base_.ang.p().z());
  Eigen::Matrix3d w_R_b = EulerConverter::GetRotationMatrixBaseToWorld(euler);
  Eigen::Vector3d desv(final_base_v_.lin.p().x(), final_base_v_.lin.p().y(), 0.0);
  Vector3d desv_rotated = w_R_b*desv;

  spline_lin->AddStartBound(kPos, {X,Y,Z}, init_pos);
  spline_lin->AddStartBound(kVel, {X,Y,Z}, final_base_v_.lin.p());
  spline_lin->AddFinalBound(kPos, params_.bounds_final_lin_pos_,   final_pos);
  spline_lin->AddFinalBound(kVel, params_.bounds_final_lin_vel_, desv_rotated);
  vars.push_back(spline_lin);

  auto spline_ang = std::make_shared<NodesVariablesAll>(n_nodes, k3D, id::base_ang_nodes);
  spline_ang->SetByLinearInterpolation(initial_base_.ang.p(), final_base_.ang.p(), params_.GetTotalTime());
  spline_ang->AddStartBound(kPos, {X,Y,Z}, initial_base_.ang.p());
  spline_ang->AddStartBound(kVel, {X,Y,Z}, final_base_v_.ang.p());
  spline_ang->AddFinalBound(kPos, params_.bounds_final_ang_pos_, final_base_.ang.p());
  spline_ang->AddFinalBound(kVel, params_.bounds_final_ang_vel_, final_base_v_.ang.p());
  vars.push_back(spline_ang);

  return vars;
}

std::vector<NodesVariablesPhaseBased::Ptr>
NlpFormulation::MakeEndeffectorVariables ()
{
  std::vector<NodesVariablesPhaseBased::Ptr> vars;

  // Endeffector Motions
  double T = params_.GetTotalTime();
  for (int ee=0; ee<params_.GetEECount(); ee++) {


    // initialize towards final footholds
    double yaw_final = final_base_.ang.p().z();
    double yaw_init = initial_base_.ang.p().z();
    Eigen::Vector3d euler_final(0.0, 0.0, yaw_final);
    Eigen::Vector3d euler_init(0.0, 0.0, yaw_init);
    Eigen::Matrix3d w_R_b_final =
        EulerConverter::GetRotationMatrixBaseToWorld(euler_final);
    Eigen::Matrix3d w_R_b_init =
        EulerConverter::GetRotationMatrixBaseToWorld(euler_init);

    double x = final_base_.lin.p().x();
    double y = final_base_.lin.p().y();
    double z = terrain_->GetHeight(x, y) -
               model_.kinematic_model_->GetNominalStanceInBase().front().z();
    Vector3d final_pos(x, y, z);


    double x2 = initial_base_.lin.p().x();
    double y2 = initial_base_.lin.p().y();
    double z2 = terrain_->GetHeight(x2,y2) - model_.kinematic_model_->GetNominalStanceInBase().front().z();
    Vector3d init_pos_base(x2, y2, z2);

    Eigen::Vector3d desv(final_base_v_.lin.p().x(), final_base_v_.lin.p().y(), 0.0);
    Vector3d desv_final = w_R_b_final*desv;


    Vector3d final_ee_pos_W =
        final_pos + w_R_b_final * (model_.kinematic_model_->GetNominalStanceInBase().at(ee));
    Vector3d init_ee_pos_W =
        init_pos_base+w_R_b_init *(model_.kinematic_model_->GetNominalStanceInBase().at(ee) );

     if(final_base_v_.lin.p().y()==0 && final_base_v_.ang.p().z()==0 ) {
//       std::cout << "test here do" << std::endl;
//       double total_t = params_.GetTotalTime();
//       params_.ee_phase_durations_.at(ee).clear();
//
//       params_.number_of_polys_per_phase_motion_.at(ee).clear();
//       params_.number_of_polys_per_phase_force_.at(ee).clear();
//       params_.number_of_polys_per_phase_decision_.at(ee).clear();
//       double current_t = 0;
//       bool contact = true;
//
//       double padding = 0.03;
//       double lrdev = 0.0;
//
//       double xasdf =
//           init_ee_pos_W.x() +
//           (final_ee_pos_W.x() - init_ee_pos_W.x()) * (current_t / total_t);
//       double yasdf =
//           init_ee_pos_W.y() +
//           (final_ee_pos_W.y() - init_ee_pos_W.y()) * (current_t / total_t);
//       double z_terrain_prev = terrain_->GetHeight(xasdf, yasdf);
//       double z_terrain = z_terrain_prev;
//
//       std::vector<double> durations;
//       std::vector<int> polys_per_phase_motion;
//       std::vector<int> polys_per_phase_force;
//       std::vector<int> polys_per_phase_decision;
//       double t_last = 0;
//       double walking_total = 0.0;
//
//       bool added_at_least_one = false;
//
//       double max_stance_time = 0.3;
//
//       while (current_t < total_t) {
//         xasdf = init_ee_pos_W.x() + (final_ee_pos_W.x() - init_ee_pos_W.x()) *
//                                         (current_t / total_t);
//         yasdf = init_ee_pos_W.y() + (final_ee_pos_W.y() - init_ee_pos_W.y()) *
//                                         (current_t / total_t);
//         z_terrain = terrain_->GetHeight(xasdf, yasdf);
//
//         if (contact) {
//
//           if (z_terrain != z_terrain_prev) {
//             contact = false;
//             double starttime = current_t - padding;
//             if (ee==0 || ee==2) starttime -= lrdev;
//             if (starttime - t_last < 0) {
//               starttime = t_last;
//             }
//             double interval_duration_s = starttime - t_last;
//             durations.emplace_back(interval_duration_s);
//             int mnodes =
//                 params_.motion_stance_nodes_per_s * interval_duration_s;
//             if (mnodes < params_.polynomials2_motion_per_stance_phase_) {
//               mnodes = params_.polynomials2_motion_per_stance_phase_;
//             }
//             int fnodes = params_.force_stance_nodes_per_s * interval_duration_s;
//             if (fnodes < params_.polynomials2_decision_per_stance_phase_) {
//               fnodes = params_.polynomials2_decision_per_stance_phase_;
//             }
//             int dnodes =
//                 params_.decision_stance_nodes_per_s * interval_duration_s;
//             if (dnodes < params_.polynomials2_force_per_stance_phase_) {
//               dnodes = params_.polynomials2_force_per_stance_phase_;
//             }
//             polys_per_phase_motion.emplace_back(mnodes);
//             polys_per_phase_force.emplace_back(fnodes);
//             polys_per_phase_decision.emplace_back(dnodes);
//             walking_total += interval_duration_s;
//             added_at_least_one = true;
//             t_last = starttime;
//           }
//
//         } else {
//           if (z_terrain == z_terrain_prev) {
//             contact = true;
//             double t_diff_min = params_.bound_phase_duration_.first;
//             double starttime = current_t + padding;
////             if (ee==0 || ee==2) starttime += 0.1;
////             else starttime += 0.2;
//             double t_diff = starttime - t_last;
//             if (t_diff < t_diff_min) {
//               t_diff = t_diff_min;
//             }
//             durations.emplace_back(t_diff);
//             polys_per_phase_motion.emplace_back(
//                 params_.polynomials2_motion_per_swing_phase_);
//             polys_per_phase_force.emplace_back(
//                 params_.polynomials2_force_per_swing_phase_);
//             polys_per_phase_decision.emplace_back(
//                 params_.polynomials2_decision_per_swing_phase_);
//             walking_total += (t_diff);
//             added_at_least_one = true;
//             t_last = t_diff + t_last;
//           }
//         }
//
//         z_terrain_prev = z_terrain;
//
//         current_t += 0.01;
//       }
//
//       if (!added_at_least_one) {
//         double interval_duration_s = total_t - walking_total;
//         durations.emplace_back(interval_duration_s);
//         walking_total += interval_duration_s;
//         int mnodes = params_.motion_stance_nodes_per_s * interval_duration_s;
//         if (mnodes < params_.polynomials2_motion_per_stance_phase_) {
//           mnodes = params_.polynomials2_motion_per_stance_phase_;
//         }
//         int fnodes = params_.force_stance_nodes_per_s * interval_duration_s;
//         if (fnodes < params_.polynomials2_decision_per_stance_phase_) {
//           fnodes = params_.polynomials2_decision_per_stance_phase_;
//         }
//         int dnodes = params_.decision_stance_nodes_per_s * interval_duration_s;
//         if (dnodes < params_.polynomials2_force_per_stance_phase_) {
//           dnodes = params_.polynomials2_force_per_stance_phase_;
//         }
//         polys_per_phase_motion.emplace_back(mnodes);
//         polys_per_phase_force.emplace_back(fnodes);
//         polys_per_phase_decision.emplace_back(dnodes);
//       }
//
//       double interval_duration_s = total_t - walking_total;
//       durations.emplace_back(interval_duration_s);
//       int mnodes = params_.motion_stance_nodes_per_s * interval_duration_s;
//       if (mnodes < 1) {
//         mnodes = 1;
//       }
//       int fnodes = params_.force_stance_nodes_per_s * interval_duration_s;
//       if (fnodes < 1) {
//         fnodes = 1;
//       }
//       int dnodes = params_.decision_stance_nodes_per_s * interval_duration_s;
//       if (dnodes < 1) {
//         dnodes = 1;
//       }
//       polys_per_phase_motion.emplace_back(mnodes);
//       polys_per_phase_force.emplace_back(fnodes);
//       polys_per_phase_decision.emplace_back(dnodes);
//
//       params_.number_of_polys_per_phase_motion_.at(ee) = polys_per_phase_motion;
//       params_.number_of_polys_per_phase_force_.at(ee) = polys_per_phase_force;
//       params_.number_of_polys_per_phase_decision_.at(ee) = polys_per_phase_decision;
//
//       params_.ee_phase_durations_.at(ee) = durations;

       // manually add durations
       std::vector<double> durations;
       durations.clear();
       std::vector<int> polys_per_phase_motion;
       std::vector<int> polys_per_phase_force;
       std::vector<int> polys_per_phase_decision;

       // For single stair case
       if (ee==0) { durations.emplace_back(0.35);durations.emplace_back(0.2);durations.emplace_back(1.85);}
       if (ee==1) { durations.emplace_back(0.55);durations.emplace_back(0.2);durations.emplace_back(1.65);}
       if (ee==2) { durations.emplace_back(1.45);durations.emplace_back(0.4);durations.emplace_back(0.55);}
       if (ee==3) { durations.emplace_back(1.65);durations.emplace_back(0.4);durations.emplace_back(0.35);}

       polys_per_phase_motion.emplace_back(std::max(int(params_.motion_stance_nodes_per_s * durations.at(0)),params_.polynomials2_motion_per_stance_phase_));
       polys_per_phase_motion.emplace_back(params_.polynomials2_motion_per_swing_phase_);
       polys_per_phase_motion.emplace_back(std::max(int(params_.motion_stance_nodes_per_s * durations.at(2)),1));

       polys_per_phase_force.emplace_back((std::max(int(params_.force_stance_nodes_per_s * durations.at(0)),params_.polynomials2_force_per_stance_phase_)));
       polys_per_phase_force.emplace_back(params_.polynomials2_force_per_swing_phase_);
       polys_per_phase_force.emplace_back((std::max(int(params_.force_stance_nodes_per_s * durations.at(2)),1)));

       polys_per_phase_decision.emplace_back((std::max(int(params_.decision_stance_nodes_per_s * durations.at(0)),params_.polynomials2_decision_per_stance_phase_)));
       polys_per_phase_decision.emplace_back(params_.polynomials2_decision_per_swing_phase_);
       polys_per_phase_decision.emplace_back((std::max(int(params_.decision_stance_nodes_per_s * durations.at(2)),1)));

       params_.ee_phase_durations_.at(ee).clear();

       params_.number_of_polys_per_phase_motion_.at(ee).clear();
       params_.number_of_polys_per_phase_force_.at(ee).clear();
       params_.number_of_polys_per_phase_decision_.at(ee).clear();

       params_.ee_phase_durations_.at(ee) = durations;
       params_.number_of_polys_per_phase_motion_.at(ee) = polys_per_phase_motion;
       params_.number_of_polys_per_phase_force_.at(ee) = polys_per_phase_force;
       params_.number_of_polys_per_phase_decision_.at(ee) = polys_per_phase_decision;

     }


    std::cout<<ee<<"   "<<params_.GetPhaseCount(ee)<<"   :"<<std::endl;
    for ( auto a: params_.ee_phase_durations_.at(ee)){
      std::cout<<a<<" ,";
    }
    std::cout<<std::endl;

    std::cout<<ee<<"  number_of_polys_per_phase_motion_  "<<params_.number_of_polys_per_phase_motion_.at(ee).size()<<"   :"<<std::endl;
    for ( auto a: params_.number_of_polys_per_phase_motion_.at(ee)){
      std::cout<<a<<" ,";
    }
    std::cout<<std::endl;
    std::cout<<ee<<"  number_of_polys_per_phase_force_  "<<params_.number_of_polys_per_phase_force_.at(ee).size()<<"   :"<<std::endl;
    for ( auto a: params_.number_of_polys_per_phase_force_.at(ee)){
      std::cout<<a<<" ,";
    }
    std::cout<<std::endl;
    std::cout<<ee<<"  number_of_polys_per_phase_decision_  "<<params_.number_of_polys_per_phase_decision_.at(ee).size()<<"   :"<<std::endl;
    for ( auto a: params_.number_of_polys_per_phase_decision_.at(ee)){
      std::cout<<a<<" ,";
    }
    std::cout<<std::endl;

    auto nodes = std::make_shared<NodesVariablesEEMotion>(
        params_.GetPhaseCount(ee),
        params_.ee_in_contact_at_start_.at(ee),
        id::EEMotionNodes(ee),        params_.number_of_polys_per_phase_motion_.at(ee));

    nodes->AdvancedInititialisationEE(
        init_ee_pos_W, final_ee_pos_W, params_.GetTotalTime(),
        params_.ee_phase_durations_.at(ee),
        params_.number_of_polys_per_phase_motion_.at(ee),
        final_base_v_.ang.p().z(),
        final_base_v_.lin.p().x(),
        final_base_v_.lin.p().y(),
        model_.kinematic_model_->GetNominalStanceInBase().at(ee), terrain_, initial_base_.ang.p().z(),
        params_.ee_in_contact_at_start_.at(ee));


    nodes->AddStartBound(kPos, {X,Y,Z}, init_ee_pos_W);
    nodes->AddStartBound(kVel, {X,Y,Z}, final_base_v_.lin.p());
    nodes->AddFinalBound(kPos, {X,Y,Z}, final_ee_pos_W);
    nodes->AddFinalBound(kVel, {X,Y,Z}, desv_final);
    vars.push_back(nodes);
  }

  return vars;
}

std::vector<NodesVariablesPhaseBased::Ptr>
NlpFormulation::MakeForceVariables () const
{
  std::vector<NodesVariablesPhaseBased::Ptr> vars;

  double T = params_.GetTotalTime();
  for (int ee=0; ee<params_.GetEECount(); ee++) {
    auto nodes = std::make_shared<NodesVariablesEEForce>(
                                              params_.GetPhaseCount(ee),
                                              params_.ee_in_contact_at_start_.at(ee),
                                              id::EEForceNodes(ee),
                                              params_.number_of_polys_per_phase_force_.at(ee));

    // initialize with mass of robot distributed equally on all legs
    double m = model_.dynamic_model_->m();
    double g = model_.dynamic_model_->g();

    Vector3d f_stance(0.0, 0.0, m*g/params_.GetEECount());
    nodes->SetByLinearInterpolation(f_stance, f_stance, T); // stay constant
    vars.push_back(nodes);
  }

  return vars;
}

std::vector<NodesVariablesPhaseBased::Ptr>
NlpFormulation::MakeDecisionVariables() const {
  std::vector<NodesVariablesPhaseBased::Ptr> vars;

  for (int ee = 0; ee < params_.GetEECount(); ee++) {
    auto nodes = std::make_shared<NodesVariablesEEDecision>(
        params_.GetPhaseCount(ee), params_.ee_in_contact_at_start_.at(ee),
        id::EEDecision(ee),  params_.number_of_polys_per_phase_decision_.at(ee));

    vars.push_back(nodes);
  }
  return vars;
}

std::vector<PhaseDurations::Ptr>
NlpFormulation::MakeContactScheduleVariables () const
{
  std::vector<PhaseDurations::Ptr> vars;

  for (int ee=0; ee<params_.GetEECount(); ee++) {
    auto var = std::make_shared<PhaseDurations>(ee,
                                                params_.ee_phase_durations_.at(ee),
                                                params_.ee_in_contact_at_start_.at(ee),
                                                params_.bound_phase_duration_.first,
                                                params_.bound_phase_duration_.second);
    vars.push_back(var);
  }

  return vars;
}

NlpFormulation::ConstraintPtrVec
NlpFormulation::GetConstraints(const SplineHolder& spline_holder) const
{
  ConstraintPtrVec constraints;
  for (auto name : params_.constraints_)
    for (auto c : GetConstraint(name, spline_holder))
      constraints.push_back(c);

  return constraints;
}

NlpFormulation::ConstraintPtrVec
NlpFormulation::GetConstraint (Parameters::ConstraintName name,
                           const SplineHolder& s) const
{
  switch (name) {
    case Parameters::Dynamic:                   return MakeDynamicConstraint(s);
    case Parameters::EndeffectorRom:            return MakeRangeOfMotionBoxConstraint(s);
    case Parameters::BaseRom:                   return MakeBaseRangeOfMotionConstraint(s);
    case Parameters::TotalTime:                 return MakeTotalTimeConstraint();
    case Parameters::Terrain:                   return MakeTerrainConstraint();
    case Parameters::Force:                     return MakeForceConstraint();
    case Parameters::Swing:                     return MakeSwingConstraint();
    case Parameters::BaseAcc:                   return MakeBaseAccConstraint(s);
    case Parameters::WheelsNonHolonomic:	return MakeWheelsNonHolonomicConstraint(s);
    case Parameters::TerrainDiscretized:        return MakeDiscretizedTerrainConstraint(s);
    case Parameters::ForceDiscretized:          return MakeDiscretizedForceConstraint(s);
    default: throw std::runtime_error("constraint not defined!");
  }
}


NlpFormulation::ConstraintPtrVec
NlpFormulation::MakeBaseRangeOfMotionConstraint (const SplineHolder& s) const
{
  return {std::make_shared<BaseMotionConstraint>(params_.GetTotalTime(),
                                                 params_.dt_constraint_base_motion_,
                                                 s)};
}

NlpFormulation::ConstraintPtrVec
NlpFormulation::MakeDynamicConstraint(const SplineHolder& s) const
{
  auto constraint = std::make_shared<DynamicConstraint>(model_.dynamic_model_,
                                                        params_.GetTotalTime(),
                                                        params_.dt_constraint_dynamic_,
                                                        s);
  return {constraint};
}

NlpFormulation::ConstraintPtrVec
NlpFormulation::MakeRangeOfMotionBoxConstraint (const SplineHolder& s) const
{
  ConstraintPtrVec c;

  for (int ee=0; ee<params_.GetEECount(); ee++) {
    auto rom = std::make_shared<RangeOfMotionConstraint>(model_.kinematic_model_,
                                                         params_.GetTotalTime(),
                                                         params_.dt_constraint_range_of_motion_,
                                                         ee,
                                                         s);
    c.push_back(rom);
  }

  return c;
}

NlpFormulation::ConstraintPtrVec
NlpFormulation::MakeTotalTimeConstraint () const
{
  ConstraintPtrVec c;
  double T = params_.GetTotalTime();

  for (int ee=0; ee<params_.GetEECount(); ee++) {
    auto duration_constraint = std::make_shared<TotalDurationConstraint>(T, ee);
    c.push_back(duration_constraint);
  }

  return c;
}

NlpFormulation::ConstraintPtrVec
NlpFormulation::MakeTerrainConstraint () const
{
  ConstraintPtrVec constraints;

  for (int ee=0; ee<params_.GetEECount(); ee++) {
    auto c = std::make_shared<TerrainConstraint>(terrain_, id::EEMotionNodes(ee));
    constraints.push_back(c);
  }

  return constraints;
}

NlpFormulation::ConstraintPtrVec
NlpFormulation::MakeForceConstraint () const
{
  ConstraintPtrVec constraints;

  for (int ee=0; ee<params_.GetEECount(); ee++) {
    auto c = std::make_shared<ForceConstraint>(terrain_,
                                               params_.force_limit_in_normal_direction_,
                                               ee);
    constraints.push_back(c);
  }

  return constraints;
}

NlpFormulation::ConstraintPtrVec
NlpFormulation::MakeSwingConstraint () const
{
  ConstraintPtrVec constraints;

  for (int ee=0; ee<params_.GetEECount(); ee++) {
    auto swing = std::make_shared<SwingConstraint>(id::EEMotionNodes(ee));
    constraints.push_back(swing);
  }

  return constraints;
}

NlpFormulation::ConstraintPtrVec
NlpFormulation::MakeBaseAccConstraint (const SplineHolder& s) const
{
  ConstraintPtrVec constraints;

  constraints.push_back(std::make_shared<SplineAccConstraint>
                        (s.base_linear_, id::base_lin_nodes));

  constraints.push_back(std::make_shared<SplineAccConstraint>
                        (s.base_angular_, id::base_ang_nodes));

  return constraints;
}

NlpFormulation::ConstraintPtrVec
NlpFormulation::GetCosts() const
{
  ConstraintPtrVec costs;
  for (const auto& pair : params_.costs_)
    for (auto c : GetCost(pair.first, pair.second))
      costs.push_back(c);

  return costs;
}

NlpFormulation::ConstraintPtrVec
NlpFormulation::MakeWheelsNonHolonomicConstraint (const SplineHolder& s) const
{
  ConstraintPtrVec c;

  for (int ee=0; ee<params_.GetEECount(); ee++) {
    auto constraint = std::make_shared<WheelsNonHolonomicConstraint>(terrain_, params_.GetTotalTime(),
                                                                     params_.dt_non_holonomic_, ee, s);
    c.push_back(constraint);
  }

  return c;
}

NlpFormulation::ConstraintPtrVec
NlpFormulation::MakeDiscretizedForceConstraint(const SplineHolder &s) const {
  ConstraintPtrVec c;

  for (int ee = 0; ee < params_.GetEECount(); ee++) {
    auto constraint = std::make_shared<ForceConstraintDiscretized>(
        terrain_, params_.GetTotalTime(), params_.dt_force_, ee,
        s, params_.force_limit_in_normal_direction_);
    c.push_back(constraint);
  }
  return c;
}

NlpFormulation::ConstraintPtrVec
NlpFormulation::MakeDiscretizedTerrainConstraint(const SplineHolder &s) const {
  ConstraintPtrVec c;

  for (int ee = 0; ee < params_.GetEECount(); ee++) {
    auto constraint = std::make_shared<TerrainConstraintDiscretized>(
        terrain_, params_.GetTotalTime(), params_.dt_terrain_discretized_, ee,
        s);
    c.push_back(constraint);
  }

  return c;
}

NlpFormulation::CostPtrVec
NlpFormulation::GetCost(const Parameters::CostName& name, double weight) const
{
  switch (name) {
    case Parameters::ForcesCostID:   return MakeForcesCost(weight);
    case Parameters::EEMotionCostID: return MakeEEMotionCost(weight);
    default: throw std::runtime_error("cost not defined!");
  }
}

NlpFormulation::CostPtrVec
NlpFormulation::MakeForcesCost(double weight) const
{
  CostPtrVec cost;

  for (int ee=0; ee<params_.GetEECount(); ee++)
    cost.push_back(std::make_shared<NodeCost>(id::EEForceNodes(ee), kPos, Z, weight));

  return cost;
}

NlpFormulation::CostPtrVec
NlpFormulation::MakeEEMotionCost(double weight) const
{
  CostPtrVec cost;

  for (int ee=0; ee<params_.GetEECount(); ee++) {
    cost.push_back(std::make_shared<NodeCost>(id::EEMotionNodes(ee), kVel, X, weight));
    cost.push_back(std::make_shared<NodeCost>(id::EEMotionNodes(ee), kVel, Y, weight));
  }

  return cost;
}

} /* namespace towr */
