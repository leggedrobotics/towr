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

  // stores these readily constructed spline
  spline_holder = SplineHolder(base_motion.at(0), // linear
                               base_motion.at(1), // angular
                               params_.GetBasePolyDurations(),
                               ee_motion,
                               ee_force,
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

  //spline_lin->SetByLinearInterpolation(initial_base_.lin.p(), final_pos, params_.GetTotalTime());

  Vector3d des_v((final_base_.lin.p().x() - initial_base_.lin.p().x()) / params_.GetTotalTime(),
                 (final_base_.lin.p().y() - initial_base_.lin.p().y()) / params_.GetTotalTime(),
                 (final_base_.lin.p().z() - initial_base_.lin.p().z()) / params_.GetTotalTime());

  double des_w = (final_base_.ang.p().z() - initial_base_.ang.p().z()) / params_.GetTotalTime();

 spline_lin->AdvancedInititialisationBase(
      init_pos, final_pos, params_.GetTotalTime(),
      params_.duration_base_polynomial_, des_w,
      des_v[0],
      des_v[1], initial_base_.ang.p().z(), terrain_, terrainID_ );



  spline_lin->AddStartBound(kPos, {X,Y,Z}, initial_base_.lin.p());
  spline_lin->AddStartBound(kVel, {X,Y,Z}, initial_base_.lin.v());
  spline_lin->AddFinalBound(kPos, params_.bounds_final_lin_pos_,   final_base_.lin.p());
  spline_lin->AddFinalBound(kVel, params_.bounds_final_lin_vel_, final_base_.lin.v());
  vars.push_back(spline_lin);

  auto spline_ang = std::make_shared<NodesVariablesAll>(n_nodes, k3D, id::base_ang_nodes);
  spline_ang->SetByLinearInterpolation(initial_base_.ang.p(), final_base_.ang.p(), params_.GetTotalTime());
  spline_ang->AddStartBound(kPos, {X,Y,Z}, initial_base_.ang.p());
  spline_ang->AddStartBound(kVel, {X,Y,Z}, initial_base_.ang.v());
  spline_ang->AddFinalBound(kPos, params_.bounds_final_ang_pos_, final_base_.ang.p());
  spline_ang->AddFinalBound(kVel, params_.bounds_final_ang_vel_, final_base_.ang.v());

  // limits all angles displacement to 30 degrees (optional)
  Vector3d ang_limit_ = Vector3d(10, 10, 10) * (M_PI/180);
  if (params_.limit_base_angles_)
	  spline_ang->AddAllNodesBounds(kPos, {Y}, -ang_limit_, ang_limit_);

  vars.push_back(spline_ang);

  return vars;
}

std::vector<NodesVariablesPhaseBased::Ptr>
NlpFormulation::MakeEndeffectorVariables () const
{
  std::vector<NodesVariablesPhaseBased::Ptr> vars;

  // Endeffector Motions
  double T = params_.GetTotalTime();

  std::vector<int> n_polys;

  for (int ee=0; ee<params_.GetEECount(); ee++) {

	bool phase_constant = params_.ee_in_contact_at_start_.at(ee);
	const std::vector<double> phase_durations = params_.ee_phase_durations_.at(ee);
	for (auto v : phase_durations) {
		if (phase_constant) {
			n_polys.push_back(floor(v/params_.dt_drive_constraint_)+1);
		}
		else {
			n_polys.push_back(params_.n_polynomials_per_swing_phase_);
		}
		phase_constant = !phase_constant;
	}

    auto nodes = std::make_shared<NodesVariablesEEMotion>(
                                              params_.GetPhaseCount(ee),
                                              params_.ee_in_contact_at_start_.at(ee),
                                              id::EEMotionNodes(ee),
                                              n_polys, !params_.use_non_holonomic_constraint_);

    // initialize towards final footholds
    double yaw = final_base_.ang.p().z();
    Eigen::Vector3d euler(0.0, 0.0, yaw);
    Eigen::Matrix3d w_R_b = EulerConverter::GetRotationMatrixBaseToWorld(euler);
    Vector3d final_ee_pos_W = final_base_.lin.p() + w_R_b*model_.kinematic_model_->GetNominalStanceInBase().at(ee);
    double x = final_ee_pos_W.x();
    double y = final_ee_pos_W.y();
    double z = terrain_->GetHeight(x,y);

    double yaw_init = initial_base_.ang.p().z();
    Eigen::Vector3d euler_init(0.0, 0.0, yaw_init);
    Eigen::Matrix3d w_R_b_init =
        EulerConverter::GetRotationMatrixBaseToWorld(euler_init);
    double x2 = initial_base_.lin.p().x();
    double y2 = initial_base_.lin.p().y();
    double z2 = terrain_->GetHeight(x2,y2) - model_.kinematic_model_->GetNominalStanceInBase().front().z();
    Vector3d init_pos_base(x2, y2, z2);
    Vector3d init_ee_pos_W =
        init_pos_base+w_R_b_init *(model_.kinematic_model_->GetNominalStanceInBase().at(ee) );

    //nodes->SetByLinearInterpolation(initial_ee_W_.at(ee), Vector3d(x,y,z), T);

    Vector3d des_v((final_base_.lin.p().x() - initial_base_.lin.p().x()) / params_.GetTotalTime(),
                   (final_base_.lin.p().y() - initial_base_.lin.p().y()) / params_.GetTotalTime(),
                   (final_base_.lin.p().z() - initial_base_.lin.p().z()) / params_.GetTotalTime());

    double des_w = (final_base_.ang.p().z() - initial_base_.ang.p().z()) / params_.GetTotalTime();

    nodes->AdvancedInititialisationEE(
        init_ee_pos_W, final_ee_pos_W, params_.GetTotalTime(),
        params_.ee_phase_durations_.at(ee),
        n_polys,
        des_w,
        des_v[0],
        des_v[1],
        model_.kinematic_model_->GetNominalStanceInBase().at(ee), terrain_,  initial_base_.ang.p().z(),
        params_.ee_in_contact_at_start_.at(ee),terrainID_);

    nodes->AddStartBound(kPos, {X,Y,Z}, initial_ee_W_.at(ee));
    nodes->AddStartBound(kVel, {X,Y,Z}, Vector3d(0, 0, 0));
    nodes->AddFinalBound(kVel, {X,Y,Z}, Vector3d(0, 0, 0));

    // bound on the final ee_pos
//    nodes->AddFinalBound(kPos, {X,Y,Z}, Vector3d(x,y,z));

    vars.push_back(nodes);

    n_polys.clear();
    std::vector<double> poly_durations = nodes->ConvertPhaseToPolyDurations(params_.ee_phase_durations_.at(ee));
  }

  return vars;
}

std::vector<NodesVariablesPhaseBased::Ptr>
NlpFormulation::MakeForceVariables () const
{
  std::vector<NodesVariablesPhaseBased::Ptr> vars;

  // force is still constant and zero on swing phase!!
  int force_polynomials_per_swing_phase_ = 1;
  std::vector<int> n_polys;

  double T = params_.GetTotalTime();
  for (int ee=0; ee<params_.GetEECount(); ee++) {

	bool phase_constant = !params_.ee_in_contact_at_start_.at(ee);
	const std::vector<double> phase_durations = params_.ee_phase_durations_.at(ee);
	for (auto v : phase_durations) {
		if (phase_constant) {
			n_polys.push_back(1);
		}
		else {
			n_polys.push_back(floor(v/params_.dt_drive_constraint_)+1);
		}
		phase_constant = !phase_constant;
	}

    auto nodes = std::make_shared<NodesVariablesEEForce>(
                                              params_.GetPhaseCount(ee),
                                              params_.ee_in_contact_at_start_.at(ee),
                                              id::EEForceNodes(ee),
											  n_polys);

    // initialize with mass of robot distributed equally on all legs
    double m = model_.dynamic_model_->m();
    double g = model_.dynamic_model_->g();

    n_polys.clear();
    std::vector<double> poly_durations = nodes->ConvertPhaseToPolyDurations(params_.ee_phase_durations_.at(ee));

    Vector3d f_stance(0.0, 0.0, m*g/params_.GetEECount());
    nodes->SetByLinearInterpolation(f_stance, f_stance, T); // stay constant
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
    case Parameters::BaseAccLimits:  			return MakeBaseAccLimitsConstraint(s);
    case Parameters::EndeffectorAcc: 			return MakeEENodesAccConstraint(s);
    case Parameters::EEAccLimits:	 			return MakeEEAccLimitsConstraint(s);
    case Parameters::WheelsLateralConstraint: 	return MakeWheelsLateralConstraint(s);
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
	auto c = std::make_shared<TerrainConstraint>(terrain_, id::EEMotionNodes(ee), params_.min_distance_above_terrain_);
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

  std::vector<int> dim = {X, Y, Z};
  constraints.push_back(std::make_shared<SplineAccConstraint>
                        (s.base_linear_, id::base_lin_nodes, dim));

  constraints.push_back(std::make_shared<SplineAccConstraint>(s.base_angular_, id::base_ang_nodes, dim));

  return constraints;
}

NlpFormulation::ConstraintPtrVec
NlpFormulation::MakeBaseAccLimitsConstraint(const SplineHolder& s) const
{
  std::vector<Vector3d> max_base_acc;
  max_base_acc.push_back(Vector3d(params_.max_base_acc_lin_.data()));  // linear acc limits
  max_base_acc.push_back(Vector3d(params_.max_base_acc_ang_.data()));  // angular acc limits

  return {std::make_shared<BaseAccLimitsConstraint>(max_base_acc,
    											    params_.GetTotalTime(),
												    params_.dt_constraint_dynamic_, s)};
}

NlpFormulation::ConstraintPtrVec
NlpFormulation::MakeEENodesAccConstraint (const SplineHolder& s) const
{
  ConstraintPtrVec constraints;

  std::vector<int> dim = {X};
  for (int ee=0; ee<params_.GetEECount(); ee++) {

//    auto ee_acc = std::make_shared<EEAccConstraint>(s.ee_motion_.at(ee), id::EEMotionNodes(ee), dim);
    auto ee_acc = std::make_shared<SplineAccConstraint>(s.ee_motion_.at(ee), id::EEMotionNodes(ee), dim);

    constraints.push_back(ee_acc);
  }

  return constraints;
}

NlpFormulation::ConstraintPtrVec
NlpFormulation::MakeEEAccLimitsConstraint (const SplineHolder& s) const
{
  ConstraintPtrVec c;

  Vector3d acc_limits (params_.max_wheels_acc_.at(X), params_.max_wheels_acc_.at(Y), params_.max_wheels_acc_.at(Z));
  for (int ee=0; ee<params_.GetEECount(); ee++) {
    auto constraint = std::make_shared<EEAccLimitsConstraint>(acc_limits, ee, s);
    c.push_back(constraint);
  }
  return c;
}

NlpFormulation::ConstraintPtrVec
NlpFormulation::MakeWheelsLateralConstraint (const SplineHolder& s) const
{
  ConstraintPtrVec c;

  for (int ee=0; ee<params_.GetEECount(); ee++) {
    auto constraint = std::make_shared<WheelsLateralConstraint>(terrain_, ee, s);
    c.push_back(constraint);
  }
  return c;
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
