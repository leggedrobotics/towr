/*
 * nlp_formulation_extended.cc
 *
 *  Created on: Mar 28, 2019
 *      Author: vivian
 */
#include <iostream>
#include <cmath>

#include "towr/nlp_formulation_drive.h"

#include <towr/variables/variable_names.h>
#include <towr/variables/nodes_variables_all.h>
#include <towr/variables/nodes_variables_ee_wheels.h>

#include <towr/constraints/base_motion_constraint.h>
#include <towr/constraints/range_of_motion_constraint.h>
#include <towr/constraints/spline_acc_constraint.h>
#include <towr/constraints/dynamic_constraint.h>
#include <towr/constraints/terrain_wheels_constraint.h>
#include <towr/constraints/force_wheels_constraint.h>
#include <towr/constraints/wheels_acc_limits_constraint.h>
#include <towr/constraints/base_acc_limits_constraint.h>
#include <towr/constraints/wheels_motion_constraint.h>
#include <towr/constraints/stability_constraint.h>
#include <towr/constraints/wheels_non_holonomic_constraint.h>

namespace towr {

NlpFormulationDrive::VariablePtrVec
NlpFormulationDrive::GetVariableSets (SplineHolderDrive& spline_holder)
{
  VariablePtrVec vars;

  auto base_motion = MakeBaseVariables();
  vars.insert(vars.end(), base_motion.begin(), base_motion.end());

  auto ee_wheels_motion = MakeEEWheelsMotionVariables();
  vars.insert(vars.end(), ee_wheels_motion.begin(), ee_wheels_motion.end());

  auto ee_wheels_force = MakeEEWheelsForceVariables();
  vars.insert(vars.end(), ee_wheels_force.begin(), ee_wheels_force.end());

  // stores these readily constructed spline
  spline_holder = SplineHolderDrive(base_motion.at(0), // linear
                               	    base_motion.at(1), // angular
									params_drive_.GetBasePolyDurations(),
									ee_wheels_motion,
									ee_wheels_force,
									params_drive_.GetEEWheelsPolyDurations());

  return vars;
}

std::vector<NodesVariables::Ptr>
NlpFormulationDrive::MakeBaseVariables () const
{
  std::vector<NodesVariables::Ptr> vars;

  int n_nodes = params_drive_.GetBasePolyDurations().size() + 1;

  auto spline_lin = std::make_shared<NodesVariablesAll>(n_nodes, k3D, id::base_lin_nodes);

  double x = final_base_.lin.p().x();
  double y = final_base_.lin.p().y();
  double z = terrain_->GetHeight(x,y) - model_.kinematic_model_->GetNominalStanceInBase().front().z();
  Vector3d final_pos(x, y, z);

  spline_lin->SetByLinearInterpolation(initial_base_.lin.p(), final_pos, params_drive_.GetTotalTime());
  spline_lin->AddStartBound(kPos, {X,Y,Z}, initial_base_.lin.p());
  spline_lin->AddStartBound(kVel, {X,Y,Z}, initial_base_.lin.v());
  spline_lin->AddFinalBound(kPos, params_drive_.bounds_final_lin_pos_, final_base_.lin.p());
  spline_lin->AddFinalBound(kVel, params_drive_.bounds_final_lin_vel_, final_base_.lin.v());
  vars.push_back(spline_lin);

  auto spline_ang = std::make_shared<NodesVariablesAll>(n_nodes, k3D, id::base_ang_nodes);
  spline_ang->SetByLinearInterpolation(initial_base_.ang.p(), final_base_.ang.p(), params_drive_.GetTotalTime());

  // limits all angles displacement to 30 degrees (optional)
  Vector3d ang_limit_ = Vector3d(30, 30, 30) * (M_PI/180);
//  spline_ang->AddAllNodesBounds(kPos, {X,Y,Z}, -ang_limit_, ang_limit_);

  spline_ang->AddStartBound(kPos, {X,Y,Z}, initial_base_.ang.p());
  spline_ang->AddStartBound(kVel, {X,Y,Z}, initial_base_.ang.v());
  spline_ang->AddFinalBound(kPos, params_drive_.bounds_final_ang_pos_, final_base_.ang.p());
  spline_ang->AddFinalBound(kVel, params_drive_.bounds_final_ang_vel_, final_base_.ang.v());
  vars.push_back(spline_ang);

  return vars;
}

std::vector<NodesVariables::Ptr>
NlpFormulationDrive::MakeEEWheelsMotionVariables () const
{
  std::vector<NodesVariables::Ptr> vars;

  int n_nodes = params_drive_.GetEEWheelsPolyDurations().size() + 1;
//  std::cout<< n_nodes << std::endl;

  // Endeffector Motions
  double T = params_drive_.GetTotalTime();
  for (int ee=0; ee<params_drive_.GetEECount(); ee++) {
    auto nodes = std::make_shared<NodesVariablesEEWheelsMotion>(n_nodes, k3D, id::EEWheelsMotionNodes(ee), terrain_);

    // initialize towards final footholds
    double yaw = final_base_.ang.p().z();
    Eigen::Vector3d euler(0.0, 0.0, yaw);
    Eigen::Matrix3d w_R_b = EulerConverter::GetRotationMatrixBaseToWorld(euler);
    Vector3d final_ee_pos_W = final_base_.lin.p() + w_R_b*model_.kinematic_model_->GetNominalStanceInBase().at(ee);
    double x = final_ee_pos_W.x();
    double y = final_ee_pos_W.y();
    double z = terrain_->GetHeight(x,y);
    nodes->SetByLinearInterpolation(initial_ee_W_.at(ee), Vector3d(x,y,z), T);

    // initial and final wheel's speed = 0 (optional)
//    nodes->AddStartBound(kVel, {X,Y,Z}, Vector3d(0.0, 0.0, 0.0));  // wheels vel zero at the beginning
//    nodes->AddFinalBound(kVel, {X,Y,Z}, Vector3d(0.0, 0.0, 0.0));  // wheels vel zero at the end

    // initial wheel's y-position
    nodes->AddStartBound(kPos, {Y}, initial_ee_W_.at(ee));

    // wheels non-holonomic constraint (v_y = 0) -> IMPORTANT!
    if (!params_drive_.use_non_holonomic_constraint_)
      nodes->AddAllNodesBounds(kVel, {Y}, Vector3d(0.0, 0.0, 0.0), Vector3d(0.0, 0.0, 0.0));

    vars.push_back(nodes);
  }

  return vars;
}

std::vector<NodesVariables::Ptr>
NlpFormulationDrive::MakeEEWheelsForceVariables () const
{
  std::vector<NodesVariables::Ptr> vars;

  int n_nodes = params_drive_.GetEEWheelsPolyDurations().size() + 1;

  for (int ee=0; ee<params_drive_.GetEECount(); ee++) {
    auto nodes = std::make_shared<NodesVariablesEEWheelsForce>(n_nodes, k3D, id::EEWheelsForceNodes(ee));

    double f_tangent = params_drive_.force_limit_in_x_direction_;
    double f_normal = (model_.dynamic_model_->m()*model_.dynamic_model_->g())/params_drive_.GetEECount();
    nodes->SetAllAtConstantValue(Vector3d(f_tangent, 0.0, f_normal));

    vars.push_back(nodes);
  }

  return vars;
}

NlpFormulationDrive::ConstraintPtrVec
NlpFormulationDrive::GetConstraints(const SplineHolderDrive& spline_holder) const
{
  ContraintPtrVec constraints;
  for (auto name : params_drive_.constraints_)
    for (auto c : GetConstraint(name, spline_holder))
      constraints.push_back(c);

  return constraints;
}

NlpFormulationDrive::ConstraintPtrVec
NlpFormulationDrive::GetConstraint (Parameters::ConstraintName name,
                           	   	    const SplineHolderDrive& s) const
{
  switch (name) {
    case Parameters::DynamicWheels:     	return MakeDynamicWheelsConstraint(s);
    case Parameters::EndeffectorRom: 		return MakeRangeOfMotionBoxConstraint(s);
    case Parameters::BaseRom:        		return MakeBaseRangeOfMotionConstraint(s);
    case Parameters::TerrainWheels:     	return MakeTerrainWheelsConstraint(s);
    case Parameters::ForceWheels:       	return MakeForceWheelsConstraint(s);
    case Parameters::BaseAcc:      			return MakeBaseAccConstraint(s);
    case Parameters::EndeffectorAcc:  		return MakeEndEffectorAccConstraint(s);
    case Parameters::BaseAccLimits:  		return MakeBaseAccLimitsConstraint(s);
    case Parameters::WheelsAccLimits:  		return MakeWheelsAccLimitsConstraint(s);
    case Parameters::WheelsMotion: 			return MakeWheelsMotionConstraint(s);
    case Parameters::Stability: 			return MakeStabilityConstraint(s);
    case Parameters::WheelsNonHolonomic:	return MakeWheelsNonHolonomicConstraint(s);
    default: throw std::runtime_error("constraint not defined!");
  }
}

NlpFormulationDrive::ConstraintPtrVec
NlpFormulationDrive::MakeRangeOfMotionBoxConstraint (const SplineHolderDrive& s) const
{
  ContraintPtrVec c;

  for (int ee=0; ee<params_drive_.GetEECount(); ee++) {
    auto rom = std::make_shared<RangeOfMotionConstraint>(model_.kinematic_model_,
                                                         params_drive_.GetTotalTime(),
														 params_drive_.dt_constraint_range_of_motion_,
                                                         ee,
                                                         s);
    c.push_back(rom);
  }

  return c;
}

NlpFormulationDrive::ContraintPtrVec
NlpFormulationDrive::MakeBaseRangeOfMotionConstraint (const SplineHolderDrive& s) const
{
  return {std::make_shared<BaseMotionConstraint>(params_drive_.GetTotalTime(),
		  	  	  	  	  	  	  	  	  	  	 params_drive_.dt_constraint_base_motion_,
                                                 s)};
}

NlpFormulationDrive::ContraintPtrVec
NlpFormulationDrive::MakeBaseAccConstraint (const SplineHolderDrive& s) const
{
  ContraintPtrVec constraints;

  constraints.push_back(std::make_shared<SplineAccConstraint>(s.base_linear_, id::base_lin_nodes));

  constraints.push_back(std::make_shared<SplineAccConstraint>(s.base_angular_, id::base_ang_nodes));

  return constraints;
}

NlpFormulationDrive::ConstraintPtrVec
NlpFormulationDrive::MakeEndEffectorAccConstraint (const SplineHolderDrive& s) const
{
  ConstraintPtrVec c;

  for (int ee=0; ee<params_drive_.GetEECount(); ee++) {
    auto constraint = std::make_shared<SplineAccConstraint>(s.ee_wheels_motion_.at(ee), id::EEWheelsMotionNodes(ee));

    c.push_back(constraint);
  }

  return c;
}

NlpFormulationDrive::ContraintPtrVec
NlpFormulationDrive::MakeDynamicWheelsConstraint(const SplineHolderDrive& s) const
{
  auto constraint = std::make_shared<DynamicConstraint>(model_.dynamic_model_,
		  	  	  	  	  	  	  	  	  	  	  	    params_drive_.GetTotalTime(),
														params_drive_.dt_constraint_dynamic_,
                                                        s);
  return {constraint};
}

NlpFormulationDrive::ContraintPtrVec
NlpFormulationDrive::MakeTerrainWheelsConstraint(const SplineHolderDrive& s) const
{
  ContraintPtrVec c;

  for (int ee=0; ee<params_drive_.GetEECount(); ee++) {
    auto constraint = std::make_shared<TerrainWheelsConstraint>(terrain_, params_drive_.GetTotalTime(),
														  	    params_drive_.dt_constraint_dynamic_,
																ee, s);
    c.push_back(constraint);
  }

  return c;
}

NlpFormulationDrive::ContraintPtrVec
NlpFormulationDrive::MakeForceWheelsConstraint(const SplineHolderDrive& s) const
{
  ContraintPtrVec c;

  for (int ee=0; ee<params_drive_.GetEECount(); ee++) {
    auto constraint = std::make_shared<ForceWheelsConstraint>(terrain_, params_drive_.force_limit_in_x_direction_,
    														  params_drive_.GetTotalTime(),
														  	  params_drive_.dt_constraint_dynamic_, ee, s);
    c.push_back(constraint);
  }
  return c;
}

NlpFormulationDrive::ContraintPtrVec
NlpFormulationDrive::MakeWheelsAccLimitsConstraint(const SplineHolderDrive& s) const
{
  ContraintPtrVec c;

  for (int ee=0; ee<params_drive_.GetEECount(); ee++) {
    auto constraint = std::make_shared<WheelsAccLimitsConstraint>(Vector3d(params_drive_.max_wheels_acc_.data()),
    														      params_drive_.GetTotalTime(),
														  	      params_drive_.dt_constraint_dynamic_, ee, s);
    c.push_back(constraint);
  }
  return c;
}

NlpFormulationDrive::ContraintPtrVec
NlpFormulationDrive::MakeBaseAccLimitsConstraint(const SplineHolderDrive& s) const
{
  std::vector<Vector3d> max_base_acc;
  max_base_acc.push_back(Vector3d(params_drive_.max_base_acc_lin_.data()));  // linear acc limits
  max_base_acc.push_back(Vector3d(params_drive_.max_base_acc_ang_.data()));  // angular acc limits

  return {std::make_shared<BaseAccLimitsConstraint>(max_base_acc,
    											    params_drive_.GetTotalTime(),
												    params_drive_.dt_constraint_dynamic_, s)};
}

NlpFormulationDrive::ConstraintPtrVec
NlpFormulationDrive::MakeWheelsMotionConstraint (const SplineHolderDrive& s) const
{
  ContraintPtrVec c;

  for (int ee=0; ee<params_drive_.GetEECount(); ee++) {
    auto constraint = std::make_shared<WheelsMotionConstraint>(model_.kinematic_model_,
                                                        	   params_drive_.GetTotalTime(),
															   params_drive_.dt_constraint_range_of_motion_,
															   ee, s);
    c.push_back(constraint);
  }

  return c;
}

NlpFormulationDrive::ConstraintPtrVec
NlpFormulationDrive::MakeStabilityConstraint (const SplineHolderDrive& s) const
{
  auto constraint = std::make_shared<StabilityConstraint>(model_.dynamic_model_,
		  	  	  	  	  	  	  	  	  	  	  	    params_drive_.GetTotalTime(),
														params_drive_.dt_constraint_dynamic_,
                                                        s);
  return {constraint};
}

NlpFormulationDrive::ConstraintPtrVec
NlpFormulationDrive::MakeWheelsNonHolonomicConstraint (const SplineHolderDrive& s) const
{
  ContraintPtrVec c;

  for (int ee=0; ee<params_drive_.GetEECount(); ee++) {
	auto constraint = std::make_shared<WheelsNonHolonomicConstraint>(terrain_, params_drive_.GetTotalTime(),
																	 params_drive_.dt_constraint_dynamic_, ee, s);
    c.push_back(constraint);
  }

  return c;
}

void
NlpFormulationDrive::UpdateDynamicsModel (const SplineHolderDrive& solution, double t) const
{
  auto com = solution.base_linear_->GetPoint(t);
  EulerConverter base_angular_ = EulerConverter(solution.base_angular_);

  Eigen::Matrix3d w_R_b = base_angular_.GetRotationMatrixBaseToWorld(t);
  Eigen::Vector3d omega = base_angular_.GetAngularVelocityInWorld(t);
  Eigen::Vector3d omega_dot = base_angular_.GetAngularAccelerationInWorld(t);

  int n_ee = model_.dynamic_model_->GetEECount();
  std::vector<Eigen::Vector3d> ee_pos;
  std::vector<Eigen::Vector3d> ee_force;
  for (int ee=0; ee<n_ee; ++ee) {
    ee_force.push_back(solution.ee_wheels_force_.at(ee)->GetPoint(t).p());
    ee_pos.push_back(solution.ee_wheels_motion_.at(ee)->GetPoint(t).p());
  }

  model_.dynamic_model_->SetCurrent(com.p(), com.a(), w_R_b, omega, omega_dot, ee_force, ee_pos);
}

void
NlpFormulationDrive::PrintSolution (const SplineHolderDrive& solution, double dt) const
{
  double t = 0.0;
  while (t<=solution.base_linear_->GetTotalTime() + 1e-5) {
	std::cout << "t=" << t << "\n";
	std::cout << "Base linear position x,y,z:   \t";
	std::cout << solution.base_linear_->GetPoint(t).p().transpose() << "\t[m]" << std::endl;

	std::cout << "Base Euler roll, pitch, yaw:  \t";
	Vector3d rad = solution.base_angular_->GetPoint(t).p();
	std::cout << (rad/M_PI*180).transpose() << "\t[deg]" << std::endl;

	for (int ee=0; ee<solution.ee_wheels_motion_.size(); ee++) {
	  std::string ee_name = ee_names.find(ee)->second;
      std::cout << "Foot position x,y,z (" << ee_name << "):\t";
	  std::cout << solution.ee_wheels_motion_.at(ee)->GetPoint(t).p().transpose() << "\t[m]" << std::endl;
	}

	for (int ee=0; ee<solution.ee_wheels_force_.size(); ee++) {
	  std::string ee_name = ee_names.find(ee)->second;
	  std::cout << "Contact force x,y,z (" << ee_name << "):\t";
	  std::cout << solution.ee_wheels_force_.at(ee)->GetPoint(t).p().transpose() << "\t[N]" << std::endl;
	}

	UpdateDynamicsModel (solution, t);
	std::pair<int, double> stab_angle = model_.dynamic_model_->GetStabilityMeasure();
	std::cout << "Stability Measure: " << stab_angle.second/M_PI*180 << " [deg]" << std::endl;
//	std::vector<double> stab_angle = model_.dynamic_model_->GetStabilityMeasure();
//	std::cout << "Stability Measure (LF,RF,RH,LH): \t";
//	for (int ee=0; ee<solution.ee_wheels_force_.size(); ee++) {
//	  std::cout << stab_angle.at(ee)/M_PI*180 << " ";
//	}
//	std::cout << "\t[deg]" << std::endl;

    std::cout << std::endl;

    t += dt;
  }
}

} /* namespace towr */
