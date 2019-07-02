/*
 * nlp_formulation_drive.h
 *
 *  Created on: Mar 28, 2019
 *      Author: vivian
 */

#ifndef TOWR_INCLUDE_TOWR_NLP_FORMULATION_DRIVE_H_
#define TOWR_INCLUDE_TOWR_NLP_FORMULATION_DRIVE_H_

#include <memory>
#include <iostream>

#include <ifopt/variable_set.h>
#include <ifopt/constraint_set.h>
#include <ifopt/cost_term.h>

#include <towr/nlp_formulation.h>
#include <towr/variables/spline_holder_drive.h>
#include <towr/parameters_drive.h>
#include <towr/models/endeffector_mappings.h>

namespace towr {

const static std::map<int, std::string> ee_names =
{ {LF, "LF"}, {RF, "RF"}, {LH, "LH"}, {RH, "RH"}, };

class NlpFormulationDrive : public NlpFormulation
{
public:
  using Base = NlpFormulation;
  using ConstraintPtrVec  = std::vector<ifopt::ConstraintSet::Ptr>;

  NlpFormulationDrive() = default;

  // variables
  VariablePtrVec GetVariableSets(SplineHolderDrive& spline_holder);

  // constraints
  ConstraintPtrVec GetConstraints(const SplineHolderDrive& spline_holder) const;

  // costs
  ConstraintPtrVec GetCosts() const;

  ParametersDrive params_drive_;

  void PrintSolution (const SplineHolderDrive& spline_holder, double dt) const;

private:

  // variables
  std::vector<NodesVariables::Ptr> MakeBaseVariables () const;
  std::vector<NodesVariables::Ptr> MakeEEWheelsMotionVariables () const;
  std::vector<NodesVariables::Ptr> MakeEEWheelsForceVariables () const;

  // constraints
  ConstraintPtrVec GetConstraint (Parameters::ConstraintName name, const SplineHolderDrive& s) const;
  ConstraintPtrVec MakeRangeOfMotionBoxConstraint (const SplineHolderDrive& s) const;
  ConstraintPtrVec MakeBaseRangeOfMotionConstraint (const SplineHolderDrive& s) const;
  ConstraintPtrVec MakeBaseAccConstraint (const SplineHolderDrive& s) const;
  ConstraintPtrVec MakeEndEffectorAccConstraint (const SplineHolderDrive& s) const;
  ConstraintPtrVec MakeDynamicWheelsConstraint(const SplineHolderDrive& s) const;
  ConstraintPtrVec MakeTerrainWheelsConstraint(const SplineHolderDrive& s) const;
  ConstraintPtrVec MakeForceWheelsConstraint(const SplineHolderDrive& s) const;
  ConstraintPtrVec MakeWheelsAccLimitsConstraint(const SplineHolderDrive& s) const;
  ConstraintPtrVec MakeBaseAccLimitsConstraint(const SplineHolderDrive& s) const;
  ConstraintPtrVec MakeWheelsNonHolonomicConstraint(const SplineHolderDrive& s) const;
  ConstraintPtrVec MakeWheelsMotionConstraint (const SplineHolderDrive& s) const;
  ConstraintPtrVec MakeStabilityConstraint (const SplineHolderDrive& s) const;

  // costs
  CostPtrVec GetCost(const Parameters::CostName& id, double weight) const;
  CostPtrVec MakeTorqueCost(double weight) const;

  // for printing the stability angles
  void UpdateDynamicsModel (const SplineHolderDrive& solution, double t) const;

};

}

#endif /* TOWR_INCLUDE_TOWR_NLP_FORMULATION_DRIVE_H_ */
