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

#ifndef TOWR_NLP_FACTORY_H_
#define TOWR_NLP_FACTORY_H_

#include <ifopt/variable_set.h>
#include <ifopt/constraint_set.h>
#include <ifopt/cost_term.h>

#include <towr/variables/spline_holder.h>
#include <towr/models/robot_model.h>
#include <towr/terrain/height_map.h>
#include <towr/parameters.h>


#include <towr/variables/variable_names.h>
#include <towr/variables/phase_durations.h>

#include <towr/constraints/base_motion_constraint.h>
#include <towr/constraints/dynamic_constraint.h>
#include <towr/constraints/force_constraint.h>
#include <towr/constraints/range_of_motion_constraint.h>
#include <towr/constraints/swing_constraint.h>
#include <towr/constraints/terrain_constraint.h>
#include <towr/constraints/total_duration_constraint.h>
#include <towr/constraints/spline_acc_constraint.h>
#include <towr/constraints/wheels_non_holonomic_constraint.h>
#include <towr/constraints/terrain_constraint_discretized.h>
#include <towr/constraints/force_constraint_discretized.h>



#include <towr/costs/node_cost.h>
#include <towr/variables/nodes_variables_all.h>

#include <iostream>

namespace towr {

/**
 * @defgroup Constraints
 * @brief Constraints of the trajectory optimization problem.
 *
 * These are the constraint sets that characterize legged locomotion.
 *
 * Folder: @ref include/towr/constraints
 */

/**
 * @defgroup Costs
 * @brief Costs of the trajectory optimization problem.
 *
 * These the the cost terms that prioritize certain solutions to the
 * legged locomotion problem.
 *
 * Folder: @ref include/towr/costs
 */

/**
 *
 * @brief A sample combination of variables, cost and constraints.
 *
 * This is _one_ example of how to combine the variables, constraints and costs
 * provided by this library. Additional variables or constraints can be added
 * to the NLP, or existing elements replaced to find a more powerful/general
 * formulation. This formulation was used to generate the motions described
 * in this paper: https://ieeexplore.ieee.org/document/8283570/
 */
class NlpFormulation {
public:
  using VariablePtrVec   = std::vector<ifopt::VariableSet::Ptr>;
  using ConstraintPtrVec = std::vector<ifopt::ConstraintSet::Ptr>;
  using CostPtrVec       = std::vector<ifopt::CostTerm::Ptr>;
  using EEPos            = std::vector<Eigen::Vector3d>;
  using Vector3d         = Eigen::Vector3d;

  NlpFormulation ();
  virtual ~NlpFormulation () = default;

  /**
   * @brief The ifopt variable sets that will be optimized over.
   * @param[in/out] builds fully-constructed splines from the variables.
   */
  VariablePtrVec GetVariableSets(SplineHolder& spline_holder);

  /**
   * @brief The ifopt constraints that enforce feasible motions.
   * @param[in] uses the fully-constructed splines for initialization of constraints.
   */
  ConstraintPtrVec GetConstraints(const SplineHolder& spline_holder) const;

  /** @brief The ifopt costs to tune the motion. */
  ConstraintPtrVec GetCosts() const;


  BaseState initial_base_;
  BaseState final_base_;
  BaseState final_base_v_;
  EEPos  initial_ee_W_;
  RobotModel model_;
  HeightMap::Ptr terrain_;
  Parameters params_;

  std::vector<std::vector<int>> swing_nodes_;

  std::vector<int> endeffector_swing0;
  std::vector<int> endeffector_swing1;
  std::vector<int> endeffector_swing2;
  std::vector<int> endeffector_swing3;

private:
  ConstraintPtrVec MakeWheelsNonHolonomicConstraint(const SplineHolder& s) const;
  ConstraintPtrVec MakeDiscretizedTerrainConstraint(const SplineHolder& s) const;
  ConstraintPtrVec MakeDiscretizedForceConstraint(const SplineHolder& s) const;

  // variables
  std::vector<NodesVariables::Ptr> MakeBaseVariables() const;
  std::vector<NodesVariablesPhaseBased::Ptr> MakeEndeffectorVariables();
  std::vector<NodesVariablesPhaseBased::Ptr> MakeForceVariables() const;
  std::vector<NodesVariablesPhaseBased::Ptr> MakeDecisionVariables() const;
  std::vector<PhaseDurations::Ptr> MakeContactScheduleVariables() const;

  // constraints
  ConstraintPtrVec GetConstraint(Parameters::ConstraintName name,
                                const SplineHolder& splines) const;
  ConstraintPtrVec MakeDynamicConstraint(const SplineHolder& s) const;
  ConstraintPtrVec MakeRangeOfMotionBoxConstraint(const SplineHolder& s) const;
  ConstraintPtrVec MakeTotalTimeConstraint() const;
  ConstraintPtrVec MakeTerrainConstraint() const;
  ConstraintPtrVec MakeForceConstraint() const;
  ConstraintPtrVec MakeSwingConstraint() const;
  ConstraintPtrVec MakeBaseRangeOfMotionConstraint(const SplineHolder& s) const;
  ConstraintPtrVec MakeBaseAccConstraint(const SplineHolder& s) const;

  // costs
  CostPtrVec GetCost(const Parameters::CostName& id, double weight) const;
  CostPtrVec MakeForcesCost(double weight) const;
  CostPtrVec MakeEEMotionCost(double weight) const;
};

} /* namespace towr */

#endif /* TOWR_NLP_FACTORY_H_ */
