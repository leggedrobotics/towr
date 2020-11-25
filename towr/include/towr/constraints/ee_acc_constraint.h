/*
 * ee_acc_constraint.h
 *
 *  Created on: Jan 16, 2020
 *      Author: vivian
 */

#ifndef TOWR_CONSTRAINTS_EE_ACC_CONSTRAINT_H_
#define TOWR_CONSTRAINTS_EE_ACC_CONSTRAINT_H_

#include <ifopt/constraint_set.h>

#include <towr/variables/node_spline.h>
#include <towr/variables/nodes_variables_phase_based.h>

namespace towr {

/**
 *  @brief Ensures continuous accelerations between ee polynomials during stance phase.
 *
 * @ingroup Constraints
 */
class EEAccConstraint : public ifopt::ConstraintSet {
public:
  EEAccConstraint(const NodeSpline::Ptr& spline, std::string name, const std::vector<int>& dimensions);
  virtual ~EEAccConstraint() = default;

  void InitVariableDependedQuantities(const VariablesPtr& x) override;

  VectorXd GetValues() const override;
  VecBound GetBounds() const override;
  void FillJacobianBlock (std::string var_set, Jacobian&) const override;

private:
  NodeSpline::Ptr spline_;        ///< a spline comprised of polynomials
  std::string ee_motion_id_; /// polynomial parameterized node values

  NodesVariablesPhaseBased::Ptr ee_motion_; ///< the position of the endeffector at each node.
  std::vector<int> stance_node_ids_; ///< the indices of nodes in stance phase.
  std::vector<int> stance_poly_ids_; ///< the indices of the polynomials in stance phase.

  std::vector<int> swing_node_ids_; ///< the indices of nodes in stance phase.
  std::vector<int> swing_poly_ids_; ///< the indices of the polynomials in stance phase.

  int n_junctions_;       ///< number of polynomials junctions in spline.
  std::vector<double> T_; ///< Duration of each polynomial in spline.

  std::vector<int> dimensions_;  ///< dimensions that are constrained
  int n_dim_;             ///< number of constrained dimensions
};

} /* namespace towr */


#endif /* TOWR_CONSTRAINTS_EE_ACC_CONSTRAINT_H_ */
