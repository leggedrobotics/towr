/*
 * wheels_lateral_constraint.h
 *
 *  Created on: Feb 13, 2020
 *      Author: vivian
 */

#ifndef TOWR_CONSTRAINTS_WHEELS_LATERAL_CONSTRAINT_H_
#define TOWR_CONSTRAINTS_WHEELS_LATERAL_CONSTRAINT_H_

#include <ifopt/constraint_set.h>

#include <towr/variables/spline.h>
#include <towr/variables/spline_holder.h>
#include <towr/variables/euler_converter.h>

#include <towr/variables/nodes_variables_phase_based.h>
#include <towr/terrain/height_map.h>

namespace towr {

class WheelsLateralConstraint : public ifopt::ConstraintSet {
public:
  using EE = uint;
  using Vector3d = Eigen::Vector3d;
  using Jacobian = Eigen::SparseMatrix<double, Eigen::RowMajor>;
  using VecDurations = std::vector<double>;

  WheelsLateralConstraint (const HeightMap::Ptr& terrain, const EE& ee,
		  	  	  	  	   const SplineHolder& spline_holder);
  virtual ~WheelsLateralConstraint () = default;

  void InitVariableDependedQuantities(const VariablesPtr& x) override;

  VectorXd GetValues() const override;
  VecBound GetBounds() const override;
  void FillJacobianBlock (std::string var_set, Jacobian&) const override;

private:

  std::vector<double> GetVectorGlobalTimeAtNodes (VecDurations poly);
  Eigen::Matrix3d GetJacobianTerrainBasis(HeightMap::Direction basis, double x, double y) const;
  Jacobian SkewSymmetricMatrix (const Vector3d v) const;

  NodesVariablesPhaseBased::Ptr ee_motion_; ///< the position of the endeffector (nodes).

  NodeSpline::Ptr base_linear_;     ///< the linear position of the base.
  EulerConverter base_angular_; ///< the orientation of the base.
  NodeSpline::Ptr ee_spline_; ///the position of the endeffector (spline).

  HeightMap::Ptr terrain_;    ///< the height map of the current terrain.

  EE ee_;
  int n_constraints_per_node_; ///< number of constraint for each node.

  VecDurations poly_durations_;
  std::vector<double> T_;

  std::string ee_motion_id_;  ///< the name of the endeffector variable set.
  std::vector<int> pure_stance_node_ids_; ///< the indices of the nodes constrained.
};

} /* namespace towr */



#endif /* TOWR_CONSTRAINTS_WHEELS_LATERAL_CONSTRAINT_H_ */
