/*
 * wheels_non_holonomic_constraint.h
 *
 *  Created on: Apr 3, 2019
 *      Author: vivian
 */

#ifndef TOWR_INCLUDE_TOWR_CONSTRAINTS_WHEELS_NON_HOLONOMIC_CONSTRAINT_H_
#define TOWR_INCLUDE_TOWR_CONSTRAINTS_WHEELS_NON_HOLONOMIC_CONSTRAINT_H_

#include <towr/variables/spline.h>
#include <towr/variables/spline_holder.h>
#include <towr/variables/cartesian_dimensions.h>
#include <towr/variables/euler_converter.h>
#include <iostream>
#include <towr/variables/variable_names.h>

#include <towr/terrain/height_map.h>

#include "time_discretization_constraint.h"

namespace towr {

// Implementation for pure driving motion: constrain zero velocity in the lateral direction (y)
class WheelsNonHolonomicConstraint : public TimeDiscretizationConstraint {
public:
  using EE = uint;
  using Vector3d = Eigen::Vector3d;
  using Jacobian = Eigen::SparseMatrix<double, Eigen::RowMajor>;

  WheelsNonHolonomicConstraint (const HeightMap::Ptr& terrain, double T, double dt, const EE& ee,
						        const SplineHolder& spline_holder);
  virtual ~WheelsNonHolonomicConstraint () = default;

private:
  NodeSpline::Ptr base_linear_;     	///< the linear position of the base.
  EulerConverter base_angular_; 		///< the orientation of the base.
  NodeSpline::Ptr ee_wheels_motion_;    ///< the linear position of the wheels.
  PhaseDurations::Ptr durations_;
  int n_constraints_per_node_; 		  	///< number of constraint for each node.
  HeightMap::Ptr terrain_;    			///< the height map of the current terrain.
  EE ee_;

  void UpdateConstraintAtInstance(double t, int k, VectorXd& g) const override;
  void UpdateBoundsAtInstance(double t, int k, VecBound& bounds) const override;
  void UpdateJacobianAtInstance(double t, int k, std::string var_set, Jacobian& jac) const override;

  Eigen::Matrix3d GetJacobianTerrainBasis(HeightMap::Direction basis, double x, double y) const;
  int GetCol(int node, int dim) const;
  Jacobian SkewSymmetricMatrix (const Vector3d v) const;

};

} /* namespace towr */



#endif /* TOWR_INCLUDE_TOWR_CONSTRAINTS_WHEELS_NON_HOLONOMIC_CONSTRAINT_H_ */
