/*
 * force_discretized.cc
 *
 *  Created on: jun, 2020
 *  Author: cla
 */

#ifndef TOWR_INCLUDE_TOWR_CONSTRAINTS_FORCE_CONSTRAINT_DISCRETIZED_H_
#define TOWR_INCLUDE_TOWR_CONSTRAINTS_FORCE_CONSTRAINT_DISCRETIZED_H_

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
class ForceConstraintDiscretized : public TimeDiscretizationConstraint {
public:
  using EE = uint;
  using Vector3d = Eigen::Vector3d;
  using Jacobian = Eigen::SparseMatrix<double, Eigen::RowMajor>;

    ForceConstraintDiscretized (const HeightMap::Ptr& terrain, double T, double dt, const EE& ee,
						        const SplineHolder& spline_holder, double force_limit);
  virtual ~ForceConstraintDiscretized () = default;

private:
  NodeSpline::Ptr ee_wheels_motion_;    ///< the linear position of the wheels.
  NodeSpline::Ptr ee_force_;  ///< the current xyz foot forces.
  HeightMap::Ptr terrain_;    			///< the height map of the current terrain.
  EE ee_;
  int n_constraints_per_node_; 		  	///< number of constraint for each node.
  NodeSpline::Ptr decision_;    ///< the contact state of the wheels.
  double mu_;              ///< friction coeff between robot feet and terrain.
  double fn_max_;          ///< force limit in normal direction.
  void UpdateConstraintAtInstance(double t, int k, VectorXd& g) const override;
  void UpdateBoundsAtInstance(double t, int k, VecBound& bounds) const override;
  void UpdateJacobianAtInstance(double t, int k, std::string var_set, Jacobian& jac) const override;
};
} /* namespace towr */

#endif /* TOWR_INCLUDE_TOWR_CONSTRAINTS_WHEELS_NON_HOLONOMIC_CONSTRAINT_H_ */
