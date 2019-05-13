/*
 * force_wheels_constraint.h
 *
 *  Created on: Apr 1, 2019
 *      Author: vivian
 */

#ifndef TOWR_INCLUDE_TOWR_CONSTRAINTS_FORCE_WHEELS_CONSTRAINT_H_
#define TOWR_INCLUDE_TOWR_CONSTRAINTS_FORCE_WHEELS_CONSTRAINT_H_

#include <towr/variables/spline.h>
#include <towr/variables/spline_holder_drive.h>

#include <towr/terrain/height_map.h>

#include "time_discretization_constraint.h"

namespace towr {

// Implementation for pure driving motion: forces constraints must be ensured at a regular time interval,
// which makes this constraint a time discretization constraint
class ForceWheelsConstraint : public TimeDiscretizationConstraint {
public:
  using EE = uint;
  using Vector3d = Eigen::Vector3d;

  ForceWheelsConstraint (const HeightMap::Ptr& terrain, double ft_max,
		  	  	  	  	 double T, double dt, const EE& ee,
						 const SplineHolderDrive& spline_holder);
  virtual ~ForceWheelsConstraint () = default;

private:
  NodeSpline::Ptr ee_wheels_motion_;  ///< endeffector position in world frame.
  NodeSpline::Ptr ee_wheels_force_;   ///< endeffector forces in world frame.
  HeightMap::Ptr terrain_;    		  ///< the height map of the current terrain.
  double ft_max_;          			  ///< force limit in tangent direction.
  double mu_;              			  ///< friction coeff between robot feet and terrain.
  int n_constraints_per_node_; 		  ///< number of constraint for each node.
  EE ee_;                  			  ///< The endeffector force to be constrained.

  void UpdateConstraintAtInstance(double t, int k, VectorXd& g) const override;
  void UpdateBoundsAtInstance(double t, int k, VecBound& bounds) const override;
  void UpdateJacobianAtInstance(double t, int k, std::string var_set, Jacobian& jac) const override;

  Eigen::Matrix3d GetJacobianTerrainBasis(HeightMap::Direction basis, double x, double y) const;
  int GetCol(int node, int dim) const;

};

} /* namespace towr */



#endif /* TOWR_INCLUDE_TOWR_CONSTRAINTS_FORCE_WHEELS_CONSTRAINT_H_ */
