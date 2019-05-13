/*
 * terrain_wheels_constraint.h
 *
 *  Created on: Mar 29, 2019
 *      Author: vivian
 */

#ifndef TOWR_INCLUDE_TOWR_CONSTRAINTS_TERRAIN_WHEELS_CONSTRAINT_H_
#define TOWR_INCLUDE_TOWR_CONSTRAINTS_TERRAIN_WHEELS_CONSTRAINT_H_

#include <towr/variables/spline.h>
#include <towr/variables/spline_holder_drive.h>

#include <towr/terrain/height_map.h>

#include "time_discretization_constraint.h"

namespace towr {

// Implementation for pure driving motion: contact must be ensured at regular time interval,
// which makes this constraint a time discretization constraint
class TerrainWheelsConstraint : public TimeDiscretizationConstraint {
public:
  using EE = uint;
  using Vector3d = Eigen::Vector3d;

  TerrainWheelsConstraint (const HeightMap::Ptr& terrain,
		  	  	  	  	   double T, double dt, const EE& ee,
						   const SplineHolderDrive& spline_holder);
  virtual ~TerrainWheelsConstraint () = default;

private:
  NodeSpline::Ptr ee_wheels_motion_; ///< endeffector position in world frame.
  HeightMap::Ptr terrain_;    ///< the height map of the current terrain.
  EE ee_;

  void UpdateConstraintAtInstance(double t, int k, VectorXd& g) const override;
  void UpdateBoundsAtInstance(double t, int k, VecBound& bounds) const override;
  void UpdateJacobianAtInstance(double t, int k, std::string var_set, Jacobian& jac) const override;

  Jacobian GetJacobianTerrainHeight(double x, double y) const;

};

} /* namespace towr */

#endif /* TOWR_INCLUDE_TOWR_CONSTRAINTS_TERRAIN_WHEELS_CONSTRAINT_H_ */
