/*
 * wheels_motion_constraint.h
 *
 *  Created on: Apr 4, 2019
 *      Author: vivian
 */

#ifndef TOWR_INCLUDE_TOWR_CONSTRAINTS_WHEELS_MOTION_CONSTRAINT_H_
#define TOWR_INCLUDE_TOWR_CONSTRAINTS_WHEELS_MOTION_CONSTRAINT_H_

#include <towr/variables/spline_holder_drive.h>
#include <towr/variables/spline.h>
#include <towr/variables/euler_converter.h>

#include <towr/models/kinematic_model.h>

#include "time_discretization_constraint.h"

namespace towr {

class WheelsMotionConstraint : public TimeDiscretizationConstraint {
public:
  using EE = uint;
  using Vector3d = Eigen::Vector3d;

  WheelsMotionConstraint (const KinematicModel::Ptr& model,
		    			  double T, double dt, const EE& ee,
						  const SplineHolderDrive& spline_holder);
  virtual ~WheelsMotionConstraint () = default;

  void UpdateConstraintAtInstance (double t, int k, VectorXd& g) const override;
  void UpdateBoundsAtInstance (double t, int k, VecBound&) const override;
  void UpdateJacobianAtInstance(double t, int k, std::string, Jacobian&) const override;

private:
  NodeSpline::Ptr base_linear_;     ///< the linear position of the base.
  EulerConverter base_angular_; 	///< the orientation of the base.
  NodeSpline::Ptr ee_wheels_motion_;       ///< the linear position of the wheels.

  double max_leg_extension_;
  double min_leg_extension_;
  Vector3d nominal_pos_hip_B_;
  Vector3d max_deviation_from_nominal_;
  Vector3d nominal_ee_pos_B_;
  int n_constraints_per_node_;
  EE ee_;

  Vector3d NormDerivative (const Vector3d& v) const;
  int GetCol(int node, int dim) const;
};

} /* namespace towr */



#endif /* TOWR_INCLUDE_TOWR_CONSTRAINTS_WHEELS_MOTION_CONSTRAINT_H_ */
