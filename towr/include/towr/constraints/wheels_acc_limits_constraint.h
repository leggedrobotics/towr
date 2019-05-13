/*
 * wheels_acc_limits_constraint.h
 *
 *  Created on: Apr 2, 2019
 *      Author: vivian
 */

#ifndef TOWR_INCLUDE_TOWR_CONSTRAINTS_WHEELS_ACC_LIMITS_CONSTRAINT_H_
#define TOWR_INCLUDE_TOWR_CONSTRAINTS_WHEELS_ACC_LIMITS_CONSTRAINT_H_

#include <towr/variables/spline.h>
#include <towr/variables/spline_holder_drive.h>
#include <towr/variables/cartesian_dimensions.h>

#include "time_discretization_constraint.h"

namespace towr {

// Implementation for pure driving motion: limits the acceleration on the wheels
class WheelsAccLimitsConstraint : public TimeDiscretizationConstraint {
public:
  using EE = uint;
  using Vector3d = Eigen::Vector3d;

  WheelsAccLimitsConstraint (Vector3d acc_max, double T, double dt, const EE& ee,
						     const SplineHolderDrive& spline_holder);
  virtual ~WheelsAccLimitsConstraint () = default;

private:
  NodeSpline::Ptr ee_wheels_motion_; ///< endeffector position in world frame.
  Vector3d acc_max_;
  EE ee_;

  int GetRow (int node, int dim) const;

  void UpdateConstraintAtInstance(double t, int k, VectorXd& g) const override;
  void UpdateBoundsAtInstance(double t, int k, VecBound& bounds) const override;
  void UpdateJacobianAtInstance(double t, int k, std::string var_set, Jacobian& jac) const override;

};

} /* namespace towr */



#endif /* TOWR_INCLUDE_TOWR_CONSTRAINTS_WHEELS_ACC_LIMITS_CONSTRAINT_H_ */
