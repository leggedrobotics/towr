/*
 * ee_acc_limits_constraint.h
 *
 *  Created on: Jan 17, 2020
 *      Author: vivian
 */

#ifndef TOWR_INCLUDE_TOWR_CONSTRAINTS_EE_ACC_LIMITS_CONSTRAINT_H_
#define TOWR_INCLUDE_TOWR_CONSTRAINTS_EE_ACC_LIMITS_CONSTRAINT_H_

#include <towr/variables/spline.h>
#include <towr/variables/spline_holder.h>
#include <towr/variables/cartesian_dimensions.h>

#include "time_discretization_constraint.h"

namespace towr {

// Limits the acceleration on the wheels
class EEAccLimitsConstraint : public ifopt::ConstraintSet {
public:
  using EE = uint;
  using Vector3d = Eigen::Vector3d;
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

  EEAccLimitsConstraint (const Vector3d& acc_max, EE ee,
						 const SplineHolder& spline_holder);
  virtual ~EEAccLimitsConstraint () = default;

  VectorXd GetValues() const override;
  VecBound GetBounds() const override;
  void FillJacobianBlock (std::string var_set, Jacobian&) const override;

private:
  NodeSpline::Ptr ee_motion_; ///< endeffector position in world frame.
  Vector3d acc_max_;
  EE ee_;

  int n_polys_;       ///< number of polynomials in spline.
  std::vector<double> T_; ///< Duration of each polynomial in spline.

  int GetRow (int node, int dim) const;

};

} /* namespace towr */



#endif /* TOWR_INCLUDE_TOWR_CONSTRAINTS_EE_ACC_LIMITS_CONSTRAINT_H_ */
