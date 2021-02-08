/*
 * base_acc_limits_constraint.h
 *
 *  Created on: Apr 3, 2019
 *      Author: vivian
 */

#ifndef TOWR_CONSTRAINTS_BASE_ACC_LIMITS_CONSTRAINT_H_
#define TOWR_CONSTRAINTS_BASE_ACC_LIMITS_CONSTRAINT_H_

#include <towr/variables/spline_holder.h>
#include <towr/variables/spline.h>
#include <towr/variables/cartesian_dimensions.h>

#include "time_discretization_constraint.h"

namespace towr {

// Limits the acceleration of the base
class BaseAccLimitsConstraint : public TimeDiscretizationConstraint {
public:
  using Vector3d = Eigen::Vector3d;
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
  BaseAccLimitsConstraint (const std::vector<Vector3d>& acc_max,
		  	  	  	  	   double T, double dt,
						   const SplineHolder& spline_holder);
  virtual ~BaseAccLimitsConstraint () = default;

private:
  NodeSpline::Ptr base_linear_;
  NodeSpline::Ptr base_angular_;
  Vector3d acc_max_lin_;
  Vector3d acc_max_ang_;

  int GetRow (int node, int dim) const;

  void UpdateConstraintAtInstance(double t, int k, VectorXd& g) const override;
  void UpdateBoundsAtInstance(double t, int k, VecBound& bounds) const override;
  void UpdateJacobianAtInstance(double t, int k, std::string var_set, Jacobian& jac) const override;

};

} /* namespace towr */



#endif /* TOWR_INCLUDE_TOWR_CONSTRAINTS_BASE_ACC_CONSTRAINT_H_ */
