/*
 * drift_constraint.h
 *
 *  Created on: Apr 25, 2019
 *      Author: dominic landolf
 */


#ifndef TOWR_TOWR_INCLUDE_TOWR_CONSTRAINTS_DRIFT_CONSTRAINT_H_
#define TOWR_TOWR_INCLUDE_TOWR_CONSTRAINTS_DRIFT_CONSTRAINT_H_

#include <towr/variables/spline.h>
#include <towr/variables/spline_holder.h>
#include <towr/variables/euler_converter.h>

#include <towr/models/kinematic_model.h>

#include "time_discretization_constraint.h"
#include <towr/parameters.h>
#include <Eigen/Dense>
using namespace Eigen;

namespace towr {

// constrains the tangential force and the total velocity to point in opposite directions during drift.


class DriftConstraint : public TimeDiscretizationConstraint {
public:
  using EE = uint;
  using Vector3d = Eigen::Vector3d;

  DriftConstraint(const KinematicModel::Ptr& robot_model,
                          double T, double dt,
                          const EE& ee,
                          const SplineHolder& spline_holder);
  virtual ~DriftConstraint() = default;

private:
  NodeSpline::Ptr base_linear_;     ///< the linear position of the base.
  EulerConverter base_angular_; 	///< the orientation of the base.
  NodeSpline::Ptr base_angular_for_angles_;
  NodeSpline::Ptr ee_motion_;       ///< the linear position of the endeffectors.
  NodeSpline::Ptr ee_force_;		///< forces at the endeffectors
  Parameters params_;

  int n_constraints_per_node_;

  Vector3d NormDerivative (const Vector3d& v) const;

  EE ee_;

  // see TimeDiscretizationConstraint for documentation
  void UpdateConstraintAtInstance (double t, int k, VectorXd& g) const override;
  void UpdateBoundsAtInstance (double t, int k, VecBound&) const override;
  void UpdateJacobianAtInstance(double t, int k, std::string, Jacobian&) const override;

  VecBound node_bounds_;     ///< same bounds for each discretized node
  int GetRow(int node, int dimension) const;
};

} /* namespace towr */

#endif /* TOWR_TOWR_INCLUDE_TOWR_CONSTRAINTS_DRIVE_CONSTRAINT_H_ */
