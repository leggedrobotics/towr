/*
 * drive_constraint.h
 *
 *  Created on: Apr 25, 2019
 *      Author: dominic
 */

#ifndef TOWR_TOWR_INCLUDE_TOWR_CONSTRAINTS_DRIVE_CONSTRAINT_H_
#define TOWR_TOWR_INCLUDE_TOWR_CONSTRAINTS_DRIVE_CONSTRAINT_H_

#include <towr/variables/spline.h>
#include <towr/variables/spline_holder.h>
#include <towr/variables/euler_converter.h>

#include <towr/models/kinematic_model.h>

#include "time_discretization_constraint.h"
#include <towr/parameters.h>

namespace towr {

/** @brief Constrains an endeffector to lie in a box around the nominal stance.
  *
  * These constraints are necessary to avoid configurations
  * that are outside the kinematic reach of the robot. The constraint
  * is defined by Cartesian estimates of the reachability of each endeffector.
  *
  * This constraint calculates the position of of the contact expressed in the
  * current CoM frame and constrains it to lie in a box around the nominal/
  * natural contact position for that leg.
  *
  * @ingroup Constraints
  */
class DriveConstraint : public TimeDiscretizationConstraint {
public:
  using EE = uint;
  using Vector3d = Eigen::Vector3d;

  /**
   * @brief Constructs a constraint instance.
   * @param robot_model   The kinematic restrictions of the robot.
   * @param T   The total duration of the optimization.
   * @param dt  the discretization intervall at which to enforce constraints.
   * @param ee            The endeffector for which to constrain the range.
   * @param spline_holder Pointer to the current variables.
   */
  DriveConstraint(const KinematicModel::Ptr& robot_model,
                          double T, double dt,
                          const EE& ee,
                          const SplineHolder& spline_holder);
  virtual ~DriveConstraint() = default;

private:
  NodeSpline::Ptr base_linear_;     ///< the linear position of the base.
  EulerConverter base_angular_; 	///< the orientation of the base.
  NodeSpline::Ptr ee_motion_;       ///< the linear position of the endeffectors.
  Parameters params_;

  int n_constraints_per_node_;

  EE ee_;
  std::vector<double> T_;

  // see TimeDiscretizationConstraint for documentation
  void UpdateConstraintAtInstance (double t, int k, VectorXd& g) const override;
  void UpdateBoundsAtInstance (double t, int k, VecBound&) const override;
  void UpdateJacobianAtInstance(double t, int k, std::string, Jacobian&) const override;

  VecBound node_bounds_;     ///< same bounds for each discretized node
  int GetRow(int node, int dimension) const;
};

} /* namespace towr */

#endif /* TOWR_TOWR_INCLUDE_TOWR_CONSTRAINTS_DRIVE_CONSTRAINT_H_ */
