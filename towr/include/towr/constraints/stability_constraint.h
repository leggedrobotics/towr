/*
 * stability_constraint.h
 *
 *  Created on: Apr 12, 2019
 *      Author: vivian
 */

#ifndef TOWR_INCLUDE_TOWR_CONSTRAINTS_STABILITY_CONSTRAINT_H_
#define TOWR_INCLUDE_TOWR_CONSTRAINTS_STABILITY_CONSTRAINT_H_

#include <towr/variables/spline.h>
#include <towr/variables/spline_holder.h>
#include <towr/variables/spline_holder_drive.h>
#include <towr/variables/euler_converter.h>

#include <towr/models/dynamic_model.h>

#include "time_discretization_constraint.h"

namespace towr {

class StabilityConstraint : public TimeDiscretizationConstraint {
public:
  using MatrixXd = Eigen::MatrixXd;
  using Vector3d = Eigen::Vector3d;

  StabilityConstraint (const DynamicModel::Ptr& model,
                       double T, double dt,
                       const SplineHolderDrive& spline_holder);
  virtual ~StabilityConstraint () = default;

private:
  NodeSpline::Ptr base_linear_;   ///< lin. base pos/vel/acc in world frame
  EulerConverter base_angular_;        ///< angular base state
  std::vector<NodeSpline::Ptr> ee_forces_; ///< endeffector forces in world frame.
  std::vector<NodeSpline::Ptr> ee_motion_; ///< endeffector position in world frame.

  mutable DynamicModel::Ptr model_;    ///< the dynamic model (e.g. Centroidal)

  void UpdateModel(double t) const;

  void UpdateConstraintAtInstance(double t, int k, VectorXd& g) const override;
  void UpdateBoundsAtInstance(double t, int k, VecBound& bounds) const override;
  void UpdateJacobianAtInstance(double t, int k, std::string, Jacobian&) const override;
};

} /* namespace towr */

#endif /* TOWR_INCLUDE_TOWR_CONSTRAINTS_STABILITY_CONSTRAINT_H_ */
