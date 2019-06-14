/*
 * parameters_driving.h
 *
 *  Created on: Mar 29, 2019
 *      Author: vivian
 */

#ifndef TOWR_INCLUDE_TOWR_PARAMETERS_DRIVE_H_
#define TOWR_INCLUDE_TOWR_PARAMETERS_DRIVE_H_

#include "towr/parameters.h"
#include <iostream>
#include <memory>

namespace towr {

class ParametersDrive : public Parameters
{
public:
  using Ptr  = std::shared_ptr<ParametersDrive>;
  using Base = Parameters;

  ParametersDrive();
  virtual ~ParametersDrive() = default;

  /// Number of end-effectors
  int n_ee_;

  /// Total time of the driving motion
  double total_time_;

  /// True if the motion is pure driving
  bool is_pure_driving_motion_ = true;

  /// Fixed duration of each cubic polynomial describing the ee motion in driving condition.
  double duration_ee_polynomial_;

  /// The maximum allowable force [N] in tangent direction
  double force_limit_in_x_direction_;

  /// The maximum allowable torque [Nm] on each wheel
  double wheels_torque_limit_;

  /// The wheels radius [m]
  double wheels_radius_;

  /// Non-holonomic constraint
  double use_non_holonomic_constraint_;

  /// Maximum acceleration of the wheels [m/s^2] (x, y, z)
  std::vector<double> max_wheels_acc_;

  /// Maximum acceleration of the base (x, y, z)
  std::vector<double> max_base_acc_lin_;

  /// Maximum acceleration of the base (roll, pitch, yaw)
  std::vector<double> max_base_acc_ang_;

  /// The durations of each ee polynomial in the spline (motion and force).
  VecTimes GetEEWheelsPolyDurations() const;

  /// The durations of each base polynomial in the spline (lin+ang).
  VecTimes GetBasePolyDurations() const override;

  /// Total duration [s] of the motion.
  double GetTotalTime() const override;

  /// The number of end-effectors.
  int GetEECount() const override;

  /// Clear the constraints initialized by the Base class
  void DeleteAllConstraints();

  void SetWheelsMotionConstraint ();
  void SetEndeffectorRomConstraint ();
  void SetNonHolonomicConstraint ();

};

} // namespace towr


#endif /* TOWR_INCLUDE_TOWR_PARAMETERS_DRIVE_H_ */
