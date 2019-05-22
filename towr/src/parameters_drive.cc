/*
 * parameters_drive.cc
 *
 *  Created on: Mar 29, 2019
 *      Author: vivian
 */

#include "towr/parameters_drive.h"
#include <towr/variables/cartesian_dimensions.h>

namespace towr {

ParametersDrive::ParametersDrive()
{
  duration_ee_polynomial_ = duration_base_polynomial_;
  wheels_torque_limit_ = 32.0;
  max_wheels_acc_ = {6.0, 0.0, 6.0};
  max_base_acc_lin_ = {8.0, 8.0, 8.0};
  max_base_acc_ang_ = {5.0, 5.0, 5.0};
  wheels_radius_ = 0.08;
  force_limit_in_x_direction_ = 400.0;

//  bounds_final_lin_pos_ = {X,Y,Z};   // adds the bound on the final z-position of the base

  DeleteAllConstraints();  // clear the constraints initialized in the base class

  // constraints for driving motions
  constraints_.push_back(TerrainWheels);   		// ensures terrain contact during the entire motion
  constraints_.push_back(DynamicWheels);   		// ensures that the dynamic model is fullfilled at discrete times.
  constraints_.push_back(BaseAcc); 		   		// so accelerations don't jump between polynomials
  constraints_.push_back(EndeffectorAcc);  		// so accelerations don't jump between polynomials
  constraints_.push_back(ForceWheels); 	   		// ensures unilateral forces and inside the friction cone
  constraints_.push_back(WheelsAccLimits); 		// constrain the acceleration on the wheels
  constraints_.push_back(BaseAccLimits);   		// constrain the acceleration of the base
  constraints_.push_back(Stability);       		// ensure stability margin (weird roll angles)
//  constraints_.push_back(WheelsNonHolonomic);  	// non-holonomic driving constraint

  // only one of these constraints are defined according to the config file
//  constraints_.push_back(WheelsMotion);    // constrain minimum and maximum legs extension
//  constraints_.push_back(EndeffectorRom);  // ensures that the range of motion of the wheels is respected at discrete times.

}

void
ParametersDrive::SetWheelsMotionConstraint () {
  constraints_.push_back(WheelsMotion);
}

void
ParametersDrive::SetEndeffectorRomConstraint () {
  constraints_.push_back(EndeffectorRom);
}

ParametersDrive::VecTimes
ParametersDrive::GetEEWheelsPolyDurations() const
{
  std::vector<double> ee_spline_timings_;
  double dt = duration_ee_polynomial_;
  double t_left = total_time_;

  double eps = 1e-10; // since repeated subtraction causes inaccuracies
  while (t_left > eps) {
    double duration = t_left>dt?  dt : t_left;
    ee_spline_timings_.push_back(duration);

    t_left -= dt;
  }

  return ee_spline_timings_;
}

ParametersDrive::VecTimes
ParametersDrive::GetBasePolyDurations () const
{
  std::vector<double> base_spline_timings_;
  double dt = duration_base_polynomial_;
  double t_left = total_time_;

  double eps = 1e-10; // since repeated subtraction causes inaccuracies
  while (t_left > eps) {
    double duration = t_left>dt?  dt : t_left;
    base_spline_timings_.push_back(duration);

    t_left -= dt;
  }

  return base_spline_timings_;
}

double
ParametersDrive::GetTotalTime () const
{
  return total_time_;
}

int
ParametersDrive::GetEECount() const
{
  return n_ee_;
}

void
ParametersDrive::DeleteAllConstraints()
{
  constraints_.clear();
}

} // namespace towr
