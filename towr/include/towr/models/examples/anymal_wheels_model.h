/*
 * anymal_wheels_model.h
 *
 *  Created on: Apr 4, 2019
 *      Author: vivian
 */

#ifndef TOWR_INCLUDE_TOWR_MODELS_EXAMPLES_ANYMAL_WHEELS_MODEL_H_
#define TOWR_INCLUDE_TOWR_MODELS_EXAMPLES_ANYMAL_WHEELS_MODEL_H_

#include <towr/models/kinematic_model.h>
#include <towr/models/single_rigid_body_dynamics.h>
#include <towr/models/endeffector_mappings.h>

namespace towr {

/**
 * @brief The Kinematics of the quadruped robot ANYmal with wheels.
 */
class AnymalWheelsKinematicModel : public KinematicModel {
public:
  AnymalWheelsKinematicModel () : KinematicModel(4)
  {

    const double x_nominal_b = 0.353088;
    const double y_nominal_b = 0.146229;
    const double z_nominal_b = -0.5; //-0.60174258;
    //const double z_nominal_b = -0.55; //viv's parameter

    double offsets = 0.095;// 0.1;

    nominal_stance_.at(LF) <<  x_nominal_b + offsets,   y_nominal_b, z_nominal_b;
    nominal_stance_.at(RF) <<  x_nominal_b - offsets,  -y_nominal_b, z_nominal_b;
    nominal_stance_.at(LH) << -x_nominal_b + offsets,   y_nominal_b, z_nominal_b;
    nominal_stance_.at(RH) << -x_nominal_b - offsets,  -y_nominal_b, z_nominal_b;

    //max_dev_from_nominal_ << 0.225-offsets, 0.095, 0.095; // viv's parameter
    max_dev_from_nominal_ << 0.25-offsets, 0.095, 0.15;

    const double x_nominal_hip = 0.3405;
    const double y_nominal_hip = y_nominal_b; //0.1710;
    const double z_nominal_hip = 0.0;
  }

};

/**
 * @brief The Dynamics of the quadruped robot ANYmal with wheels.
 */
class AnymalWheelsDynamicModel : public SingleRigidBodyDynamics {
public:
  AnymalWheelsDynamicModel()
      : SingleRigidBodyDynamics(46.6237, 1.72911, 4.48562, 4.38949,
                                0.0022872, 0.0752196, -0.0141115, 4) {}
};

} // namespace towr



#endif /* TOWR_INCLUDE_TOWR_MODELS_EXAMPLES_ANYMAL_WHEELS_MODEL_H_ */
