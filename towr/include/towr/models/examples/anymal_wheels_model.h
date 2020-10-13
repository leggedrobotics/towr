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
    const double z_nominal_b = -0.60174258;

    double offsets = 0.0;// 0.1;

    nominal_stance_.at(LF) <<  x_nominal_b + offsets,   y_nominal_b, z_nominal_b;
    nominal_stance_.at(RF) <<  x_nominal_b - offsets,  -y_nominal_b, z_nominal_b;
    nominal_stance_.at(LH) << -x_nominal_b + offsets,   y_nominal_b, z_nominal_b;
    nominal_stance_.at(RH) << -x_nominal_b - offsets,  -y_nominal_b, z_nominal_b;

    max_dev_from_nominal_ << 0.225-offsets, 0.095, 0.095;
//    max_dev_from_nominal_ << 0.3, 0.25, 0.25;  // if want cross over stairs, try increasing constraint box

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
      : SingleRigidBodyDynamics(25.586298, 0.17825624, 1.7073739, 1.7424951,
                                -0.0038845614, -0.042678002, 0.002850078, 4) {}
};

} // namespace towr



#endif /* TOWR_INCLUDE_TOWR_MODELS_EXAMPLES_ANYMAL_WHEELS_MODEL_H_ */
