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
    const double x_nominal_b = 0.35;
    const double y_nominal_b = 0.235; // 0.21;
    const double z_nominal_b = -0.44;

    nominal_stance_.at(LF) <<  x_nominal_b,   y_nominal_b, z_nominal_b;
    nominal_stance_.at(RF) <<  x_nominal_b,  -y_nominal_b, z_nominal_b;
    nominal_stance_.at(LH) << -x_nominal_b,   y_nominal_b, z_nominal_b;
    nominal_stance_.at(RH) << -x_nominal_b,  -y_nominal_b, z_nominal_b;

    max_dev_from_nominal_ << 0.15, 0.1, 0.10;

    const double x_nominal_hip = 0.3405;
    const double y_nominal_hip = y_nominal_b; //0.1710;
    const double z_nominal_hip = 0.0;

   //hip_pos_wrt_base_.at(LF) <<  x_nominal_hip,  y_nominal_hip, z_nominal_hip;
   //hip_pos_wrt_base_.at(RF) <<  x_nominal_hip, -y_nominal_hip, z_nominal_hip;
   //hip_pos_wrt_base_.at(LH) << -x_nominal_hip,  y_nominal_hip, z_nominal_hip;
   //hip_pos_wrt_base_.at(RH) << -x_nominal_hip, -y_nominal_hip, z_nominal_hip;
  }

};

/**
 * @brief The Dynamics of the quadruped robot ANYmal with wheels.
 */
class AnymalWheelsDynamicModel : public SingleRigidBodyDynamics {
public:
  AnymalWheelsDynamicModel()
  : SingleRigidBodyDynamics(19.642, 0.248057547486776, 0.650151658461464, 0.620944203624185,
		  	  	  	  	    0.001097838800893, -0.003945011648535, -0.002135691054868, 4) {}
};

} // namespace towr



#endif /* TOWR_INCLUDE_TOWR_MODELS_EXAMPLES_ANYMAL_WHEELS_MODEL_H_ */
