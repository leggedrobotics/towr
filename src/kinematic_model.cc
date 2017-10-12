/**
 @file    kinematic_model.cc
 @author  Alexander W. Winkler (winklera@ethz.ch)
 @date    Sep 19, 2017
 @brief   Brief description
 */

#include <xpp/models/kinematic_model.h>

#include <Eigen/Dense>

namespace xpp {
namespace opt {

KinematicModel::KinematicModel (int n_ee)
{
  nominal_stance_.SetCount(n_ee);
  max_dev_from_nominal_.setZero();
}

KinematicModel::~KinematicModel ()
{
  // TODO Auto-generated destructor stub
}

} /* namespace opt */
} /* namespace xpp */