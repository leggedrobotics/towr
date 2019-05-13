/*
 * terrain_wheels_constraint.cc
 *
 *  Created on: Mar 29, 2019
 *      Author: vivian
 */

#include <iostream>
#include <towr/constraints/terrain_wheels_constraint.h>
#include <towr/variables/variable_names.h>

namespace towr {

TerrainWheelsConstraint::TerrainWheelsConstraint (const HeightMap::Ptr& terrain,
        										  double T, double dt, const EE& ee,
												  const SplineHolderDrive& spline_holder)
	:TimeDiscretizationConstraint(T, dt, "terrain-wheels-" + std::to_string(ee))
{
  ee_wheels_motion_ = spline_holder.ee_wheels_motion_.at(ee);
  terrain_ = terrain;
  ee_ = ee;

  SetRows(GetNumberOfNodes());   // one constraint per node
}

void
TerrainWheelsConstraint::UpdateConstraintAtInstance (double t, int k, VectorXd& g) const
{
  Vector3d pos_ee_W = ee_wheels_motion_->GetPoint(t).p();

  g(k) = pos_ee_W.z() - terrain_->GetHeight(pos_ee_W.x(), pos_ee_W.y());
}

void
TerrainWheelsConstraint::UpdateBoundsAtInstance (double t, int k, VecBound& bounds) const
{
  bounds.at(k) = ifopt::BoundZero;
}

TerrainWheelsConstraint::Jacobian
TerrainWheelsConstraint::GetJacobianTerrainHeight(double x, double y) const
{
  Jacobian jac = Eigen::Vector3d::Ones().transpose().sparseView();

  jac.coeffRef(0, X_) = -terrain_->GetDerivativeOfHeightWrt(X_, x, y);
  jac.coeffRef(0, Y_) = -terrain_->GetDerivativeOfHeightWrt(Y_, x, y);

  return jac;
}

void
TerrainWheelsConstraint::UpdateJacobianAtInstance(double t, int k, std::string var_set, Jacobian& jac) const
{
  if (var_set == id::EEWheelsMotionNodes(ee_)) {
    Vector3d pos = ee_wheels_motion_->GetPoint(t).p();
    jac.row(k) = GetJacobianTerrainHeight(pos.x(), pos.y()) * ee_wheels_motion_->GetJacobianWrtNodes(t, kPos);
  }
}

} /* namespace towr */



