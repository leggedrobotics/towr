/*
 * terrain_constraint_discretized.cc
 *
 *  Created on: mai, 2020
 *  Author: cla
 */
#include <towr/constraints/terrain_constraint_discretized.h>

namespace towr {

TerrainConstraintDiscretized::TerrainConstraintDiscretized(
    const HeightMap::Ptr &terrain, double T, double dt, const EE &ee,
    const SplineHolder &spline_holder)
    : TimeDiscretizationConstraint(
          T, dt, "terrain-discretized-" + std::to_string(ee)) {
  ee_ = ee;
  ee_wheels_motion_ = spline_holder.ee_motion_.at(ee_);
  terrain_ = terrain;
  decision_ = spline_holder.ee_decision_.at(ee_);
  n_constraints_per_node_ = 2; // lateral velocity and acceleration
  SetRows(GetNumberOfNodes() * n_constraints_per_node_);
}

void TerrainConstraintDiscretized::UpdateConstraintAtInstance(
    double t, int k, VectorXd &g) const {
  int row = k * n_constraints_per_node_;
  Vector3d ee_d = decision_->GetPoint(t).p();
  Vector3d p = ee_wheels_motion_->GetPoint(t).p();
  double terrain_h =
      terrain_->GetHeight(p.x(),
                          p.y());
  g(row++) = ee_d.x() * (p.z() - terrain_h);
  g(row++) =
      (1 - ee_d.x()) * (p.z() - terrain_h);
}

void TerrainConstraintDiscretized::UpdateBoundsAtInstance(
    double t, int k, VecBound &bounds) const {
  int row = k * n_constraints_per_node_;
  bounds.at(row++) = ifopt::Bounds(-0.01, 0.01);
  bounds.at(row++) = ifopt::BoundGreaterZero;
}

void TerrainConstraintDiscretized::UpdateJacobianAtInstance(
    double t, int k, std::string var_set, Jacobian &jac) const {

  Vector3d p = ee_wheels_motion_->GetPoint(t).p();
  Vector3d ee_d = decision_->GetPoint(t).p();

  double terrain_h =
      terrain_->GetHeight(p.x(),
                          p.y());

  double dterrain_h_dx = terrain_->GetDerivativeOfHeightWrt(X_, p.x(), p.y());
  double dterrain_h_dy = terrain_->GetDerivativeOfHeightWrt(Y_, p.x(), p.y());

  int row = k * n_constraints_per_node_;

  if (var_set == id::EEMotionNodes(ee_)) {
    jac.row(row++) =
        ee_d.x() *
        ((ee_wheels_motion_->GetJacobianWrtNodes(t, kPos)).row(Z) -
            dterrain_h_dx *
             (ee_wheels_motion_->GetJacobianWrtNodes(t, kPos)).row(X) -
            dterrain_h_dy *
             (ee_wheels_motion_->GetJacobianWrtNodes(t, kPos)).row(Y));
    jac.row(row++) =
        ((1 - ee_d.x()) *
         ((ee_wheels_motion_->GetJacobianWrtNodes(t, kPos)).row(Z) -
             dterrain_h_dx *
              (ee_wheels_motion_->GetJacobianWrtNodes(t, kPos)).row(X) -
             dterrain_h_dy *
              (ee_wheels_motion_->GetJacobianWrtNodes(t, kPos)).row(Y)));
  }

  if (var_set == id::EEDecision(ee_)) {
    //this values are not optimized over generally, but the jacobian has been put here still if someone were to change that
    jac.row(row++) = decision_->GetJacobianWrtNodes(t, kPos).row(X) *
                     (p.z() - terrain_h);

    jac.row(row++) = -decision_->GetJacobianWrtNodes(t, kPos).row(X) *
                     (p.z() - terrain_h);
  }


  if (var_set == id::EESchedule(ee_)) {
    jac.row(row++) =(
        ee_d.x() *
        ((ee_wheels_motion_->GetJacobianOfPosWrtDurations(t)).row(Z) -
            dterrain_h_dx *
         (ee_wheels_motion_->GetJacobianOfPosWrtDurations(t)).row(X) -
            dterrain_h_dy *
         (ee_wheels_motion_->GetJacobianOfPosWrtDurations(t)).row(Y)))
             +(decision_->GetJacobianWrtNodes(t, kPos).row(X) *
               (p.z() - terrain_h));
    jac.row(row++) =(        ((1 - ee_d.x()) *
         ((ee_wheels_motion_->GetJacobianOfPosWrtDurations(t)).row(Z) -
             dterrain_h_dx *
          (ee_wheels_motion_->GetJacobianOfPosWrtDurations(t)).row(X) -
             dterrain_h_dy *
          (ee_wheels_motion_->GetJacobianOfPosWrtDurations(t)).row(Y))))
          +(-decision_->GetJacobianWrtNodes(t, kPos).row(X) *
            (p.z() - terrain_h));
  }




}

} /* namespace towr */
