/*
 * terrain_constraint_discretized.cc
 *
 *  Created on: mai, 2020
 *  Author: cla
 */
#include <towr/constraints/force_constraint_discretized.h>

namespace towr {

ForceConstraintDiscretized::ForceConstraintDiscretized(
    const HeightMap::Ptr &terrain, double T, double dt, const EE &ee,
    const SplineHolder &spline_holder, double force_limit)
    : TimeDiscretizationConstraint(
          T, dt, "force-discretized-" + std::to_string(ee)) {
  ee_ = ee;
  ee_wheels_motion_ = spline_holder.ee_motion_.at(ee_);
  ee_force_ = spline_holder.ee_force_.at(ee_);
  terrain_ = terrain;

  decision_ = spline_holder.ee_decision_.at(ee_);

  n_constraints_per_node_ = 5; // lateral velocity and acceleration
  mu_ = terrain->GetFrictionCoeff();
  fn_max_ = force_limit;
  SetRows(GetNumberOfNodes() * n_constraints_per_node_);
}

void ForceConstraintDiscretized::UpdateConstraintAtInstance(double t, int k,
                                                            VectorXd &g) const {
  Vector3d ee_d = decision_->GetPoint(t).p();
  int row = k * n_constraints_per_node_;
  Vector3d p = ee_wheels_motion_->GetPoint(t).p(); // doesn't change during stance phase
  Vector3d n = terrain_->GetNormalizedBasis(HeightMap::Normal, p.x(), p.y());
  Vector3d f = ee_force_->GetPoint(t).p();


  // unilateral force
  g(row++) = ee_d.x() * (f.x() * n.x() + f.y() * n.y() +
                         f.z() * n.z()); // >0 (unilateral forces)

  // frictional pyramid
  Vector3d t1 = terrain_->GetNormalizedBasis(HeightMap::Tangent1, p.x(), p.y());
  g(row++) = ee_d.x() *
             (f.x() * (t1.x() - mu_ * n.x()) + f.y() * (t1.y() - mu_ * n.y()) +
              f.z() * (t1.z() - mu_ * n.z())); // t1 < mu*n
  g(row++) = ee_d.x() *
             (f.x() * (t1.x() + mu_ * n.x()) + f.y() * (t1.y() + mu_ * n.y()) +
              f.z() * (t1.z() + mu_ * n.z())); // t1 > -mu*n

  Vector3d t2 = terrain_->GetNormalizedBasis(HeightMap::Tangent2, p.x(), p.y());
  g(row++) = ee_d.x() *
             (f.x() * (t2.x() - mu_ * n.x()) + f.y() * (t2.y() - mu_ * n.y()) +
              f.z() * (t2.z() - mu_ * n.z())); // t2 < mu*n
  g(row++) = ee_d.x() *
             (f.x() * (t2.x() + mu_ * n.x()) + f.y() * (t2.y() + mu_ * n.y()) +
              f.z() * (t2.z() + mu_ * n.z())); // t2 > -mu*n

} // namespace towr

void ForceConstraintDiscretized::UpdateBoundsAtInstance(
    double t, int k, VecBound &bounds) const {
  int row = k * n_constraints_per_node_;
  bounds.at(row++) = ifopt::Bounds(0.0, fn_max_); // unilateral forces
  bounds.at(row++) = ifopt::BoundSmallerZero;     // f_t1 <  mu*n
  bounds.at(row++)  = ifopt::BoundGreaterZero;     // f_t1 > -mu*n
  bounds.at(row++)  = ifopt::BoundSmallerZero;     // f_t2 <  mu*n
  bounds.at(row++)  = ifopt::BoundGreaterZero;     // f_t2 > -mu*n
}

void ForceConstraintDiscretized::UpdateJacobianAtInstance(double t, int k,
                                                          std::string var_set,
                                                          Jacobian &jac) const {


  int n_jac = jac.cols();
  int row = k * n_constraints_per_node_;
  Vector3d ee_d = decision_->GetPoint(t).p();
  Vector3d p = ee_wheels_motion_->GetPoint(t).p(); // doesn't change during stance phase
  Vector3d n = terrain_->GetNormalizedBasis(HeightMap::Normal, p.x(), p.y());
  Vector3d f = ee_force_->GetPoint(t).p();

  Vector3d t1 = terrain_->GetNormalizedBasis(HeightMap::Tangent1, p.x(), p.y());
  Vector3d t2 = terrain_->GetNormalizedBasis(HeightMap::Tangent2, p.x(), p.y());




  Vector3d dndx = terrain_->GetDerivativeOfNormalizedBasisWrt(
      HeightMap::Normal, X_, p.x(), p.y());
  Vector3d dndy = terrain_->GetDerivativeOfNormalizedBasisWrt(
      HeightMap::Normal, Y_, p.x(), p.y());


  Vector3d dt1dx = terrain_->GetDerivativeOfNormalizedBasisWrt(
      HeightMap::Tangent1, X_, p.x(), p.y());
  Vector3d dt1dy = terrain_->GetDerivativeOfNormalizedBasisWrt(
      HeightMap::Tangent1, Y_, p.x(), p.y());

  Vector3d dt2dx = terrain_->GetDerivativeOfNormalizedBasisWrt(
      HeightMap::Tangent2, X_, p.x(), p.y());
  Vector3d dt2dy = terrain_->GetDerivativeOfNormalizedBasisWrt(
      HeightMap::Tangent2, Y_, p.x(), p.y());


  if (var_set == id::EEForceNodes(ee_)) {

    Jacobian dforce_dforcenode_x(1, n_jac);
    Jacobian dforce_dforcenode_y(1, n_jac);
    Jacobian dforce_dforcenode_z(1, n_jac);
    dforce_dforcenode_x = ee_force_->GetJacobianWrtNodes(t, kPos).row(X);
    dforce_dforcenode_y = ee_force_->GetJacobianWrtNodes(t, kPos).row(Y);
    dforce_dforcenode_z = ee_force_->GetJacobianWrtNodes(t, kPos).row(Z);

    // unilateral force
    jac.row(row++) =
        ee_d.x() * (dforce_dforcenode_x * n.x() + dforce_dforcenode_y * n.y() +
                    dforce_dforcenode_z * n.z()); // >0 (unilateral forces)

    // frictional pyramid
    jac.row(row++) =
        ee_d.x() * (dforce_dforcenode_x * (t1.x() - mu_ * n.x()) +
                    dforce_dforcenode_y * (t1.y() - mu_ * n.y()) +
                    dforce_dforcenode_z * (t1.z() - mu_ * n.z())); // t1 < mu*n
    jac.row(row++) =
        ee_d.x() * (dforce_dforcenode_x * (t1.x() + mu_ * n.x()) +
                    dforce_dforcenode_y * (t1.y() + mu_ * n.y()) +
                    dforce_dforcenode_z * (t1.z() + mu_ * n.z())); // t1 > -mu*n

    jac.row(row++) =
        ee_d.x() * (dforce_dforcenode_x * (t2.x() - mu_ * n.x()) +
                    dforce_dforcenode_y * (t2.y() - mu_ * n.y()) +
                    dforce_dforcenode_z * (t2.z() - mu_ * n.z())); // t2 < mu*n
    jac.row(row++) =
        ee_d.x() * (dforce_dforcenode_x * (t2.x() + mu_ * n.x()) +
                    dforce_dforcenode_y * (t2.y() + mu_ * n.y()) +
                    dforce_dforcenode_z * (t2.z() + mu_ * n.z())); // t2 > -mu*n
  }


  if (var_set == id::EEMotionNodes(ee_)) {


    Jacobian dmotion_dmotionnode_x(1, n_jac);
    Jacobian dmotion_dmotionnode_y(1, n_jac);
    dmotion_dmotionnode_x = ee_wheels_motion_->GetJacobianWrtNodes(t, kPos).row(X);
    dmotion_dmotionnode_y = ee_wheels_motion_->GetJacobianWrtNodes(t, kPos).row(Y);

    Jacobian dn_dmotionnode_x(1, n_jac);
    Jacobian dn_dmotionnode_y(1, n_jac);
    Jacobian dn_dmotionnode_z(1, n_jac);
    dn_dmotionnode_x =
        dndx.x() * dmotion_dmotionnode_x + dndy.x() * dmotion_dmotionnode_y;
    dn_dmotionnode_y =
        dndx.y() * dmotion_dmotionnode_x + dndy.y() * dmotion_dmotionnode_y;
    dn_dmotionnode_z =
        dndx.z() * dmotion_dmotionnode_x + dndy.z() * dmotion_dmotionnode_y;




    Jacobian dt1_dmotionnode_x(1, n_jac);
    Jacobian dt1_dmotionnode_y(1, n_jac);
    Jacobian dt1_dmotionnode_z(1, n_jac);
    dt1_dmotionnode_x =
        dt1dx.x() * dmotion_dmotionnode_x + dt1dy.x() * dmotion_dmotionnode_y;
    dt1_dmotionnode_y =
        dt1dx.y() * dmotion_dmotionnode_x + dt1dy.y() * dmotion_dmotionnode_y;
    dt1_dmotionnode_z =
        dt1dx.z() * dmotion_dmotionnode_x + dt1dy.z() * dmotion_dmotionnode_y;

    Jacobian dt2_dmotionnode_x(1, n_jac);
    Jacobian dt2_dmotionnode_y(1, n_jac);
    Jacobian dt2_dmotionnode_z(1, n_jac);
    dt2_dmotionnode_x =
       dt2dx.x() * dmotion_dmotionnode_x + dt2dy.x() * dmotion_dmotionnode_y;
    dt2_dmotionnode_y =
        dt2dx.y() * dmotion_dmotionnode_x + dt2dy.y() * dmotion_dmotionnode_y;
    dt2_dmotionnode_z =
        dt2dx.z() * dmotion_dmotionnode_x + dt2dy.z() * dmotion_dmotionnode_y;

    // unilateral force
    jac.row(row++) =
        ee_d.x() * (dn_dmotionnode_x * f.x() + dn_dmotionnode_y * f.y() +
                    dn_dmotionnode_z * f.z()); // >0 (unilateral forces)

    // frictional pyramid
    jac.row(row++) =
        ee_d.x() * (f.x() * (dt1_dmotionnode_x - mu_ * dn_dmotionnode_x) +
                    f.y() * (dt1_dmotionnode_y - mu_ * dn_dmotionnode_y) +
                    f.z() * (dt1_dmotionnode_z - mu_ * dn_dmotionnode_z));
    jac.row(row++) =
        ee_d.x() * (f.x() * (dt1_dmotionnode_x + mu_ * dn_dmotionnode_x) +
                    f.y() * (dt1_dmotionnode_y + mu_ * dn_dmotionnode_y) +
                    f.z() * (dt1_dmotionnode_z + mu_ * dn_dmotionnode_z));
    jac.row(row++) =
        ee_d.x() * (f.x() * (dt2_dmotionnode_x - mu_ * dn_dmotionnode_x) +
                    f.y() * (dt2_dmotionnode_y - mu_ * dn_dmotionnode_y) +
                    f.z() * (dt2_dmotionnode_z - mu_ * dn_dmotionnode_z));
    jac.row(row++) =
        ee_d.x() * (f.x() * (dt2_dmotionnode_x + mu_ * dn_dmotionnode_x) +
                    f.y() * (dt2_dmotionnode_y + mu_ * dn_dmotionnode_y) +
                    f.z() * (dt2_dmotionnode_z + mu_ * dn_dmotionnode_z));

    }
  if (var_set == id::EEDecision(ee_)) {
    //this values are not optimized over generally, but the jacobian has been put here still

    // unilateral force
    jac.row(row++) = decision_->GetJacobianWrtNodes(t, kPos).row(X) *
                     (f.x() * n.x() + f.y() * n.y() +
                      f.z() * n.z()); // >0 (unilateral forces)

    // frictional pyramid
    jac.row(row++) =
        decision_->GetJacobianWrtNodes(t, kPos).row(X) *
        static_cast<double>(f.transpose() * (t1 - mu_ * n)); // t1 < mu*n
    jac.row(row++) =
        decision_->GetJacobianWrtNodes(t, kPos).row(X) *
        static_cast<double>(f.transpose() * (t1 + mu_ * n)); // t1 > -mu*n

    jac.row(row++) =
        decision_->GetJacobianWrtNodes(t, kPos).row(X) *
        static_cast<double>(f.transpose() * (t2 - mu_ * n)); // t2 < mu*n
    jac.row(row++) =
        decision_->GetJacobianWrtNodes(t, kPos).row(X) *
        static_cast<double>(f.transpose() * (t2 + mu_ * n)); // t2 > -mu*n
  }
  if (var_set == id::EESchedule(ee_)) {






    Jacobian dforce_dduration_x(1, n_jac);
    Jacobian dforce_dduration_y(1, n_jac);
    Jacobian dforce_dduration_z(1, n_jac);
    dforce_dduration_x = ee_force_->GetJacobianOfPosWrtDurations(t).row(X);
    dforce_dduration_y = ee_force_->GetJacobianOfPosWrtDurations(t).row(Y);
    dforce_dduration_z = ee_force_->GetJacobianOfPosWrtDurations(t).row(Z);

    Jacobian dmotion_dduration_x(1, n_jac);
    Jacobian dmotion_dduration_y(1, n_jac);
    dmotion_dduration_x = ee_wheels_motion_->GetJacobianOfPosWrtDurations(t).row(X);
    dmotion_dduration_y = ee_wheels_motion_->GetJacobianOfPosWrtDurations(t).row(Y);

    Jacobian dn_dduration_x(1, n_jac);
    Jacobian dn_dduration_y(1, n_jac);
    Jacobian dn_dduration_z(1, n_jac);
    dn_dduration_x =
        dndx.x() * dmotion_dduration_x + dndy.x() * dmotion_dduration_y;
    dn_dduration_y =
        dndx.y() * dmotion_dduration_x + dndy.y() * dmotion_dduration_y;
    dn_dduration_z =
        dndx.z() * dmotion_dduration_x + dndy.z() * dmotion_dduration_y;


    Jacobian dt1_dduration_x(1, n_jac);
    Jacobian dt1_dduration_y(1, n_jac);
    Jacobian dt1_dduration_z(1, n_jac);
    dt1_dduration_x =
        dt1dx.x() * dmotion_dduration_x + dt1dy.x() * dmotion_dduration_y;
    dt1_dduration_y =
        dt1dx.y() * dmotion_dduration_x + dt1dy.y() * dmotion_dduration_y;
    dt1_dduration_z =
        dt1dx.z() * dmotion_dduration_x + dt1dy.z() * dmotion_dduration_y;

    Jacobian dt2_dduration_x(1, n_jac);
    Jacobian dt2_dduration_y(1, n_jac);
    Jacobian dt2_dduration_z(1, n_jac);
    dt2_dduration_x =
        dt2dx.x() * dmotion_dduration_x + dt2dy.x() * dmotion_dduration_y;
    dt2_dduration_y =
        dt2dx.y() * dmotion_dduration_x + dt2dy.y() * dmotion_dduration_y;
    dt2_dduration_z =
        dt2dx.z() * dmotion_dduration_x + dt2dy.z() * dmotion_dduration_y;

    // unilateral force
    jac.row(row++) =
        (ee_d.x() * (dforce_dduration_x * n.x() + dforce_dduration_y * n.y() +
                     dforce_dduration_z * n.z())) +
        (ee_d.x() * (dn_dduration_x * f.x() + dn_dduration_y * f.y() +
            dn_dduration_z * f.z())) +
        (decision_->GetJacobianWrtNodes(t, kPos).row(X) *
         (f.x() * n.x() + f.y() * n.y() +
          f.z() * n.z())); // >0 (unilateral forces)

    // frictional pyramid
    jac.row(row++) =
        (ee_d.x() * (dforce_dduration_x * (t1.x() - mu_ * n.x()) +
                     dforce_dduration_y * (t1.y() - mu_ * n.y()) +
                     dforce_dduration_z * (t1.z() - mu_ * n.z()))) +
        (ee_d.x() * (f.x() * (dt1_dduration_x - mu_ * dn_dduration_x) +
                     f.y() * (dt1_dduration_y - mu_ * dn_dduration_y) +
                     f.z() * (dt1_dduration_z - mu_ * dn_dduration_z))) +
        (decision_->GetJacobianWrtNodes(t, kPos).row(X) *
         static_cast<double>(f.transpose() * (t1 - mu_ * n))); // t1 < mu*n
    jac.row(row++) =
        (ee_d.x() * (dforce_dduration_x * (t1.x() + mu_ * n.x()) +
                     dforce_dduration_y * (t1.y() + mu_ * n.y()) +
                     dforce_dduration_z * (t1.z() + mu_ * n.z()))) +
        (ee_d.x() * (f.x() * (dt1_dduration_x + mu_ * dn_dduration_x) +
                     f.y() * (dt1_dduration_y + mu_ * dn_dduration_y) +
                     f.z() * (dt1_dduration_z + mu_ * dn_dduration_z))) +
        (decision_->GetJacobianWrtNodes(t, kPos).row(X) *
         static_cast<double>(f.transpose() * (t1 + mu_ * n))); // t1 > -mu*n

    jac.row(row++) =
        (ee_d.x() * (dforce_dduration_x * (t2.x() - mu_ * n.x()) +
                     dforce_dduration_y * (t2.y() - mu_ * n.y()) +
                     dforce_dduration_z * (t2.z() - mu_ * n.z()))) +
        (ee_d.x() * (f.x() * (dt2_dduration_x - mu_ * dn_dduration_x) +
                     f.y() * (dt2_dduration_y - mu_ * dn_dduration_y) +
                     f.z() * (dt2_dduration_z - mu_ * dn_dduration_z))) +
        (decision_->GetJacobianWrtNodes(t, kPos).row(X) *
         static_cast<double>(f.transpose() * (t2 - mu_ * n))); // t2 < mu*n
    jac.row(row++) =
        (ee_d.x() * (dforce_dduration_x * (t2.x() + mu_ * n.x()) +
                     dforce_dduration_y * (t2.y() + mu_ * n.y()) +
                     dforce_dduration_z * (t2.z() + mu_ * n.z()))) +
        (ee_d.x() * (f.x() * (dt2_dduration_x + mu_ * dn_dduration_x) +
                     f.y() * (dt2_dduration_y + mu_ * dn_dduration_y) +
                     f.z() * (dt2_dduration_z + mu_ * dn_dduration_z))) +
        (decision_->GetJacobianWrtNodes(t, kPos).row(X) *
         static_cast<double>(f.transpose() * (t2 + mu_ * n))); // t2 > -mu*n
  }

}

} /* namespace towr */
