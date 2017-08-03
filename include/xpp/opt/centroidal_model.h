/**
 @file    centroidal_model.h
 @author  Alexander W. Winkler (winklera@ethz.ch)
 @date    May 19, 2017
 @brief   Brief description
 */

#ifndef XPP_OPT_INCLUDE_XPP_OPT_CENTROIDAL_MODEL_H_
#define XPP_OPT_INCLUDE_XPP_OPT_CENTROIDAL_MODEL_H_

#include "dynamic_model.h"

namespace xpp {
namespace opt {

/**
 * @brief Centroidal Dynamics for the 6-DoF Base to model the system.
 */
class CentroidalModel : public DynamicModel {
public:
  CentroidalModel (double mass, const Eigen::Matrix3d& inertia, int ee_count);
  virtual ~CentroidalModel ();

  virtual BaseAcc GetBaseAcceleration() const override;

  virtual Jacobian GetJacobianOfAccWrtBaseLin(const Jacobian& jac_base_lin_pos) const override;
  virtual Jacobian GetJacobianOfAccWrtBaseAng(const Jacobian& jac_base_ang_pos) const override;
  virtual Jacobian GetJacobianofAccWrtForce(const Jacobian& jac_force,
                                            EndeffectorID) const override;
  virtual Jacobian GetJacobianofAccWrtEEPos(const Jacobian& jac_ee_pos,
                                            EndeffectorID) const override;

private:
  double m_;          /// mass of robot
  Eigen::Matrix3d I_; /// inertia tensor of robot
  Eigen::SparseMatrix<double, Eigen::RowMajor> I_inv_;


  static Jacobian
  BuildCrossProductMatrix(const Vector3d& in)
  {
    Jacobian out(3,3);

    out.coeffRef(0,1) = -in(2); out.coeffRef(0,2) =  in(1);
    out.coeffRef(1,0) =  in(2); out.coeffRef(1,2) = -in(0);
    out.coeffRef(2,0) = -in(1); out.coeffRef(2,1) =  in(0);

    return out;
  }
};

} /* namespace opt */
} /* namespace xpp */

#endif /* XPP_OPT_INCLUDE_XPP_OPT_CENTROIDAL_MODEL_H_ */