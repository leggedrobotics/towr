/******************************************************************************
Copyright (c) 2018, Alexander W. Winkler. All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

* Neither the name of the copyright holder nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
******************************************************************************/

#include <towr/models/single_rigid_body_dynamics.h>
#include <towr/models/endeffector_mappings.h>

#include <towr/variables/cartesian_dimensions.h>

namespace towr {

// some Eigen helper functions
static Eigen::Matrix3d BuildInertiaTensor( double Ixx, double Iyy, double Izz,
                                           double Ixy, double Ixz, double Iyz)
{
  Eigen::Matrix3d I;
  I <<  Ixx, -Ixy, -Ixz,
       -Ixy,  Iyy, -Iyz,
       -Ixz, -Iyz,  Izz;
  return I;
}

// builds a cross product matrix out of "in", so in x v = X(in)*v
SingleRigidBodyDynamics::Jac
Cross(const Eigen::Vector3d& in)
{
  SingleRigidBodyDynamics::Jac out(3,3);

  out.coeffRef(0,1) = -in(2); out.coeffRef(0,2) =  in(1);
  out.coeffRef(1,0) =  in(2); out.coeffRef(1,2) = -in(0);
  out.coeffRef(2,0) = -in(1); out.coeffRef(2,1) =  in(0);

  return out;
}

SingleRigidBodyDynamics::SingleRigidBodyDynamics (double mass,
                                  double Ixx, double Iyy, double Izz,
                                  double Ixy, double Ixz, double Iyz,
                                  int ee_count)
   : SingleRigidBodyDynamics(mass,
                     BuildInertiaTensor(Ixx, Iyy, Izz, Ixy, Ixz, Iyz),
                     ee_count)
{
}

SingleRigidBodyDynamics::SingleRigidBodyDynamics (double mass, const Eigen::Matrix3d& inertia_b,
                                  int ee_count)
    :DynamicModel(mass, ee_count)
{
  I_b = inertia_b.sparseView();
}

SingleRigidBodyDynamics::BaseAcc
SingleRigidBodyDynamics::GetDynamicViolation () const
{
  // https://en.wikipedia.org/wiki/Newton%E2%80%93Euler_equations

  Vector3d f_sum, tau_sum;
  f_sum.setZero(); tau_sum.setZero();

  for (int ee=0; ee<ee_pos_.size(); ++ee) {
    Vector3d f = ee_force_.at(ee);
    tau_sum += f.cross(com_pos_ - ee_pos_.at(ee));
    f_sum   += f;
  }

  // express inertia matrix in world frame based on current body orientation
  Jac I_w = w_R_b_.sparseView() * I_b * w_R_b_.transpose().sparseView();

  BaseAcc acc;
  acc.segment(AX, k3D) = I_w*omega_dot_
                         + Cross(omega_)*(I_w*omega_)
                         - tau_sum;
  acc.segment(LX, k3D) = m()*com_acc_
                         - f_sum
                         - Vector3d(0.0, 0.0, -m()*g()); // gravity force
  return acc;
}

SingleRigidBodyDynamics::Jac
SingleRigidBodyDynamics::GetJacobianWrtBaseLin (const Jac& jac_pos_base_lin,
                                        		const Jac& jac_acc_base_lin) const
{
  // build the com jacobian
  int n = jac_pos_base_lin.cols();

  Jac jac_tau_sum(k3D, n);
  for (const Vector3d& f : ee_force_) {
    Jac jac_tau = Cross(f)*jac_pos_base_lin;
    jac_tau_sum += jac_tau;
  }

  Jac jac(k6D, n);
  jac.middleRows(AX, k3D) = -jac_tau_sum;
  jac.middleRows(LX, k3D) = m()*jac_acc_base_lin;

  return jac;
}

SingleRigidBodyDynamics::Jac
SingleRigidBodyDynamics::GetJacobianWrtBaseAng (const EulerConverter& base_euler,
                                        		double t) const
{
  Jac I_w = w_R_b_.sparseView() * I_b * w_R_b_.transpose().sparseView();

  // Derivative of R*I_b*R^T * wd
  // 1st term of product rule (derivative of R)
  Vector3d v11 = I_b*w_R_b_.transpose()*omega_dot_;
  Jac jac11 = base_euler.DerivOfRotVecMult(t, v11, false);

  // 2nd term of product rule (derivative of R^T)
  Jac jac12 = w_R_b_.sparseView()*I_b*base_euler.DerivOfRotVecMult(t, omega_dot_, true);

  // 3rd term of product rule (derivative of wd)
  Jac jac_ang_acc = base_euler.GetDerivOfAngAccWrtEulerNodes(t);
  Jac jac13 = I_w * jac_ang_acc;
  Jac jac1 = jac11 + jac12 + jac13;


  // Derivative of w x Iw
  // w x d_dn(R*I_b*R^T*w) -(I*w x d_dnw)
  // right derivative same as above, just with velocity instead acceleration
  Vector3d v21 = I_b*w_R_b_.transpose()*omega_;
  Jac jac21 = base_euler.DerivOfRotVecMult(t, v21, false);

  // 2nd term of product rule (derivative of R^T)
  Jac jac22 = w_R_b_.sparseView()*I_b*base_euler.DerivOfRotVecMult(t, omega_, true);

  // 3rd term of product rule (derivative of omega)
  Jac jac_ang_vel = base_euler.GetDerivOfAngVelWrtEulerNodes(t);
  Jac jac23 = I_w * jac_ang_vel;

  Jac jac2 = Cross(omega_)*(jac21+jac22+jac23) - Cross(I_w*omega_)*jac_ang_vel;


  // Combine the two to get sensitivity to I_w*w + w x (I_w*w)
  int n = jac_ang_vel.cols();
  Jac jac(k6D, n);
  jac.middleRows(AX, k3D) = jac1 + jac2;

  return jac;
}

SingleRigidBodyDynamics::Jac
SingleRigidBodyDynamics::GetJacobianWrtForce (const Jac& jac_force, EE ee) const
{
  Vector3d r = com_pos_ - ee_pos_.at(ee);
  Jac jac_tau = -Cross(r)*jac_force;

  int n = jac_force.cols();
  Jac jac(k6D, n);
  jac.middleRows(AX, k3D) = -jac_tau;
  jac.middleRows(LX, k3D) = -jac_force;

  return jac;
}

SingleRigidBodyDynamics::Jac
SingleRigidBodyDynamics::GetJacobianWrtEEPos (const Jac& jac_ee_pos, EE ee) const
{
  Vector3d f = ee_force_.at(ee);
  Jac jac_tau = Cross(f)*(-jac_ee_pos);

  Jac jac(k6D, jac_tau.cols());
  jac.middleRows(AX, k3D) = -jac_tau;

  // linear dynamics don't depend on endeffector position.
  return jac;
}

// right now this method only works for quadrupedal robots
//std::vector<double>
std::pair<int, double>
SingleRigidBodyDynamics::GetStabilityMeasure () const
{
  // https://en.wikipedia.org/wiki/Newton%E2%80%93Euler_equations

  int n_ee = ee_pos_.size();

  EEPos ee_pos_cw (n_ee+1); // re-ordered clockwise from above
  ee_pos_cw.at(0) = ee_pos_.at(LF);
  ee_pos_cw.at(1) = ee_pos_.at(RF);
  ee_pos_cw.at(2) = ee_pos_.at(RH);
  ee_pos_cw.at(3) = ee_pos_.at(LH);
  ee_pos_cw.at(4) = ee_pos_.at(LF);

  Jac I_w = w_R_b_.sparseView() * I_b * w_R_b_.transpose().sparseView();

  // compute tipover axes
  std::vector<double> theta(n_ee);  // tipover angles
  std::vector<double> beta(n_ee);   // force-angle stability measure

  for (int ee=0; ee<n_ee; ++ee)
  {
	Vector3d a = (ee_pos_cw.at(ee+1) - ee_pos_cw.at(ee));

	Matrix3d proj_a = Matrix3d::Identity() - a.normalized() * a.normalized().transpose();
	Vector3d l = proj_a * (ee_pos_cw.at(ee+1) - com_pos_);
	Vector3d l_n = l.normalized();

    Vector3d f = proj_a * (m()*com_acc_ + Vector3d(0.0, 0.0, -m()*g()));

    Vector3d n = (Matrix3d::Identity() - proj_a) * (I_w*omega_dot_ + Cross(omega_)*(I_w*omega_));

    // I am neglecting the contributions from angular loads because they considerable small and would make defining
    // the Jacobian analytically far more difficult!
    Vector3d ft = f; //+ l_n.cross(n)/l.norm();
    Vector3d ft_n = ft.normalized();

    Vector3d d = -l + (l.transpose()*ft_n)*ft_n;

    double tmp = ft_n.cross(l_n).transpose() * a.normalized();
    int sgn = tmp > 0 ? 1 : -1;
    theta.at(ee) = sgn * std::acos(ft_n.transpose()*l_n);
    beta.at(ee) = theta.at(ee) * d.norm() * ft.norm();
  }

  auto it = std::min_element(std::begin(theta), std::end(theta));
  auto min_idx = std::distance(std::begin(theta), it);

  return std::make_pair((int) min_idx, *it);

  //return theta;
}

SingleRigidBodyDynamics::Vector3d
SingleRigidBodyDynamics::GetStabilityDerivativeWrtBaseLinPos () const
{
  std::pair<int, double> stab = GetStabilityMeasure();
  int tipover_axis = stab.first;
  double measure = stab.second;
  std::vector<int> ee = GetEndEffectorsAdjacentTo(tipover_axis);

  Vector3d a = ee_pos_.at(ee.at(1)) - ee_pos_.at(ee.at(0));
  Matrix3d proj_a = Matrix3d::Identity() - a.normalized() * a.normalized().transpose();
  Vector3d l = proj_a * (ee_pos_.at(ee.at(1)) - com_pos_);
  Vector3d f = proj_a * (m()*com_acc_ + Vector3d(0.0, 0.0, -m()*g()));
  double arg = f.normalized().transpose()*l.normalized();
  int sgn = measure > 0 ? 1 : -1;
  double acos_deriv = GetDerivativeOfAcos(arg);

  Matrix3d diff_lnorm_pos = GetDerivativeOfNormalizedVector(l) * (-proj_a);
  Vector3d d_stab_pos = sgn*acos_deriv * ( f.normalized().transpose()*diff_lnorm_pos ).transpose() + Vector3d::Constant(1e-10);

  return d_stab_pos;
}

SingleRigidBodyDynamics::Vector3d
SingleRigidBodyDynamics::GetStabilityDerivativeWrtBaseLinAcc () const
{
  std::pair<int, double> stab = GetStabilityMeasure();
  int tipover_axis = stab.first;
  double measure = stab.second;
  std::vector<int> ee = GetEndEffectorsAdjacentTo(tipover_axis);

  Vector3d a = ee_pos_.at(ee.at(1)) - ee_pos_.at(ee.at(0));
  Matrix3d proj_a = Matrix3d::Identity() - a.normalized() * a.normalized().transpose();
  Vector3d l = proj_a * (ee_pos_.at(ee.at(1)) - com_pos_);
  Vector3d f = proj_a * (m()*com_acc_ + Vector3d(0.0, 0.0, -m()*g()));
  double arg = f.normalized().transpose()*l.normalized();
  int sgn = measure > 0 ? 1 : -1;
  double acos_deriv = GetDerivativeOfAcos(arg);

  Matrix3d diff_fnorm_acc = GetDerivativeOfNormalizedVector(f) * (m()*proj_a);

  Vector3d d_stab_acc = sgn*acos_deriv * ( diff_fnorm_acc.transpose()*l.normalized() ) + Vector3d::Constant(1e-10);

  return d_stab_acc;
}

SingleRigidBodyDynamics::MatrixXd
SingleRigidBodyDynamics::GetStabilityDerivativeWrtEEPos () const
{
  std::pair<int, double> stab = GetStabilityMeasure();
  int tipover_axis = stab.first;
  double measure = stab.second;
  std::vector<int> ee = GetEndEffectorsAdjacentTo(tipover_axis);

  Vector3d a = ee_pos_.at(ee.at(1)) - ee_pos_.at(ee.at(0));
  Matrix3d proj_a = Matrix3d::Identity() - a.normalized() * a.normalized().transpose();
  Vector3d l = proj_a * (ee_pos_.at(ee.at(1)) - com_pos_);
  Vector3d f = proj_a * (m()*com_acc_ + Vector3d(0.0, 0.0, -m()*g()));
  double arg = f.normalized().transpose()*l.normalized();
  int sgn = measure > 0 ? 1 : -1;
  double acos_deriv = GetDerivativeOfAcos(arg);

  // derivative wrt ee_pos_
  Vector3d diff_ee0;
  Vector3d diff_ee1;
  for (int dim : {X, Y, Z})
  {
	Vector3d diff_fnorm_eepos = GetDerivativeOfNormalizedVector(f) *
								GetDerivativeOfProjectionMatrixWrtVector(a,-Vector3d::Unit(dim)) *
								(m()*com_acc_ + Vector3d(0.0, 0.0, -m()*g()));

	Vector3d diff_lnorm_eepos = GetDerivativeOfNormalizedVector(l) *
								GetDerivativeOfProjectionMatrixWrtVector(a,-Vector3d::Unit(dim)) *
								(ee_pos_.at(ee.at(1)) - com_pos_);

	double d1 = diff_fnorm_eepos.transpose()*l.normalized();
	double d2 = f.normalized().transpose()*diff_lnorm_eepos;
	diff_ee0(dim) = d1 + d2;

	diff_fnorm_eepos = GetDerivativeOfNormalizedVector(f) *
					   GetDerivativeOfProjectionMatrixWrtVector(a,Vector3d::Unit(dim)) *
					   (m()*com_acc_ + Vector3d(0.0, 0.0, -m()*g()));

	diff_lnorm_eepos = GetDerivativeOfNormalizedVector(l) *
					   ( GetDerivativeOfProjectionMatrixWrtVector(a,Vector3d::Unit(dim))*(ee_pos_.at(ee.at(1))-com_pos_) +
					   proj_a*Vector3d::Unit(dim) );

	d1 = diff_fnorm_eepos.transpose()*l.normalized();
	d2 = f.normalized().transpose()*diff_lnorm_eepos;
	diff_ee1(dim) = d1 + d2;

  }

  MatrixXd d_stab_ee = MatrixXd::Zero(3,4);
  d_stab_ee.col(ee.at(0)) = sgn*acos_deriv*diff_ee0 + Vector3d::Constant(1e-10);
  d_stab_ee.col(ee.at(1)) = sgn*acos_deriv*diff_ee1 + Vector3d::Constant(1e-10);

  return d_stab_ee;
}

SingleRigidBodyDynamics::Matrix3d
SingleRigidBodyDynamics::GetDerivativeOfProjectionMatrixWrtVector (const Vector3d& v, const Vector3d& dv) const
{
  Matrix3d dP = (dv*v.transpose() + v*dv.transpose())*v.squaredNorm() -
		        (v*v.transpose())*2*(v.transpose()*dv);

  return -dP/(v.squaredNorm()*v.squaredNorm());
}

std::vector<int>
SingleRigidBodyDynamics::GetEndEffectorsAdjacentTo (int tipover_axis) const
{
  std::vector<int> ee;

  switch (tipover_axis) {
	case 0: ee = {0, 1}; break;
	case 1: ee = {1, 3}; break;
	case 2: ee = {3, 2}; break;
	case 3: ee = {2, 0}; break;
	default: assert(false); break;
  }

  return ee;
}

double
SingleRigidBodyDynamics::GetDerivativeOfAcos (const double x) const
{
  return -1.0/(sqrt(1-pow(x,2)));
}

SingleRigidBodyDynamics::Matrix3d
SingleRigidBodyDynamics::GetDerivativeOfNormalizedVector (const Vector3d& v) const
{
  // see http://blog.mmacklin.com/2012/05/
  double aux = 1/v.squaredNorm();

  return aux * (v.norm() * Matrix3d::Identity() - v*(v.normalized().transpose()));
}

} /* namespace towr */
