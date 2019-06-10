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

#include <towr/constraints/spline_acc_constraint.h>

namespace towr {

SplineAccConstraint::SplineAccConstraint (const NodeSpline::Ptr& spline,
                                          std::string node_variable_name)
    :ConstraintSet(kSpecifyLater, "splineacc-" + node_variable_name)
{
  spline_ = spline;
  node_variables_id_ = node_variable_name;

  //n_dim_ anders für ee_splines als base splines?
  n_dim_       = spline->GetPoint(0.0).p().rows();
//  n_dim_ = 3;


  //TODO: different poly durations for ee and for base!!

  if (node_variables_id_ == "base-lin" or node_variables_id_ == "base-ang"){
  	  T_ = spline->GetPolyDurations();
  	  n_junctions_ = spline->GetPolynomialCount() - 1;
  }

  else {
	  n_junctions_ = (params_.polynomials_in_first_drive_phase_ + params_.polynomials_in_drift_phase_ + params_.polynomials_in_second_drive_phase_) - 1;

	  for (int i = 0; i < (n_junctions_ + 1); i++){
		  T_.push_back(params_.duration_node_polynomial_);
	  }
  }


  SetRows(n_dim_*n_junctions_);
}

Eigen::VectorXd
SplineAccConstraint::GetValues () const
{
  VectorXd g(GetRows());

  //for forces, "velocity" needs to be constraint:
  if (node_variables_id_ == "ee-force_0" or node_variables_id_ == "ee-force_1" or node_variables_id_ == "ee-force_2" or node_variables_id_ == "ee-force_3"){
	  for (int j=0; j<n_junctions_; ++j) {
	      int p_prev = j; // id of previous polynomial
	      VectorXd vel_prev = spline_->GetPoint(p_prev, T_.at(p_prev)).v();

	      int p_next = j+1;
	      VectorXd vel_next = spline_->GetPoint(p_next, 0.0).v();

	      g.segment(j*n_dim_, n_dim_) = vel_prev - vel_next;
	    }
  }

  else {
  for (int j=0; j<n_junctions_; ++j) {
    int p_prev = j; // id of previous polynomial
    VectorXd acc_prev = spline_->GetPoint(p_prev, T_.at(p_prev)).a();

    int p_next = j+1;
    VectorXd acc_next = spline_->GetPoint(p_next, 0.0).a();

    g.segment(j*n_dim_, n_dim_) = acc_prev - acc_next;
  }
  }

  return g;
}

void
SplineAccConstraint::FillJacobianBlock (std::string var_set, Jacobian& jac) const
{
  if (var_set == node_variables_id_) {
	  if (node_variables_id_ == "ee-force_0" or node_variables_id_ == "ee-force_1" or node_variables_id_ == "ee-force_2" or node_variables_id_ == "ee-force_3"){
		  for (int j=0; j<n_junctions_; ++j) {
		        int p_prev = j; // id of previous polynomial
		        Jacobian vel_prev = spline_->GetJacobianWrtNodes(p_prev, T_.at(p_prev), kVel, false);

		        int p_next = j+1;
		        Jacobian vel_next = spline_->GetJacobianWrtNodes(p_next, 0.0, kVel, false);

		        jac.middleRows(j*n_dim_, n_dim_) = vel_prev - vel_next;
		  }
	  }
	else {
    for (int j=0; j<n_junctions_; ++j) {
      int p_prev = j; // id of previous polynomial
      Jacobian acc_prev = spline_->GetJacobianWrtNodes(p_prev, T_.at(p_prev), kAcc, false);

      int p_next = j+1;
      Jacobian acc_next = spline_->GetJacobianWrtNodes(p_next, 0.0, kAcc, false);

      jac.middleRows(j*n_dim_, n_dim_) = acc_prev - acc_next;
    }
	}
  }
}

SplineAccConstraint::VecBound
SplineAccConstraint::GetBounds () const
{
  return VecBound(GetRows(), ifopt::BoundZero);
}

} /* namespace xpp */
