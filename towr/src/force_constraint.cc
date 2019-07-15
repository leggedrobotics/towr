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

#include <towr/constraints/force_constraint.h>

#include <towr/variables/variable_names.h>

#include <towr/constraints/swing_constraint.h>
#include <towr/constraints/base_motion_constraint.h>
#include <towr/variables/cartesian_dimensions.h>
#include <iostream>
#include <math.h>

using namespace std;

namespace towr {

ForceConstraint::ForceConstraint (const HeightMap::Ptr& terrain,
                                  double force_limit,
                                  EE ee,
								  const SplineHolder& spline_holder)
    :ifopt::ConstraintSet(kSpecifyLater, "force-" + id::EEForceNodes(ee))
{
  terrain_ = terrain;
  fn_max_  = force_limit;
  mu_      = terrain->GetFrictionCoeff();
  ee_      = ee;
  base_angular_ = EulerConverter(spline_holder.base_angular_);

  n_constraints_per_node_ = 1 + 2*k2D;
  n_constraints_drift_node_ = 2;

}

void
ForceConstraint::InitVariableDependedQuantities (const VariablesPtr& x)
{
  ee_force_  = x->GetComponent<NodesVariablesPhaseBased>(id::EEForceNodes(ee_));
  ee_motion_ = x->GetComponent<NodesVariablesPhaseBased>(id::EEMotionNodes(ee_));

  int n_nodes_tot = params_.polynomials_in_first_drive_phase_ + params_.polynomials_in_second_drive_phase_ + params_.polynomials_in_drift_phase_ + 1;
  for (int i = 0; i<n_nodes_tot; i++){
  		node_ids_.push_back(i);
  }

  int constraint_count = 0;

  if (params_.just_drive_){
  	constraint_count = node_ids_.size()*n_constraints_per_node_;
  }
  else {
  		if (ee_ == 0 or ee_ == 1){
  		  	constraint_count = node_ids_.size()*n_constraints_per_node_;
  		}
  		if (ee_ == 2 or ee_ == 3){
  		  	constraint_count = (params_.polynomials_in_first_drive_phase_ + params_.polynomials_in_second_drive_phase_+1)*n_constraints_per_node_ + (params_.polynomials_in_drift_phase_)*n_constraints_drift_node_;
  		}
  }

  SetRows(constraint_count);

}

Eigen::VectorXd
ForceConstraint::GetValues () const
{
  VectorXd g(GetRows());

  int row=0;
  auto force_nodes = ee_force_->GetNodes();
  auto nodes = ee_motion_->GetNodes();
  for (int f_node_id : node_ids_) {
    int phase  = ee_force_->GetPhase(f_node_id, ee_);

    Vector3d p = nodes.at(f_node_id).p();

    Vector3d n = terrain_->GetNormalizedBasis(HeightMap::Normal, p.x(), p.y());
    Vector3d f = force_nodes.at(f_node_id).p();

    // unilateral force
    g(row++) = f.transpose() * n; // >0 (unilateral forces)

    Vector3d t1 = terrain_->GetNormalizedBasis(HeightMap::Tangent1, p.x(), p.y());
    Vector3d t2 = terrain_->GetNormalizedBasis(HeightMap::Tangent2, p.x(), p.y());

    // drive phase: friction rectangle constraint
    if (phase == 0){
    	g(row++) = f.transpose() * (t1 - mu_*n); // t1 < mu*n
    	g(row++) = f.transpose() * (t1 + mu_*n); // t1 > -mu*n
    	g(row++) = f.transpose() * (t2 - mu_*n); // t2 < mu*n
    	g(row++) = f.transpose() * (t2 + mu_*n); // t2 > -mu*n
    }

    double f_tang = sqrt(pow(f.transpose()*t1, 2) + pow(f.transpose()*t2, 2));

    // drift phase: friction circle constraint
    if (phase == 1){
    	g(row++) = f_tang - f.transpose() * (mu_*n);
    }
  }

  return g;
}

ForceConstraint::VecBound
ForceConstraint::GetBounds () const
{
  VecBound bounds;

  for (int f_node_id : node_ids_) {
	int phase  = ee_force_->GetPhase(f_node_id, ee_);
//	std::cout << "BOUNDS ee: " << ee_ << " Phase: " << phase << std::endl;

	if (phase == 0){
		bounds.push_back(ifopt::Bounds(0.0, fn_max_)); // unilateral forces
		bounds.push_back(ifopt::BoundSmallerZero); // f_t1 >  mu*n
		bounds.push_back(ifopt::BoundGreaterZero); // f_t1 < -mu*n
		bounds.push_back(ifopt::BoundSmallerZero); // f_t2 <  mu*n
		bounds.push_back(ifopt::BoundGreaterZero); // f_t2 > -mu*n
    }

    else if (phase == 1) {
    	bounds.push_back(ifopt::Bounds(0.0, fn_max_)); // unilateral forces
    	bounds.push_back(ifopt::BoundZero); // f_t1 =  mu*n
    }
  }

  return bounds;
}

void
ForceConstraint::FillJacobianBlock (std::string var_set,
                                    Jacobian& jac) const
{

  if (var_set == ee_force_->GetName()) {
      int row = 0;
      for (int f_node_id : node_ids_) {
        // unilateral force
    	auto nodes = ee_motion_->GetNodes();
    	auto force_nodes = ee_force_->GetNodes();
    	Vector3d p = nodes.at(f_node_id).p();
        int phase   = ee_force_->GetPhase(f_node_id, ee_);
        Vector3d n  = terrain_->GetNormalizedBasis(HeightMap::Normal,   p.x(), p.y());
        Vector3d t1 = terrain_->GetNormalizedBasis(HeightMap::Tangent1, p.x(), p.y());
        Vector3d t2 = terrain_->GetNormalizedBasis(HeightMap::Tangent2, p.x(), p.y());
        Vector3d f = force_nodes.at(f_node_id).p();

        for (auto dim : {X,Y,Z}) {
          int idx = ee_force_->GetOptIndex(NodesVariables::NodeValueInfo(f_node_id, kPos, dim));

          int row_reset=row;

          if (phase == 0){
        	  jac.coeffRef(row_reset++, idx) = n(dim);              // unilateral force
        	  jac.coeffRef(row_reset++, idx) = t1(dim)-mu_*n(dim);  // f_t1 <  mu*n
        	  jac.coeffRef(row_reset++, idx) = t1(dim)+mu_*n(dim);  // f_t1 > -mu*n
        	  jac.coeffRef(row_reset++, idx) = t2(dim)-mu_*n(dim);  // f_t2 <  mu*n
        	  jac.coeffRef(row_reset++, idx) = t2(dim)+mu_*n(dim);  // f_t2 > -mu*n
          }

          else if (phase == 1){
              jac.coeffRef(row_reset++, idx) = n(dim);              // unilateral force
              if (dim == X or dim == Y){
            	  if (f.transpose()(0) == 0 and f.transpose()(1) == 0)
                    	jac.coeffRef(row_reset++, idx) = 0;
                  else
                    	jac.coeffRef(row_reset++, idx) = f.transpose()(dim)/sqrt(pow(f.transpose()(0),2) + pow(f.transpose()(1),2)) - mu_*n(dim);
              }
              if (dim == Z) {
                  jac.coeffRef(row_reset++, idx) = -mu_;
              }
          }
        }

        if (phase == 1)
        	row += n_constraints_drift_node_;
        else
        	row += n_constraints_per_node_;
      }
    }

    if (var_set == ee_motion_->GetName()) {
      int row = 0;
      auto force_nodes = ee_force_->GetNodes();
      auto nodes = ee_motion_->GetNodes();
      for (int f_node_id : node_ids_) {
        int phase  = ee_force_->GetPhase(f_node_id, ee_);

        Vector3d p = nodes.at(f_node_id).p();
        Vector3d f = force_nodes.at(f_node_id).p();

        for (auto dim : {X_,Y_}) {
          Vector3d dn  = terrain_->GetDerivativeOfNormalizedBasisWrt(HeightMap::Normal, dim, p.x(), p.y());
          Vector3d dt1 = terrain_->GetDerivativeOfNormalizedBasisWrt(HeightMap::Tangent1, dim, p.x(), p.y());
          Vector3d dt2 = terrain_->GetDerivativeOfNormalizedBasisWrt(HeightMap::Tangent2, dim, p.x(), p.y());

          int idx = ee_motion_->GetOptIndex(NodesVariables::NodeValueInfo(f_node_id, kPos, dim));
          int row_reset=row;

          // unilateral force
          jac.coeffRef(row_reset++, idx) = f.transpose()*dn;

          if (phase == 0){
        	  jac.coeffRef(row_reset++, idx) = f.transpose()*(dt1-mu_*dn);
        	  jac.coeffRef(row_reset++, idx) = f.transpose()*(dt1+mu_*dn);
        	  jac.coeffRef(row_reset++, idx) = f.transpose()*(dt2-mu_*dn);
        	  jac.coeffRef(row_reset++, idx) = f.transpose()*(dt2+mu_*dn);
          }
          if (phase == 1){
        	  jac.coeffRef(row_reset++, idx) = 0;
          }
        }
        if (phase == 1)
                row += n_constraints_drift_node_;
        else
        		row += n_constraints_per_node_;
      }
  }
}

} /* namespace towr */
