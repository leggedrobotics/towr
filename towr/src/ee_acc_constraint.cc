/*
 * ee_acc_constraint.cc
 *
 *  Created on: Jan 16, 2020
 *      Author: vivian
 */

#include <towr/constraints/ee_acc_constraint.h>

#include <towr/variables/variable_names.h>
#include <towr/variables/cartesian_dimensions.h>

#include <iostream>

namespace towr {

EEAccConstraint::EEAccConstraint (const NodeSpline::Ptr& spline,
                                  std::string node_variable_name,
								  const std::vector<int>& dimensions)
    :ConstraintSet(kSpecifyLater, "acc-" + node_variable_name)
{
  spline_ = spline;
  ee_motion_id_ = node_variable_name;

  n_junctions_ = spline->GetPolynomialCount()-1;
  T_ = spline->GetPolyDurations();

  dimensions_  = dimensions;
  n_dim_       = dimensions.size();

}

void
EEAccConstraint::InitVariableDependedQuantities (const VariablesPtr& x)
{
  ee_motion_ = x->GetComponent<NodesVariablesPhaseBased>(ee_motion_id_);

  // select stance nodes
  auto nodes = ee_motion_->GetNodes();
  for (int id=0; id<nodes.size(); ++id) {
    if (ee_motion_->AdjacentPolyIdsConstant(id)) {
    	stance_node_ids_.push_back(id);
    	stance_poly_ids_.push_back(id-1);
    	stance_poly_ids_.push_back(id);
    }
  }

  std::cout << "ee" << ee_motion_id_[10] << " ACC stance nodes: ";
  for (int id : stance_node_ids_)
	  std::cout << id << " ";
  std::cout << std::endl;

  std::cout << "ee" << ee_motion_id_[10] << " ACC stance polys: ";
  for (int id : stance_poly_ids_)
	  std::cout << id << " ";
  std::cout << std::endl;

  swing_node_ids_ = ee_motion_->GetIndicesOfNonConstantNodes();
  std::cout << "ee" << ee_motion_id_[10] << " ACC swing nodes: ";
  for (int id : swing_node_ids_)
	  std::cout << id << " ";
  std::cout << std::endl;

//  SetRows(n_junctions_ + swing_node_ids_.size());
//  SetRows(stance_node_ids_.size() + swing_node_ids_.size());
  SetRows(stance_node_ids_.size()*n_dim_);
//  SetRows(n_junctions_*n_dim_);

}

Eigen::VectorXd
EEAccConstraint::GetValues () const
{
  VectorXd g(GetRows());

  int row = 0;

//  // constraint in x-direction for all nodes (TODO: change for all stance nodes)
//  for (int j=0; j<n_junctions_; ++j) {
//	int p_prev = j; // id of previous polynomial
//	VectorXd acc_prev = spline_->GetPoint(p_prev, T_.at(p_prev)).a();
//
//	int p_next = j+1;
//	VectorXd acc_next = spline_->GetPoint(p_next, 0.0).a();
//
//    for (auto dim : dimensions_)
//    	g(row++) = acc_prev(dim) - acc_next(dim);
//  }

  // constraint in x-direction for stance nodes
  for (auto node : stance_node_ids_) {
	int p_prev = node-1;
	VectorXd acc_prev = spline_->GetPoint(p_prev, T_.at(p_prev)).a();

	int p_next = node;
	VectorXd acc_next = spline_->GetPoint(p_next, 0.0).a();

    for (auto dim : dimensions_)
    	g(row++) = acc_prev(dim) - acc_next(dim);
  }

//  // constraint in z-direction for swing nodes
//  for (auto node : swing_node_ids_) {
//	int p_prev = node-1;
//	VectorXd acc_prev = spline_->GetPoint(p_prev, T_.at(p_prev)).a();
//
//	int p_next = node;
//	VectorXd acc_next = spline_->GetPoint(p_next, 0.0).a();
//
//	g(row++) = acc_prev(Z) - acc_next(Z);
//  }

  return g;
}

void
EEAccConstraint::FillJacobianBlock (std::string var_set, Jacobian& jac) const
{
  if (var_set == ee_motion_id_) {
	int row = 0;
//	for (int j=0; j<n_junctions_; ++j) {
//	  int p_prev = j; // id of previous polynomial
//	  Jacobian acc_prev = spline_->GetJacobianWrtNodes(p_prev, T_.at(p_prev), kAcc);
//
//	  int p_next = j+1;
//      Jacobian acc_next = spline_->GetJacobianWrtNodes(p_next, 0.0, kAcc);
//
//      for (auto dim : dimensions_)
//      	jac.row(row++) = acc_prev.row(dim) - acc_next.row(dim);
//	}

	// constraint in x-direction for stance nodes
	for (auto node : stance_node_ids_) {
	  int p_prev = node-1;
	  Jacobian acc_prev = spline_->GetJacobianWrtNodes(p_prev, T_.at(p_prev), kAcc);

	  int p_next = node;
	  Jacobian acc_next = spline_->GetJacobianWrtNodes(p_next, 0.0, kAcc);

      for (auto dim : dimensions_)
      	jac.row(row++) = acc_prev.row(dim) - acc_next.row(dim);
	}

//	for (auto node : swing_node_ids_) {
//	  int p_prev = node-1;
//	  Jacobian acc_prev = spline_->GetJacobianWrtNodes(p_prev, T_.at(p_prev), kAcc);
//
//	  int p_next = node;
//	  Jacobian acc_next = spline_->GetJacobianWrtNodes(p_next, 0.0, kAcc);
//
//	  jac.row(row++) = acc_prev.row(Z) - acc_next.row(Z);
//	}

  }
}

EEAccConstraint::VecBound
EEAccConstraint::GetBounds () const
{
  return VecBound(GetRows(), ifopt::BoundZero);
}

} /* namespace towr */



