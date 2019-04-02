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
#include <towr/variables/cartesian_dimensions.h>
#include <iostream>

using namespace std;

namespace towr {

ForceConstraint::ForceConstraint (const HeightMap::Ptr& terrain,
                                  double force_limit,
                                  EE ee)
    :ifopt::ConstraintSet(kSpecifyLater, "force-" + id::EEForceNodes(ee))
{
  terrain_ = terrain;
  fn_max_  = force_limit;
  mu_      = terrain->GetFrictionCoeff();
  ee_      = ee;

  n_constraints_per_node_ = 1 + 2*k2D; // positive normal force + 4 friction pyramid constraints

  //new: plus pos and vel of every node: Node::n_derivatives*k2D
//  n_constraints_per_node_ = 1 + 2*k2D + Node::n_derivatives*k2D;
//  n_constraints_per_node_ = 1 + 2*k2D + 1;

//  cout << "n_constraints_per_node: " << n_constraints_per_node_ << endl;
//  cout << "ee_: " << ee_ << endl;

}

//SwingConstraint::SwingConstraint (std::string ee_motion)
//    :ConstraintSet(kSpecifyLater, "swing-" + ee_motion)
//{
//  ee_motion_id_ = ee_motion;
//}

void
ForceConstraint::InitVariableDependedQuantities (const VariablesPtr& x)
{
  ee_force_  = x->GetComponent<NodesVariablesPhaseBased>(id::EEForceNodes(ee_));
  ee_motion_ = x->GetComponent<NodesVariablesPhaseBased>(id::EEMotionNodes(ee_));
//  ee_motion_ = x->GetComponent<NodesVariablesPhaseBased>(ee_motion_id_);

  //take all nodes because we have pure driving:
  	pure_stance_force_node_ids_ = ee_force_->GetIndicesOfAllNodes();
//  	stance_node_ids_ = {0,2,4,6};
//  	motion_node_ids_ = {1,3,5,7};
//  pure_stance_force_node_ids_ = ee_motion_->GetIndicesOfAllNodes();
//  pure_stance_force_node_ids_ = ee_force_->GetIndicesOfNonConstantNodes();
//  pure_stance_force_node_ids_.push_back(ee_force_->GetIndicesOfConstantNodes());

//  	pure_stance_force_node_ids_ = {0,1,2,3};

  int constraint_count = pure_stance_force_node_ids_.size()*n_constraints_per_node_;
//  int constraint_count = 4*5 + 4*1; //4 force nodes und 4 motion nodes

  SetRows(constraint_count);
  cout << "number of TOTAL constraints (per node*nodes): " << constraint_count << endl;
//  for (int f_node_id : pure_stance_force_node_ids_) {
//	  cout << "node ids used: " << pure_stance_force_node_ids_[f_node_id] << ", ";
//  }
}

Eigen::VectorXd
ForceConstraint::GetValues () const
{
  VectorXd g(GetRows());

  int row=0;
  auto force_nodes = ee_force_->GetNodes();
  auto nodes = ee_motion_->GetNodes();	//new
  for (int f_node_id : pure_stance_force_node_ids_) {
//  for (int f_node_id : stance_node_ids_) {
    int phase  = ee_force_->GetPhase(f_node_id);
//    Vector3d p = ee_motion_->GetValueAtStartOfPhase(phase); // doesn't change during stance phase

//    if (f_node_id == 0 or f_node_id == 2 or f_node_id == 4 or f_node_id == 6){
    //new
    Vector3d p = nodes.at(f_node_id).p(); //nicht alle nodes existieren..!

    Vector3d n = terrain_->GetNormalizedBasis(HeightMap::Normal, p.x(), p.y());
    Vector3d f = force_nodes.at(f_node_id).p();

    // unilateral force
    g(row++) = f.transpose() * n; // >0 (unilateral forces)

    // frictional pyramid
    //new force restriction in t1 direction, so that movement allowed!
    Vector3d t1 = terrain_->GetNormalizedBasis(HeightMap::Tangent1, p.x(), p.y());
    g(row++) = f.transpose() * (t1 - mu_*n); // t1 < mu*n
    g(row++) = f.transpose() * (t1 + mu_*n); // t1 > -mu*n

    Vector3d t2 = terrain_->GetNormalizedBasis(HeightMap::Tangent2, p.x(), p.y());
    g(row++) = f.transpose() * (t2 - mu_*n); // t2 < mu*n
    g(row++) = f.transpose() * (t2 + mu_*n); // t2 > -mu*n
//      g(row++) = f.transpose() * t2;
//    }

////    else {
////    for (int m_node_id : motion_node_ids_) {
//    //new
//    // assumes three splines per drive phase (and starting and ending in stance??)
//        auto curr = nodes.at(f_node_id);
//
//        //TODO no velocity in Y direction during drive!
////        Vector2d prev;
////        Vector2d next;
//
//        Vector2d prev = nodes.at(f_node_id-1).p().topRows<2>();	//1 for just x-direction
////        prev.row(1).setZero();
//
//        Vector2d next = nodes.at(f_node_id+1).p().topRows<2>();
////        next.row(1).setZero();
//
////        Vector2d ybound = nodes.at(f_node_id).p().topRows<2>();
//
//
////        des_vel_center.row(1).setZero();
////        des_vel_center(Y)=0;
////        xy_center(Y)=0;
//
////        for (auto dim : {X,Y}) {
////          g(row++) = curr.p()(dim) - xy_center(dim);
////          g(row++) = curr.v()(dim) - des_vel_center(dim);
////        }
//        if (f_node_id == 0){
//        	Vector2d distance_xy    = next - curr.p().topRows<2>();
//        	        Vector2d xy_center      = curr.p().topRows<2>() + distance_xy;
//        	        Vector2d des_vel_center = distance_xy/(t_drive_/3); // linear interpolation not accurate
////        	        g(row++) = curr.p()(X) - xy_center(X);
//        	        g(row++) = curr.v()(X) - des_vel_center(X);
//        	//        g(row++) = curr.p()(Y) - prev(Y);
////        	        g(row++) = curr.v()(Y);
//        }
//        else {
//        	Vector2d distance_xy    = next - prev;
//        	Vector2d xy_center      = prev + distance_xy;
//        	Vector2d des_vel_center = distance_xy/(t_drive_/3); // linear interpolation not accurate
////        	g(row++) = curr.p()(X) - xy_center(X);
//        	g(row++) = curr.v()(X) - des_vel_center(X);
//        	//        g(row++) = curr.p()(Y) - prev(Y);
////        	g(row++) = curr.v()(Y);
//        }
//    }
  }

//  cout << "n_constraints: " << g;

  return g;
}

ForceConstraint::VecBound
ForceConstraint::GetBounds () const
{
  VecBound bounds;

  for (int f_node_id : pure_stance_force_node_ids_) {
//	  if (f_node_id == 0 or f_node_id == 2 or f_node_id == 4 or f_node_id == 6){
    bounds.push_back(ifopt::Bounds(0.0, fn_max_)); // unilateral forces
    bounds.push_back(ifopt::BoundSmallerZero); // f_t1 >  mu*n
    bounds.push_back(ifopt::BoundGreaterZero); // f_t1 < -mu*n
    bounds.push_back(ifopt::BoundSmallerZero); // f_t2 <  mu*n
    bounds.push_back(ifopt::BoundGreaterZero); // f_t2 > -mu*n
//	  }

//	  else {
//    bounds.push_back(ifopt::BoundZero);	//bound for x-velocity
//    bounds.push_back(ifopt::BoundZero);
//    bounds.push_back(ifopt::BoundZero); //bound for y-pos (=previous pos)
//    bounds.push_back(ifopt::BoundZero); //bound for y-vel
//	  }
  }

  return bounds;
}

void
ForceConstraint::FillJacobianBlock (std::string var_set,
                                    Jacobian& jac) const
{
//  if (var_set == ee_force_->GetName()) {
//    int row = 0;
//    auto nodes = ee_motion_->GetNodes();	//new
//    for (int f_node_id : pure_stance_force_node_ids_) {
////      // unilateral force
//      int phase   = ee_force_->GetPhase(f_node_id);
////      Vector3d p  = ee_motion_->GetValueAtStartOfPhase(phase); // doesn't change during phase
//      Vector3d p = nodes.at(f_node_id).p();
//      Vector3d n  = terrain_->GetNormalizedBasis(HeightMap::Normal,   p.x(), p.y());
//      Vector3d t1 = terrain_->GetNormalizedBasis(HeightMap::Tangent1, p.x(), p.y());
//      Vector3d t2 = terrain_->GetNormalizedBasis(HeightMap::Tangent2, p.x(), p.y());
////
//////      //new: added different constraints for x direction
//////            int idx = ee_force_->GetOptIndex(NodesVariables::NodeValueInfo(f_node_id, kPos, X));
//////
//////                    int row_reset=row;
//////
//////                    jac.coeffRef(row_reset++, idx) = n(X);              // unilateral force
//////                    jac.coeffRef(row_reset++, idx) = -t1(X)+mu_*n(X);  // f_t1 >  mu*n	//new
//////                    jac.coeffRef(row_reset++, idx) = -t1(X)-mu_*n(X);  // f_t1 < -mu*n	//new
//////                    jac.coeffRef(row_reset++, idx) = -t2(X)+mu_*n(X);  // f_t2 <  mu*n
//////                    jac.coeffRef(row_reset++, idx) = -t2(X)-mu_*n(X);  // f_t2 > -mu*n
////
//      for (auto dim : {X,Y,Z}) {
//        int idx = ee_force_->GetOptIndex(NodesVariables::NodeValueInfo(f_node_id, kPos, dim));
//
//        int row_reset=row;
//
//        jac.coeffRef(row_reset++, idx) = n(dim);              // unilateral force
//        jac.coeffRef(row_reset++, idx) = t1(dim)-mu_*n(dim);  // f_t1 <  mu*n	//new
//        jac.coeffRef(row_reset++, idx) = t1(dim)+mu_*n(dim);  // f_t1 > -mu*n	//new
//        jac.coeffRef(row_reset++, idx) = t2(dim)-mu_*n(dim);  // f_t2 <  mu*n
//        jac.coeffRef(row_reset++, idx) = t2(dim)+mu_*n(dim);  // f_t2 > -mu*n
//
//////        auto curr = nodes.at(f_node_id);
//////
//////                              //TODO no velocity in Y direction during drive!
//////                      //        Vector2d prev;
//////                      //        Vector2d next;
//////
//////                              Vector2d prev = nodes.at(f_node_id-1).p().topRows<2>();	//1 for just x-direction
//////                      //        prev.row(1).setZero();
//////
//////                              Vector2d next = nodes.at(f_node_id+1).p().topRows<2>();
//////                      //        next.row(1).setZero();
//////
//////                      //        Vector2d ybound = nodes.at(f_node_id).p().topRows<2>();
//////
//////                              Vector2d distance_xy    = next - prev;
//////                              Vector2d xy_center      = prev + 0.5*distance_xy;
//////                              Vector2d des_vel_center = distance_xy/t_swing_avg_;
//////        if (dim == X){
//////        	jac.coeffRef(row_reset++, idx) = curr.p()(X) - xy_center(X);
//////        	jac.coeffRef(row_reset++, idx) = curr.v()(X) - des_vel_center(X);
//////        	jac.coeffRef(row_reset++, idx) = 0;
//////        	jac.coeffRef(row_reset++, idx) = 0;
//////        }
//////
//////        if (dim == Y){
//////               	jac.coeffRef(row_reset++, idx) = 0;
//////               	jac.coeffRef(row_reset++, idx) = 0;
//////               	jac.coeffRef(row_reset++, idx) = curr.p()(Y) - prev(Y);
//////               	jac.coeffRef(row_reset++, idx) = curr.v()(Y);
//////               }
//////
//////        if (dim == Z){
//////                       	jac.coeffRef(row_reset++, idx) = 0;
//////                       	jac.coeffRef(row_reset++, idx) = 0;
//////                       	jac.coeffRef(row_reset++, idx) = 0;
//////                       	jac.coeffRef(row_reset++, idx) = 0;
//////                       }
//////                              if (dim == X){
//////
//////
//////                              // position constraint
//////                                                          jac.coeffRef(row_reset++, ee_motion_->GetOptIndex(NodesVariables::NodeValueInfo(f_node_id,   kPos, X_))) =  curr.p()(X) - xy_center(X);  // current node
//////                              //                            jac.coeffRef(row_reset, ee_motion_->GetOptIndex(NodesVariables::NodeValueInfo(f_node_id+1, kPos, X_))) = -0.5;  // next node
//////                              //                            jac.coeffRef(row_reset, ee_motion_->GetOptIndex(NodesVariables::NodeValueInfo(f_node_id-1, kPos, X_))) = -0.5;  // previous node
//////                              //                            row_reset++;
//////
//////                                            // velocity constraint
//////                                                          jac.coeffRef(row_reset++, ee_motion_->GetOptIndex(NodesVariables::NodeValueInfo(f_node_id,   kVel, X_))) =  curr.v()(X) - des_vel_center(X);              // current node
//////                              //                            jac.coeffRef(row_reset, ee_motion_->GetOptIndex(NodesVariables::NodeValueInfo(f_node_id+1, kPos, X_))) = -1.0/t_swing_avg_; // next node
//////                              //                            jac.coeffRef(row_reset, ee_motion_->GetOptIndex(NodesVariables::NodeValueInfo(f_node_id-1, kPos, X_))) = +1.0/t_swing_avg_; // previous node
//////                              //                            row_reset++;
//////
//////                                                          //pos und vel in y rtg bei "X constraint" muss 0 sein
//////                                                          jac.coeffRef(row_reset++, ee_motion_->GetOptIndex(NodesVariables::NodeValueInfo(f_node_id,   kPos, X_))) =  0;
//////                                                          jac.coeffRef(row_reset++, ee_motion_->GetOptIndex(NodesVariables::NodeValueInfo(f_node_id,   kVel, X_))) =  0;
//////                              }
//////                              if (dim == Y){
//////
//////
//////                                                            // position constraint
//////                                                                                        jac.coeffRef(row_reset++, ee_motion_->GetOptIndex(NodesVariables::NodeValueInfo(f_node_id,   kPos, Y_))) =  0;  // current node
//////                                                            //                            jac.coeffRef(row_reset, ee_motion_->GetOptIndex(NodesVariables::NodeValueInfo(f_node_id+1, kPos, X_))) = -0.5;  // next node
//////                                                            //                            jac.coeffRef(row_reset, ee_motion_->GetOptIndex(NodesVariables::NodeValueInfo(f_node_id-1, kPos, X_))) = -0.5;  // previous node
//////                                                            //                            row_reset++;
//////
//////                                                                          // velocity constraint
//////                                                                                        jac.coeffRef(row_reset++, ee_motion_->GetOptIndex(NodesVariables::NodeValueInfo(f_node_id,   kVel, Y_))) =  0;              // current node
//////                                                            //                            jac.coeffRef(row_reset, ee_motion_->GetOptIndex(NodesVariables::NodeValueInfo(f_node_id+1, kPos, X_))) = -1.0/t_swing_avg_; // next node
//////                                                            //                            jac.coeffRef(row_reset, ee_motion_->GetOptIndex(NodesVariables::NodeValueInfo(f_node_id-1, kPos, X_))) = +1.0/t_swing_avg_; // previous node
//////                                                            //                            row_reset++;
//////
//////                                                                                        //pos und vel in y rtg bei "X constraint" muss 0 sein
//////                                                                                        jac.coeffRef(row_reset++, ee_motion_->GetOptIndex(NodesVariables::NodeValueInfo(f_node_id,   kPos, Y_))) =  curr.p()(Y) - prev(Y);  // current node;
//////                                                                                        jac.coeffRef(row_reset++, ee_motion_->GetOptIndex(NodesVariables::NodeValueInfo(f_node_id,   kVel, Y_))) =  curr.v()(Y);
//////                                                            }
//////                              if (dim == Z){
//////
//////
//////                                                            // position constraint
//////                                                                                        jac.coeffRef(row_reset++, ee_motion_->GetOptIndex(NodesVariables::NodeValueInfo(f_node_id,   kPos, Z))) =  0;  // current node
//////                                                            //                            jac.coeffRef(row_reset, ee_motion_->GetOptIndex(NodesVariables::NodeValueInfo(f_node_id+1, kPos, X_))) = -0.5;  // next node
//////                                                            //                            jac.coeffRef(row_reset, ee_motion_->GetOptIndex(NodesVariables::NodeValueInfo(f_node_id-1, kPos, X_))) = -0.5;  // previous node
//////                                                            //                            row_reset++;
//////
//////                                                                          // velocity constraint
//////                                                                                        jac.coeffRef(row_reset++, ee_motion_->GetOptIndex(NodesVariables::NodeValueInfo(f_node_id,   kVel, Z))) =  0;              // current node
//////                                                            //                            jac.coeffRef(row_reset, ee_motion_->GetOptIndex(NodesVariables::NodeValueInfo(f_node_id+1, kPos, X_))) = -1.0/t_swing_avg_; // next node
//////                                                            //                            jac.coeffRef(row_reset, ee_motion_->GetOptIndex(NodesVariables::NodeValueInfo(f_node_id-1, kPos, X_))) = +1.0/t_swing_avg_; // previous node
//////                                                            //                            row_reset++;
//////
//////                                                                                        //pos und vel in y rtg bei "X constraint" muss 0 sein
//////                                                                                        jac.coeffRef(row_reset++, ee_motion_->GetOptIndex(NodesVariables::NodeValueInfo(f_node_id,   kPos, Z))) =  0;
//////                                                                                        jac.coeffRef(row_reset++, ee_motion_->GetOptIndex(NodesVariables::NodeValueInfo(f_node_id,   kVel, Z))) =  0;
//////                                                            }
//      }
//
//      row += n_constraints_per_node_;
//////        row += 5;
//    }
//////    cout << "n_rows_force_Jacobian: " << row;
//////    n_rows_force_Jacobian = row;
//      }
////
//  if (var_set == ee_motion_->GetName()) {
//    int row = 0;
//    auto force_nodes = ee_force_->GetNodes();
//    auto nodes = ee_motion_->GetNodes();	//new
//    for (int f_node_id : pure_stance_force_node_ids_) {
//      int phase  = ee_force_->GetPhase(f_node_id);
//      int ee_node_id = f_node_id;		//ee_motion_->GetNodeIDAtStartOfPhase(phase);	//node ID during phase...?
//
////      Vector3d p = ee_motion_->GetValueAtStartOfPhase(phase); // doesn't change during phase
//      Vector3d p = nodes.at(f_node_id).p(); //new
//      Vector3d f = force_nodes.at(f_node_id).p();
//
//      //new: seperated for X and Y direction!
//      for (auto dim : {X_,Y_}) {
//        Vector3d dn  = terrain_->GetDerivativeOfNormalizedBasisWrt(HeightMap::Normal, dim, p.x(), p.y());
//        Vector3d dt1 = terrain_->GetDerivativeOfNormalizedBasisWrt(HeightMap::Tangent1, dim, p.x(), p.y());
//        Vector3d dt2 = terrain_->GetDerivativeOfNormalizedBasisWrt(HeightMap::Tangent2, dim, p.x(), p.y());
////
//////        //TODO what should ee_node_id be? (not only ID at Start of phase..?
//        int idx = ee_motion_->GetOptIndex(NodesVariables::NodeValueInfo(ee_node_id, kPos, dim));
//        int row_reset=row;
//
//        // unilateral force
//        jac.coeffRef(row_reset++, idx) = f.transpose()*dn;
//
//        // friction force tangent 1 derivative
//        jac.coeffRef(row_reset++, idx) = f.transpose()*(dt1-mu_*dn);
//        jac.coeffRef(row_reset++, idx) = f.transpose()*(dt1+mu_*dn);
//
//        // friction force tangent 2 derivative
//        jac.coeffRef(row_reset++, idx) = f.transpose()*(dt2-mu_*dn);
//        jac.coeffRef(row_reset++, idx) = f.transpose()*(dt2+mu_*dn);
//////
//////        //new
//////        // position constraint
//////                      jac.coeffRef(row_reset, ee_motion_->GetOptIndex(NodesVariables::NodeValueInfo(f_node_id,   kPos, dim))) =  1.0;  // current node
//////                      jac.coeffRef(row_reset, ee_motion_->GetOptIndex(NodesVariables::NodeValueInfo(f_node_id+1, kPos, dim))) = -0.5;  // next node
//////                      jac.coeffRef(row_reset, ee_motion_->GetOptIndex(NodesVariables::NodeValueInfo(f_node_id-1, kPos, dim))) = -0.5;  // previous node
//////                      row_reset++;
//////
//////        // velocity constraint
//////                      jac.coeffRef(row_reset, ee_motion_->GetOptIndex(NodesVariables::NodeValueInfo(f_node_id,   kVel, dim))) =  1.0;              // current node
//////                      jac.coeffRef(row_reset, ee_motion_->GetOptIndex(NodesVariables::NodeValueInfo(f_node_id+1, kPos, dim))) = -1.0/t_swing_avg_; // next node
//////                      jac.coeffRef(row_reset, ee_motion_->GetOptIndex(NodesVariables::NodeValueInfo(f_node_id-1, kPos, dim))) = +1.0/t_swing_avg_; // previous node
//////                      row_reset++;
//////      }
////
////              Vector3d dn  = terrain_->GetDerivativeOfNormalizedBasisWrt(HeightMap::Normal, X_, p.x(), p.y());
////              Vector3d dt1 = terrain_->GetDerivativeOfNormalizedBasisWrt(HeightMap::Tangent1, X_, p.x(), p.y());
////              Vector3d dt2 = terrain_->GetDerivativeOfNormalizedBasisWrt(HeightMap::Tangent2, X_, p.x(), p.y());
////
////              //TODO what should ee_node_id be? (not only ID at Start of phase..? same as f_node_id!
////              int idx1 = ee_motion_->GetOptIndex(NodesVariables::NodeValueInfo(ee_node_id, kPos, X_));
////              int row_reset=row;
////
////              // unilateral force
////              jac.coeffRef(row_reset++, idx1) = f.transpose()*dn;
////
////              // friction force tangent 1 derivative
////              jac.coeffRef(row_reset++, idx1) = f.transpose()*(dt1-mu_*dn);
////              jac.coeffRef(row_reset++, idx1) = f.transpose()*(dt1+mu_*dn);
////
////              // friction force tangent 2 derivative
////              jac.coeffRef(row_reset++, idx1) = f.transpose()*(dt2-mu_*dn);
////              jac.coeffRef(row_reset++, idx1) = f.transpose()*(dt2+mu_*dn);
//
//////              new
////              auto curr = nodes.at(f_node_id);
////              Vector2d prev = nodes.at(f_node_id-1).p().topRows<2>();
////              Vector2d next = nodes.at(f_node_id+1).p().topRows<2>();
////              Vector2d distance_xy    = next - prev;
////              Vector2d xy_center      = prev + 0.5*distance_xy;
////              Vector2d des_vel_center = distance_xy/t_swing_avg_;
////
////              // position constraint
////                            jac.coeffRef(row_reset++, ee_motion_->GetOptIndex(NodesVariables::NodeValueInfo(f_node_id,   kPos, X_))) =  curr.p()(X) - xy_center(X);  // current node
//////                            jac.coeffRef(row_reset, ee_motion_->GetOptIndex(NodesVariables::NodeValueInfo(f_node_id+1, kPos, X_))) = -0.5;  // next node
//////                            jac.coeffRef(row_reset, ee_motion_->GetOptIndex(NodesVariables::NodeValueInfo(f_node_id-1, kPos, X_))) = -0.5;  // previous node
//////                            row_reset++;
////
////              // velocity constraint
////                            jac.coeffRef(row_reset++, ee_motion_->GetOptIndex(NodesVariables::NodeValueInfo(f_node_id,   kVel, X_))) =  curr.v()(X) - des_vel_center(X);              // current node
//////                            jac.coeffRef(row_reset, ee_motion_->GetOptIndex(NodesVariables::NodeValueInfo(f_node_id+1, kPos, X_))) = -1.0/t_swing_avg_; // next node
//////                            jac.coeffRef(row_reset, ee_motion_->GetOptIndex(NodesVariables::NodeValueInfo(f_node_id-1, kPos, X_))) = +1.0/t_swing_avg_; // previous node
//////                            row_reset++;
////
////                            //pos und vel in y rtg bei "X constraint" muss 0 sein
////                            jac.coeffRef(row_reset++, ee_motion_->GetOptIndex(NodesVariables::NodeValueInfo(f_node_id,   kPos, X_))) =  0;
////                            jac.coeffRef(row_reset++, ee_motion_->GetOptIndex(NodesVariables::NodeValueInfo(f_node_id,   kVel, X_))) =  0;
////
////
//////              auto curr = nodes.at(f_node_id);
//////
//////                      //TODO no velocity in Y direction during drive!
//////              //        Vector2d prev;
//////              //        Vector2d next;
//////
//////                      Vector2d prev = nodes.at(f_node_id-1).p().topRows<2>();	//1 for just x-direction
//////              //        prev.row(1).setZero();
//////
//////                      Vector2d next = nodes.at(f_node_id+1).p().topRows<2>();
//////              //        next.row(1).setZero();
//////
//////              //        Vector2d ybound = nodes.at(f_node_id).p().topRows<2>();
//////
//////                      Vector2d distance_xy    = next - prev;
//////                      Vector2d xy_center      = prev + 0.5*distance_xy;
//////                      Vector2d des_vel_center = distance_xy/t_swing_avg_;
////////
//////                      jac.coeffRef(row_reset++, idx1) = curr.p()(X) - xy_center(X);
//////                      jac.coeffRef(row_reset++, idx1) = curr.v()(X) - des_vel_center(X);
//////                      jac.coeffRef(row_reset++, idx1) = 0;
//////                      jac.coeffRef(row_reset++, idx1) = 0;
////
////                            Vector3d dny  = terrain_->GetDerivativeOfNormalizedBasisWrt(HeightMap::Normal, Y_, p.x(), p.y());
////                                          Vector3d dt1y = terrain_->GetDerivativeOfNormalizedBasisWrt(HeightMap::Tangent1, Y_, p.x(), p.y());
////                                          Vector3d dt2y = terrain_->GetDerivativeOfNormalizedBasisWrt(HeightMap::Tangent2, Y_, p.x(), p.y());
////
////                                          //TODO what should ee_node_id be? (not only ID at Start of phase..?
////                                          int idx2 = ee_motion_->GetOptIndex(NodesVariables::NodeValueInfo(ee_node_id, kPos, Y_));
////                                          int row_resety=row;
////
////                                          // unilateral force
////                                          jac.coeffRef(row_resety++, idx2) = f.transpose()*dny;
////
////                                          // friction force tangent 1 derivative
////                                          jac.coeffRef(row_resety++, idx2) = f.transpose()*(dt1y-mu_*dny);
////                                          jac.coeffRef(row_resety++, idx2) = f.transpose()*(dt1y+mu_*dny);
////
////                                          // friction force tangent 2 derivative
////                                          jac.coeffRef(row_resety++, idx2) = f.transpose()*(dt2y-mu_*dny);
////                                          jac.coeffRef(row_resety++, idx2) = f.transpose()*(dt2y+mu_*dny);
////
////                                          jac.coeffRef(row_resety++, ee_motion_->GetOptIndex(NodesVariables::NodeValueInfo(f_node_id,   kPos, Y_))) =  0;
////                                          jac.coeffRef(row_resety++, ee_motion_->GetOptIndex(NodesVariables::NodeValueInfo(f_node_id,   kVel, Y_))) =  0;
//////                                          //new
//////                                          // position constraint
////                                                        jac.coeffRef(row_resety++, ee_motion_->GetOptIndex(NodesVariables::NodeValueInfo(f_node_id,   kPos, Y_))) =  curr.p()(Y) - prev(Y);  // current node
//////                                                        jac.coeffRef(row_resety, ee_motion_->GetOptIndex(NodesVariables::NodeValueInfo(f_node_id+1, kPos, Y_))) = 1.0;  // next node
//////                                                        jac.coeffRef(row_resety, ee_motion_->GetOptIndex(NodesVariables::NodeValueInfo(f_node_id-1, kPos, Y_))) = 1.0;  // previous node
//////                                                        row_resety++;
//////
//////                                          // velocity constraint
////                                                        jac.coeffRef(row_resety++, ee_motion_->GetOptIndex(NodesVariables::NodeValueInfo(f_node_id,   kVel, Y_))) =  curr.v()(Y);              // current node
//////                                                        jac.coeffRef(row_resety, ee_motion_->GetOptIndex(NodesVariables::NodeValueInfo(f_node_id+1, kPos, Y_))) = 1.0/t_swing_avg_; // next node
//////                                                        jac.coeffRef(row_resety, ee_motion_->GetOptIndex(NodesVariables::NodeValueInfo(f_node_id-1, kPos, Y_))) = 1.0/t_swing_avg_; // previous node
//////                                                        row_resety++;
////
////
////
//////                                          jac.coeffRef(row_resety++, idx2) = 0;
//////                                          jac.coeffRef(row_resety++, idx2) = 0;
//////                                          jac.coeffRef(row_resety++, idx2) = curr.p()(Y) - prev(Y);
//////                                          jac.coeffRef(row_resety++, idx2) = curr.v()(Y);
//      row += n_constraints_per_node_;
//    }
////    cout << "n_rows_motion_Jacobian: " << row;
////    n_rows_motion_Jacobian = row;
//  }
//  }
}

} /* namespace towr */
