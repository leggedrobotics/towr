///******************************************************************************
//Copyright (c) 2018, Alexander W. Winkler. All rights reserved.
//
//Redistribution and use in source and binary forms, with or without
//modification, are permitted provided that the following conditions are met:
//
//* Redistributions of source code must retain the above copyright notice, this
//  list of conditions and the following disclaimer.
//
//* Redistributions in binary form must reproduce the above copyright notice,
//  this list of conditions and the following disclaimer in the documentation
//  and/or other materials provided with the distribution.
//
//* Neither the name of the copyright holder nor the names of its
//  contributors may be used to endorse or promote products derived from
//  this software without specific prior written permission.
//
//THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
//AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
//IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
//FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
//DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
//SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
//CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
//OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
//OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//******************************************************************************/
//
//#include <towr/constraints/swing_constraint.h>
//#include <towr/variables/cartesian_dimensions.h>
//#include <towr/variables/variable_names.h>
//#include <iostream>
//using namespace std;
//
//
//namespace towr {
//
//SwingConstraint::SwingConstraint (EE ee,
//									const SplineHolder& spline_holder)
//    :ConstraintSet(kSpecifyLater, "swing-" + id::EEMotionNodes(ee))
//{
////  ee_motion_id_ = ee_motion;
//  ee_ = ee;
//  base_angular_ = EulerConverter(spline_holder.base_angular_);
////  cout << "test if swing constraint active " << endl;
//}
//
//void
//towr::SwingConstraint::InitVariableDependedQuantities (const VariablesPtr& x)
//{
//  ee_force_  = x->GetComponent<NodesVariablesPhaseBased>(id::EEForceNodes(ee_));
//  ee_motion_ = x->GetComponent<NodesVariablesPhaseBased>(id::EEMotionNodes(ee_));
//
////  pure_swing_node_ids_ = ee_motion_->GetIndicesNodesWOFirstAndLast();
//  pure_swing_node_ids_ = ee_motion_->GetIndicesOfAllNodes();
//
//
////  pure_swing_node_ids_ = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18};
////  pure_swing_node_ids_ = {1,2,3};
//
//  // constrain xy position and velocity of every swing node
////  int constraint_count =  pure_swing_node_ids_.size()*Node::n_derivatives*1;
////  int constraint_count =  pure_swing_node_ids_.size()*k2D;
//
//  int constraint_count =  0; //noch anpassen! (0,3,4 phase keine constraints)
//
////  if (phase == 1 or phase == 2 or phase == 5){
////	  constraint_count = pure_swing_node_ids_.size()*2;
////  }
//
//
//  //only constrain drive nodes!:
////  int constraint_count =  16*Node::n_derivatives*k2D;
//
////  int constraint_count =  1;
//
//  SetRows(constraint_count);
////  SetRows(0);
////  cout << "number of nodes: " << pure_swing_node_ids_.size() << endl;
////  cout << "swing constraints: " << constraint_count << endl;
//}
//
//Eigen::VectorXd
//SwingConstraint::GetValues () const
//{
//  VectorXd g(GetRows());
//
//  int row = 0;
//  auto nodes = ee_motion_->GetNodes();
//  for (int node_id : pure_swing_node_ids_) {
//    // assumes two splines per swingphase and starting and ending in stance
//    auto curr = nodes.at(node_id);
//    int phase  = ee_motion_->GetPhase(node_id, ee_);
//
//    int ee_id = ee_;
//    //    	cout << "ee_id: " << ee_id << endl;
//
////    //get times of different phases and for different ee
//////        	double time_3nodes = 2.5/params_.ee_polynomials_per_swing_phase_*3; //2.5 is total time!
////        	double t = 0.0; //get time from the current node!
////        	double t_per_node = 2.5/pure_swing_node_ids_.size();
////
////        	//TODO anpassen fuer allgemeinen Fall!!
////        	if (phase == 0 or phase == 3){
////        		t = 0.0;
////        	}
////        	if (phase == 1){
////        		t = (node_id+1)*t_per_node;
////        	}
////        	if (phase == 2 or phase == 5){
////        		t = 2.0;
////        	}
////        	//phase 4 doesnt matter!
////
////
////    //    	Spline::GetLocalTime()
////
////        	  EulerConverter::MatrixSXd b_C_w = base_angular_.GetRotationMatrixBaseToWorld(t);
////        	  EulerConverter::MatrixSXd w_C_b = base_angular_.GetRotationMatrixBaseToWorld(t).transpose();
////
////        	  Vector3d v_wrt_b = w_C_b*nodes.at(node_id).at(kVel);
////        	  Vector3d v_b_y = {0, v_wrt_b(1), 0};
////        	  Vector3d v_b_x = {v_wrt_b(0), 0, 0};
////        	  Vector3d v_x = b_C_w*v_b_x;
////        	  Vector3d v_y = b_C_w*v_b_y;
////
//
////        	  if (phase == 1 or phase == 2 or phase == 5){
////        		  g(row++) = v_y(0);	//derivatives not easy!
////        		  g(row++) = v_y(1);
////        	  }
//
//
////        	  if (phase == 0 and (ee_ == 2 or ee_ == 3)){
//        	  //    		g(row++) = nodes.at(f_node_id).v().y();
//        	  //    		g(row++) = 0;
//        	  ////    		g(row++) = nodes.at(f_node_id).v().x();
//        	  //    	}
//        	  //    	else if (f_node_id < 8){
//        	  //    		g(row++) = nodes.at(f_node_id).v().y();
//        	  //    		g(row++) = 0;
//        	  ////    		g(row++) = nodes.at(f_node_id).v().x();
//        	  ////    		g(row++) = 0;
//        	  //    	}
//        	  //    	else {
//        	  //    		g(row++) = 0;
//        	  //    		g(row++) = 0;
//        	  ////    		g(row++) = v_y(0);
//        	  ////    		g(row++) = v_y(1);
//        	  ////    		g(row++) = v_x(0);		//caution, robot cant turn drift more than 90degrees,
//        	  //    								//otherwise it will have to drive in negative x direction
//        	  //    	}
//
//        	// add y-motion constraints
////        	  int n_nodes_drive = (params_.ee_polynomials_per_swing_phase_+1)/2.5; //2.5 = total time
//
//        	  //for first drive phase of hind wheels, world y-vel is 0
//
////        	  if (phase == 0 and (ee_id == 2 or ee_id == 3)){
//////        		  cout << "the ee_id condition is met." << endl;
////        	      		g(row++) = nodes.at(node_id).v().y();
//////        	      		g(row++) = 0;
////        	  //    		g(row++) = nodes.at(f_node_id).v().x();
////        	      	}
////        	  //for first 8 nodes of phase 0 of front wheels, world y-vel is 0
////
////        	  else if (phase == 0 and node_id < n_nodes_drive and (ee_id == 0 or ee_id == 1)){
////        	      		g(row++) = nodes.at(node_id).v().y();
//////        	      		g(row++) = 0;
////        	  //    		g(row++) = nodes.at(f_node_id).v().x();
////        	  //    		g(row++) = 0;
////        	      	}
////
////        	  else if (phase == 1 or phase == 2 or (phase == 0 and node_id >= n_nodes_drive and (ee_id == 0 or ee_id == 1))){
////        	      		g(row++) = 0;
//////        	      		g(row++) = 0;
////        	  //    		g(row++) = v_y(0);
////        	  //    		g(row++) = v_y(1);
////        	  //    		g(row++) = v_x(0);		//caution, robot cant turn drift more than 90degrees,
////        	      								//otherwise it will have to drive in negative x direction
////        	      	}
//
//////    Vector2d prev = nodes.at(node_id-1).p().topRows<k2D>();
//////    Vector2d next = nodes.at(node_id+1).p().topRows<k2D>();
//////
//////        	Vector2d prev = nodes.at(0).p().topRows<k2D>();
//////        	Vector2d next = nodes.at(0).p().topRows<k2D>();
//////        	Vector2d distance_xy;
//////        	Vector2d xy_center;
//////        	Vector2d des_vel_center;
////
//////        	    if (node_id > 0 and node_id < (pure_swing_node_ids_.size()-1)){
////        	Vector2d prev = nodes.at(node_id-1).p();
////        	Vector2d next = nodes.at(node_id+1).p();
////        	Vector2d distance_xy    = next - prev;
////        	Vector2d xy_center      = prev + 0.5*distance_xy;
////        	Vector2d des_vel_center = distance_xy/time_3nodes;
////
//////        	    }
////
//////        	    else if (node_id == 0) {
//////        	    	next = nodes.at(node_id+1).p();
//////
//////        	    	distance_xy    = next - prev;
//////        	    	xy_center      = prev;
//////        	    	double time_2nodes = time_3nodes/3*2;
//////        	    	des_vel_center = distance_xy/time_2nodes;
//////        	    }
//////        	    else if (node_id == (pure_swing_node_ids_.size()-1)){
//////        	    	prev = nodes.at(node_id-1).p();
//////        	    	next = nodes.at(node_id).p();
//////        	    	distance_xy    = next - prev;
//////        	    	xy_center      = next;
//////        	    	double time_2nodes = time_3nodes/3*2;
//////        	    	des_vel_center = distance_xy/time_2nodes;
//////        	    }
////
////
//////    if (phase == 0){
//////    	des_vel_center = distance_xy/0.45; // linear interpolation not accurate
//////    }
//////    else {
//////    	des_vel_center = distance_xy/0.25;
//////    }
////
//////    g(row++) = curr.p()(X) - xy_center(X);
//////    g(row++) = curr.v()(X) - des_vel_center(X);
//////    if (phase == 0){
////
////    for (auto dim : {X,Y}) {
//////    	if (node_id > 0 and node_id < (pure_swing_node_ids_.size()-1)){
////      g(row++) = curr.p()(dim) - xy_center(dim);
////      g(row++) = curr.v()(dim) - des_vel_center(dim);
//////    	}
//////    	else {
//////    		g(row++) = 0;
////////    		g(row++) = 0;
//////    	}
////    }
//////    g(row++) = curr.v()(X);
////
//////    Vector2d distance_xy    = next - prev;
//////    Vector2d xy_center      = prev + distance_xy;
//////    Vector2d des_vel_center = distance_xy/(2.8/3); // linear interpolation not accurate
//////    g(row++) = curr.v()(X) - des_vel_center(X);
//
//  }
//
//  return g;
//}
//
//SwingConstraint::VecBound
//SwingConstraint::GetBounds () const
//{
//	VecBound bounds;
////  return VecBound(GetRows(), ifopt::BoundGreaterZero);
////	return VecBound(GetRows(), ifopt::BoundZero);
//	for (int node_id : pure_swing_node_ids_) {
//		bounds.push_back(ifopt::BoundZero); //
////		bounds.push_back(ifopt::BoundZero); //
////		bounds.push_back(ifopt::BoundZero); //
////		bounds.push_back(ifopt::BoundZero); //
//	}
//
//	return bounds;
//}
//
//void
//SwingConstraint::FillJacobianBlock (std::string var_set,
//                                    Jacobian& jac) const
//{
//  if (var_set == ee_motion_->GetName()) {
//    int row = 0;
//    for (int node_id : pure_swing_node_ids_) {
//    	int phase  = ee_motion_->GetPhase(node_id,ee_);
//    	int ee_id = ee_;
//    	    //    	cout << "ee_id: " << ee_id << endl;
//
//    	        	//get times of different phases and for different ee
////    	        	double time_3nodes = 0.0;
////    	        	if (phase == 0){
////    	        		if (node_id > 7 and node_id < 15 and (ee_id == 0 or ee_id == 1)){ //and add condition that this is only valid for front feet!
////    	        			time_3nodes = 0.375;
////    	        		}
////    	    //    		t = (node_id+1)*0.0625;
////    	        		else if (node_id > 15 and (ee_id == 0 or ee_id == 1)){
////    	        			time_3nodes = 0.375;
////    	        		}
////    	        		else if (ee_id == 2 or ee_id == 3){
////    	        			time_3nodes = 0.15;
////    	        		}
////    	        	}
////    	        	else if (phase == 2){
////    	        		time_3nodes = 0.075;
////    	        	}
////    	        	else if (phase == 1){
////    	        		time_3nodes = 0.15;
////    	        	}
////    	        	int n_nodes_drive = params_.ee_polynomials_per_swing_phase_/2.5; //2.5 = total time
////
////    	        	if (phase == 0 and (ee_id == 2 or ee_id == 3)){
////    	        		jac.coeffRef(row++, ee_motion_->GetOptIndex(NodesVariables::NodeValueInfo(node_id,   kVel, Y_))) =  1.0;
////
////    	        	        	      	}
////    	        	        	  //for first 8 nodes of phase 0 of front wheels, world y-vel is 0
////    	        	else if (phase == 0 and node_id < n_nodes_drive and (ee_id == 0 or ee_id == 1)){
////    	        	        	      		jac.coeffRef(row++, ee_motion_->GetOptIndex(NodesVariables::NodeValueInfo(node_id,   kVel, Y_))) =  1.0;
////    	        	        	      	}
////    	        	else if (phase == 1 or phase == 2 or (phase == 0 and node_id >= n_nodes_drive and (ee_id == 0 or ee_id == 1))){
////    	        	        	      		jac.coeffRef(row++, ee_motion_->GetOptIndex(NodesVariables::NodeValueInfo(node_id,   kVel, Y_))) =  0.0;
////    	        	        	      	}
//
//
////    	for (auto dim : {X,Y}) {
////        // position constraint
//////    	if (node_id > 0 and node_id < (pure_swing_node_ids_.size()-1)){
////        jac.coeffRef(row, ee_motion_->GetOptIndex(NodesVariables::NodeValueInfo(node_id,   kPos, dim))) =  1.0;  // current node
//////        jac.coeffRef(row, ee_motion_->GetOptIndex(NodesVariables::NodeValueInfo(node_id+1, kPos, dim))) = -0.5;  // next node
//////        jac.coeffRef(row, ee_motion_->GetOptIndex(NodesVariables::NodeValueInfo(node_id-1, kPos, dim))) = -0.5;  // previous node
////        row++;
////
////        // velocity constraint
////        jac.coeffRef(row, ee_motion_->GetOptIndex(NodesVariables::NodeValueInfo(node_id,   kVel, dim))) =  1.0;              // current node
//////        jac.coeffRef(row, ee_motion_->GetOptIndex(NodesVariables::NodeValueInfo(node_id+1, kPos, dim))) = -1.0/time_3nodes; // next node
//////        jac.coeffRef(row, ee_motion_->GetOptIndex(NodesVariables::NodeValueInfo(node_id-1, kPos, dim))) = +1.0/time_3nodes; // previous node
////        row++;
//////    	}
//////    	else {
//////    		jac.coeffRef(row, ee_motion_->GetOptIndex(NodesVariables::NodeValueInfo(node_id,   kPos, dim))) =  0.0;  // current node
//////    		        jac.coeffRef(row, ee_motion_->GetOptIndex(NodesVariables::NodeValueInfo(node_id+1, kPos, dim))) = 0.0;  // next node
//////    		        jac.coeffRef(row, ee_motion_->GetOptIndex(NodesVariables::NodeValueInfo(node_id-1, kPos, dim))) = 0.0;  // previous node
//////    		        row++;
//////    		        jac.coeffRef(row, ee_motion_->GetOptIndex(NodesVariables::NodeValueInfo(node_id,   kVel, dim))) =  0.0;              // current node
//////    		                jac.coeffRef(row, ee_motion_->GetOptIndex(NodesVariables::NodeValueInfo(node_id+1, kPos, dim))) = 0.0; // next node
//////    		                jac.coeffRef(row, ee_motion_->GetOptIndex(NodesVariables::NodeValueInfo(node_id-1, kPos, dim))) = 0.0; // previous node
//////    		                row++;
//////
//////    	}
////      }
////    	jac.coeffRef(row, ee_motion_->GetOptIndex(NodesVariables::NodeValueInfo(node_id,   kVel, X))) =  1.0;
////    	row++;
//    }
//  }
//}
//
//} /* namespace towr */


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

#include <towr/constraints/swing_constraint.h>
#include <towr/variables/cartesian_dimensions.h>
#include <towr/variables/variable_names.h>
#include <iostream>

//#include <towr_ros/towr_user_interface.h>
using namespace std;


namespace towr {

SwingConstraint::SwingConstraint (EE ee,
									const SplineHolder& spline_holder)
    :ConstraintSet(kSpecifyLater, "swing-" + id::EEMotionNodes(ee))
{
  ee_ = ee;
  base_angular_ = spline_holder.base_angular_;
}

void
towr::SwingConstraint::InitVariableDependedQuantities (const VariablesPtr& x)
{

}

Eigen::VectorXd
SwingConstraint::GetValues () const
{
  VectorXd g(GetRows());

  return g;
}

SwingConstraint::VecBound
SwingConstraint::GetBounds () const
{
	VecBound bounds;
	return bounds;
}

void
SwingConstraint::FillJacobianBlock (std::string var_set,
                                    Jacobian& jac) const
{

}

} /* namespace towr */
