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

#include <towr/variables/nodes_variables.h>
#include <iostream>
#include <towr/variables/node_spline.h>
#include <towr/variables/euler_converter.h>

namespace towr {


NodesVariables::NodesVariables (const std::string& name)
    : VariableSet(kSpecifyLater, name)
{
}

int
NodesVariables::GetOptIndex(const NodeValueInfo& nvi_des) const
{
  // could also cache this as map for more efficiency, but adding complexity
  for (int idx=0; idx<GetRows(); ++idx)
    for ( NodeValueInfo nvi : GetNodeValuesInfo(idx))
      if ( nvi == nvi_des )
        return idx;

  return NodeValueNotOptimized; // index representing these quantities doesn't exist
}

Eigen::VectorXd
NodesVariables::GetValues () const
{
  VectorXd x(GetRows());

  for (int idx=0; idx<x.rows(); ++idx)
    for (auto nvi : GetNodeValuesInfo(idx))
      x(idx) = nodes_.at(nvi.id_).at(nvi.deriv_)(nvi.dim_);

  return x;
}

void
NodesVariables::SetVariables (const VectorXd& x)
{
  for (int idx=0; idx<x.rows(); ++idx)
    for (auto nvi : GetNodeValuesInfo(idx))
      nodes_.at(nvi.id_).at(nvi.deriv_)(nvi.dim_) = x(idx);

  UpdateObservers();
}

void
NodesVariables::UpdateObservers() const
{
  for (auto& o : observers_)
    o->UpdateNodes();
}

void
NodesVariables::AddObserver(ObserverPtr const o)
{
   observers_.push_back(o);
}

int
NodesVariables::GetNodeId (int poly_id, Side side)
{
  return poly_id + side;
}

const std::vector<Node>
NodesVariables::GetBoundaryNodes(int poly_id) const
{
  std::vector<Node> nodes;
  nodes.push_back(nodes_.at(GetNodeId(poly_id, Side::Start)));
  nodes.push_back(nodes_.at(GetNodeId(poly_id, Side::End)));
  return nodes;
}

int
NodesVariables::GetDim() const
{
  return n_dim_;
}

int
NodesVariables::GetPolynomialCount() const
{
  return nodes_.size() - 1;
}

NodesVariables::VecBound
NodesVariables::GetBounds () const
{
  return bounds_;
}

const std::vector<Node>
NodesVariables::GetNodes() const
{
  return nodes_;
}

void
NodesVariables::AdvancedInititialisationEE(const VectorXd& initial_val,
                                          const VectorXd& final_val,
                                          double t_total,
                                          std::vector<double> timings,
                                          double des_w,
                                          double des_vx,
                                          double des_vy,
                                          double z_offset,
                                          const VectorXd& offset_full,
                                          HeightMap::Ptr  terrain,
                                          std::vector<int> poly_per_phase,
                                          double angle_init,
                                          bool incontact_start)
{
  double theta0 = atan2(des_vy, des_vx);
  double r = sqrt(des_vx * des_vx + des_vy * des_vy) / des_w;
  double t_current3 = 0.0;
  int id_prev_ = 0;

  int lastnodeadded = -10;
  double goalz_terrain_last= terrain->GetHeight(initial_val.x(),initial_val.y()) - z_offset;

  int phase = 0;
  int polycumulative= poly_per_phase[0];
  bool contact = incontact_start;

  for (int idx=0; idx<GetRows(); ++idx) {
    for (auto nvi : GetNodeValuesInfo(idx)) {
      if ( id_prev_!=nvi.id_){
        if(nvi.id_> polycumulative){
          phase+=1;
          contact=!contact;
          polycumulative+= poly_per_phase[phase];
        }
        t_current3 += (1.0/poly_per_phase[phase])*timings[phase];
        id_prev_ = nvi.id_;
      }

      double thetaT = theta0 + des_w * t_current3;
      double yaw=angle_init+ des_w * t_current3;
      double goalx;
      double goaly;
      double goalz;
      double goalvx;
      double goalvy;
      double goalvz;

      Eigen::Vector3d euler(0.0, 0.0, yaw);
      Eigen::Vector3d deswvector(0.0, 0.0, des_w);
      Eigen::Matrix3d w_R_b = EulerConverter::GetRotationMatrixBaseToWorld(euler);

      Eigen::Vector3d ofsetre = offset_full;
      Eigen::Vector3d desv(des_vx, des_vy, 0.0);

      VectorXd offset_rotated = w_R_b*offset_full;
      VectorXd desv_rotated = w_R_b*(desv+deswvector.cross(ofsetre));

      if(des_w==0){
        goalx = initial_val.x() + (t_current3*(final_val.x()-initial_val.x()))/t_total;
        goaly = initial_val.y() + (t_current3*(final_val.y()-initial_val.y()))/t_total;
        goalvx = desv_rotated.x();
        goalvy = desv_rotated.y();
      }else{
        goalx = offset_rotated.x() -r * sin(theta0) + r * sin(thetaT);
        goaly = offset_rotated.y() + r * cos(theta0) - r * cos(thetaT);
        goalvx = desv_rotated.x();
        goalvy = desv_rotated.y();
      }

      double goalz_terrain_ = terrain->GetHeight(goalx,goaly) - z_offset;

      if(contact){
        goalz = goalz_terrain_;
      }else{
        goalz = goalz_terrain_+0.1;
      }
      double dzdx = terrain->GetDerivativeOfHeightWrt(X_,goalx,goaly);
      double dzdy = terrain->GetDerivativeOfHeightWrt(Y_,goalx,goaly);
      double dzdt = dzdx * goalvx + dzdy * goalvy;
      goalvz = dzdt;//terrain gradient

      if (nvi.deriv_ == kPos) {
        Eigen::Vector3d pos;
        pos.x()= goalx;
        pos.y()= goaly;
        pos.z()= goalz;
        nodes_.at(nvi.id_).at(kPos)(nvi.dim_) = pos(nvi.dim_);
      }

      if (nvi.deriv_ == kVel) {
        Eigen::Vector3d vel;
        vel.x()= goalvx;
        vel.y()= goalvy;
        vel.z()= goalvz;
        nodes_.at(nvi.id_).at(kVel)(nvi.dim_) = vel(nvi.dim_);
      }
    }
  }
}

void
NodesVariables::AdvancedInititialisationBase(const VectorXd& initial_val,
                                          const VectorXd& final_val,
                                          double t_total,
                                          double constant_timings,
                                          double des_w,
                                          double des_vx,
                                          double des_vy,
                                          double z_offset,
                                          HeightMap::Ptr  terrain,
                                          double angle_init)
{
  double theta0 = atan2(des_vy, des_vx);
  double r = sqrt(des_vx * des_vx + des_vy * des_vy) / des_w;
  double t_current3 = 0.0;
  int id_prev_ = 0;

  for (int idx=0; idx<GetRows(); ++idx) {
    for (auto nvi : GetNodeValuesInfo(idx)) {
      if ( id_prev_!=nvi.id_){
        t_current3 += constant_timings;
        id_prev_ = nvi.id_;
      }

      double thetaT = theta0 + des_w * t_current3;
      double yaw=angle_init+ des_w * t_current3;
      double goalx;
      double goaly;
      double goalvx;
      double goalvy;

      Eigen::Vector3d euler(0.0, 0.0, yaw);
      Eigen::Vector3d deswvector(0.0, 0.0, des_w);
      Eigen::Matrix3d w_R_b = EulerConverter::GetRotationMatrixBaseToWorld(euler);

      Eigen::Vector3d desv(des_vx, des_vy, 0.0);

      VectorXd desv_rotated = w_R_b*desv;

      if(des_w==0){
        goalx = initial_val.x() + (t_current3*(final_val.x()-initial_val.x()))/t_total;
        goaly = initial_val.y() + (t_current3*(final_val.y()-initial_val.y()))/t_total;
        goalvx = desv_rotated.x();
        goalvy = desv_rotated.y();
      }else{
        goalx = -r * sin(theta0) + r * sin(thetaT);
        goaly =  r * cos(theta0) - r * cos(thetaT);
        goalvx = desv_rotated.x();
        goalvy = desv_rotated.y();
      }

      double goalz = initial_val.z() + (t_current3*(final_val.z()-initial_val.z()))/t_total;
      double goalvz =(final_val.z()-initial_val.z())/t_total;


      if (nvi.deriv_ == kPos) {
        Eigen::Vector3d pos;
        pos.x()= goalx;
        pos.y()= goaly;
        pos.z()= goalz;
        nodes_.at(nvi.id_).at(kPos)(nvi.dim_) = pos(nvi.dim_);
      }

      if (nvi.deriv_ == kVel) {
        Eigen::Vector3d vel;
        vel.x()= goalvx;
        vel.y()= goalvy;
        vel.z()= goalvz;
        nodes_.at(nvi.id_).at(kVel)(nvi.dim_) = vel(nvi.dim_);
      }
    }
  }
}

void
NodesVariables::SetByLinearInterpolation(const VectorXd& initial_val,
                                         const VectorXd& final_val,
                                         double t_total)
{
  // only set those that are part of optimization variables,
  // do not overwrite phase-based parameterization
  VectorXd dp = final_val-initial_val;
  VectorXd average_velocity = dp / t_total;
  int num_nodes = nodes_.size();

  for (int idx=0; idx<GetRows(); ++idx) {
    for (auto nvi : GetNodeValuesInfo(idx)) {

      if (nvi.deriv_ == kPos) {
        VectorXd pos = initial_val + nvi.id_/static_cast<double>(num_nodes-1)*dp;
        nodes_.at(nvi.id_).at(kPos)(nvi.dim_) = pos(nvi.dim_);
      }

      if (nvi.deriv_ == kVel) {
        nodes_.at(nvi.id_).at(kVel)(nvi.dim_) = average_velocity(nvi.dim_);
      }
    }
  }
}

void
NodesVariables::AddBounds(int node_id, Dx deriv,
                 const std::vector<int>& dimensions,
                 const VectorXd& val)
{
  for (auto dim : dimensions)
    AddBound(NodeValueInfo(node_id, deriv, dim), val(dim));
}

void
NodesVariables::AddBound (const NodeValueInfo& nvi_des, double val)
{
  for (int idx=0; idx<GetRows(); ++idx)
    for (auto nvi : GetNodeValuesInfo(idx))
      if (nvi == nvi_des)
        bounds_.at(idx) = ifopt::Bounds(val, val);
}

void
NodesVariables::AddStartBound (Dx d, const std::vector<int>& dimensions, const VectorXd& val)
{
  AddBounds(0, d, dimensions, val);
}

void
NodesVariables::AddFinalBound (Dx deriv, const std::vector<int>& dimensions,
                      const VectorXd& val)
{
  AddBounds(nodes_.size()-1, deriv, dimensions, val);
}

NodesVariables::NodeValueInfo::NodeValueInfo(int node_id, Dx deriv, int node_dim)
{
  id_    = node_id;
  deriv_ = deriv;
  dim_   = node_dim;
}

int
NodesVariables::NodeValueInfo::operator==(const NodeValueInfo& right) const
{
  return (id_    == right.id_)
      && (deriv_ == right.deriv_)
      && (dim_   == right.dim_);
};

} /* namespace towr */
