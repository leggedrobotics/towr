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

#include <towr/constraints/terrain_constraint.h>


namespace towr {


TerrainConstraint::TerrainConstraint (const HeightMap::Ptr& terrain,
                                      std::string ee_motion)
    :ConstraintSet(kSpecifyLater, "terrain-" + ee_motion)
{
  ee_motion_id_ = ee_motion;
  terrain_ = terrain;
}

void
TerrainConstraint::InitVariableDependedQuantities (const VariablesPtr& x)
{
  ee_motion_ = x->GetComponent<NodesVariablesPhaseBased>(ee_motion_id_);

  // skip first node, b/c already constrained by initial stance
  for (int id=1; id<ee_motion_->GetNodes().size(); ++id)
    node_ids_.push_back(id);

  int constraint_count = 2*node_ids_.size();
  SetRows(constraint_count);
}

Eigen::VectorXd
TerrainConstraint::GetValues () const
{
  VectorXd g(GetRows());

  auto nodes = ee_motion_->GetNodes();
  int row = 0;
  for (int id : node_ids_) {
    Vector3d p = nodes.at(id).p();
    Vector3d v = nodes.at(id).v();
    g(row++) = p.z() - terrain_->GetHeight(p.x(), p.y());
    Vector3d terrain_normal = terrain_->GetNormalizedBasis(HeightMap::Direction::Normal, p.x(), p.y());
    g(row++) = terrain_normal.x() * v.x() + terrain_normal.y() * v.y() + terrain_normal.z() * v.z();
  }

  return g;
}

TerrainConstraint::VecBound
TerrainConstraint::GetBounds () const
{
  VecBound bounds(GetRows());
  double max_distance_above_terrain = 1e20; // [m]
  double max_velocity_normal_to_terrain = 1e20; // [m/s]

  int row = 0;
  for (int id : node_ids_) {
    if (ee_motion_->IsContactNode(id)) {
      bounds.at(row++) = ifopt::BoundZero;
      bounds.at(row++) = ifopt::BoundZero;
    } else {
      bounds.at(row++) = ifopt::Bounds(0.0, max_distance_above_terrain);
      bounds.at(row++) = ifopt::Bounds(-max_velocity_normal_to_terrain, max_velocity_normal_to_terrain);
    }
  }

  return bounds;
}

void
TerrainConstraint::FillJacobianBlock (std::string var_set, Jacobian& jac) const
{
  if (var_set == ee_motion_->GetName()) {
    auto nodes = ee_motion_->GetNodes();
    int row = 0;
    for (int id : node_ids_) {
      int idx = ee_motion_->GetOptIndex(NodesVariables::NodeValueInfo(id, kPos, Z));
      jac.coeffRef(row, idx) = 1.0;

      Vector3d p = nodes.at(id).p();
      for (auto dim : {X,Y}) {
        int idx = ee_motion_->GetOptIndex(NodesVariables::NodeValueInfo(id, kPos, dim));
        jac.coeffRef(row, idx) = -terrain_->GetDerivativeOfHeightWrt(To2D(dim), p.x(), p.y());
      }
      row++;

      Vector3d v = nodes.at(id).v();
      Vector3d terrain_normal = terrain_->GetNormalizedBasis(HeightMap::Direction::Normal, p.x(), p.y());
      int idvz = ee_motion_->GetOptIndex(NodesVariables::NodeValueInfo(id, kVel, X));
      jac.coeffRef(row, idvz) =  terrain_normal.x();
      int idvy = ee_motion_->GetOptIndex(NodesVariables::NodeValueInfo(id, kVel, Y));
      jac.coeffRef(row, idvy) =  terrain_normal.y();
      int idvx = ee_motion_->GetOptIndex(NodesVariables::NodeValueInfo(id, kVel, Z));
      jac.coeffRef(row, idvx) =  terrain_normal.z();
      for (auto dim : {X,Y}) {
        int idx = ee_motion_->GetOptIndex(NodesVariables::NodeValueInfo(id, kPos, dim));
        Vector3d terrain_normal_derivative = terrain_->GetDerivativeOfNormalizedBasisWrt(HeightMap::Direction::Normal, To2D(dim), p.x(), p.y());
        jac.coeffRef(row, idx) = terrain_normal_derivative.x() * v.x() + terrain_normal_derivative.y() * v.y() + terrain_normal_derivative.z() * v.z();
      }
      row++;
    }
  }
}

} /* namespace towr */
