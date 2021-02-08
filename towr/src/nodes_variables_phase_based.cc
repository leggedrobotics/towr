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

#include <towr/variables/nodes_variables_phase_based.h>
#include <towr/variables/cartesian_dimensions.h>

#include <iostream>

namespace towr {


std::vector<NodesVariablesPhaseBased::PolyInfo>
BuildPolyInfos (int phase_count, bool first_phase_constant,
        		std::vector<int> n_polys_per_phase)
{
  using PolyInfo = NodesVariablesPhaseBased::PolyInfo;
  std::vector<PolyInfo> polynomial_info;

  bool phase_constant = first_phase_constant;

  for (int i=0; i<phase_count; ++i) {
    int n_polys = n_polys_per_phase.at(i);
	if (phase_constant) {
	  for (int j=0; j<n_polys; ++j)
	    polynomial_info.push_back(PolyInfo(i,j,n_polys, true, false));
    }
    else {
      for (int j=0; j<n_polys; ++j)
        polynomial_info.push_back(PolyInfo(i,j,n_polys, false, false));
    }

    phase_constant = !phase_constant; // constant and non-constant phase alternate
  }

  return polynomial_info;
}

NodesVariablesPhaseBased::NodesVariablesPhaseBased (int phase_count,
                                                    bool first_phase_constant,
                                                    const std::string& name,
													std::vector<int> n_polys_per_phase)
    : NodesVariables(name)
{
  polynomial_info_ = BuildPolyInfos(phase_count, first_phase_constant, n_polys_per_phase);

  n_dim_ = k3D;
  int n_nodes = polynomial_info_.size()+1;
  nodes_  = std::vector<Node>(n_nodes, Node(n_dim_));
}

NodesVariablesPhaseBased::VecDurations
NodesVariablesPhaseBased::ConvertPhaseToPolyDurations(const VecDurations& phase_durations) const
{
  VecDurations poly_durations;

  for (int i=0; i<GetPolynomialCount(); ++i) {
    auto info = polynomial_info_.at(i);
    poly_durations.push_back(phase_durations.at(info.phase_)/info.n_polys_in_phase_);
  }

  return poly_durations;
}

double
NodesVariablesPhaseBased::GetDerivativeOfPolyDurationWrtPhaseDuration (int poly_id) const
{
  int n_polys_in_phase = polynomial_info_.at(poly_id).n_polys_in_phase_;
  return 1./n_polys_in_phase;
}

int
NodesVariablesPhaseBased::GetNumberOfPrevPolynomialsInPhase (int poly_id) const
{
  return polynomial_info_.at(poly_id).poly_in_phase_;
}

bool
NodesVariablesPhaseBased::IsConstantNode (int node_id) const
{
  bool is_constant = false;

  // node is considered constant if either left or right polynomial
  // belongs to a constant phase
  for (int poly_id : GetAdjacentPolyIds(node_id))
    if (IsInConstantPhase(poly_id))
      is_constant = true;

  return is_constant;
}

bool
NodesVariablesPhaseBased::AdjacentPolyIdsConstant (int node_id) const
{
  bool is_constant = true;

  // node is considered constant if left AND right polynomial
  // belongs to a constant phase
  std::vector<int> polys = GetAdjacentPolyIds(node_id);
  if (polys.size() > 1) {
    for (int poly_id : polys) {
	  is_constant = is_constant & IsInConstantPhase(poly_id);
    }
  }
  else
	  is_constant = false;

  return is_constant;
}

bool
NodesVariablesPhaseBased::IsInConstantPhase(int poly_id) const
{
  return polynomial_info_.at(poly_id).is_constant_;
}

NodesVariablesPhaseBased::NodeIds
NodesVariablesPhaseBased::GetIndicesOfNonConstantNodes() const
{
  NodeIds node_ids;

  for (int id=0; id<GetNodes().size(); ++id)
    if (!IsConstantNode(id))
      node_ids.push_back(id);

  return node_ids;
}

int
NodesVariablesPhaseBased::GetPhase (int node_id) const
{
//  assert(!IsConstantNode(node_id)); // because otherwise it has two phases

  int poly_id = GetAdjacentPolyIds(node_id).front();
  return polynomial_info_.at(poly_id).phase_;
}

int
NodesVariablesPhaseBased::GetPolyIDAtStartOfPhase (int phase) const
{
  int poly_id=0;
  for (int i=0; i<polynomial_info_.size(); ++i)
    if (polynomial_info_.at(i).phase_ == phase)
      return i;
}

Eigen::Vector3d
NodesVariablesPhaseBased::GetValueAtStartOfPhase (int phase) const
{
  int node_id = GetNodeIDAtStartOfPhase(phase);
  return GetNodes().at(node_id).p();
}

int
NodesVariablesPhaseBased::GetNodeIDAtStartOfPhase (int phase) const
{
  int poly_id=GetPolyIDAtStartOfPhase(phase);
  return GetNodeId(poly_id, Side::Start);
}

std::vector<int>
NodesVariablesPhaseBased::GetAdjacentPolyIds (int node_id) const
{
  std::vector<int> poly_ids;
  int last_node_id = GetNodes().size()-1;

  if (node_id==0)
    poly_ids.push_back(0);
  else if (node_id==last_node_id)
    poly_ids.push_back(last_node_id-1);
  else {
    poly_ids.push_back(node_id-1);
    poly_ids.push_back(node_id);
  }

  return poly_ids;
}

NodesVariablesPhaseBased::PolyInfo::PolyInfo(int phase, int poly_id_in_phase,
                               int num_polys_in_phase, bool is_constant, bool is_driving)
    :phase_(phase),
     poly_in_phase_(poly_id_in_phase),
     n_polys_in_phase_(num_polys_in_phase),
     is_constant_(is_constant),
	 is_driving_(is_driving)
{
}

void
NodesVariablesPhaseBased::SetNumberOfVariables(int n_variables)
{
  bounds_ = VecBound(n_variables, ifopt::NoBound);
  SetRows(n_variables);
}

void
NodesVariablesPhaseBased::PrintOptIdxMap ()
{
  for(auto it = index_to_node_value_info_.cbegin(); it != index_to_node_value_info_.cend(); ++it)
  {
	  std::cout << it->first << ": ";
	  auto nodes = it->second;
	  for(auto nvi : nodes) {
		  std::cout << "node " << nvi.id_ << ", dim: " << (char)(88+nvi.dim_) << ", deriv: " << nvi.deriv_;
		  std::cout << ", IsConstantNode? " << IsConstantNode(nvi.id_) << std::endl;
	  }
  }
}

NodesVariablesEEMotion::NodesVariablesEEMotion(int phase_count,
                                               bool is_in_contact_at_start,
                                               const std::string& name,
											   std::vector<int> n_polys_per_phase,
											   bool y_lateral_constraint)
    :NodesVariablesPhaseBased(phase_count,
                              is_in_contact_at_start,
                              name, n_polys_per_phase)
{
  y_lateral_constraint_ = y_lateral_constraint;

  index_to_node_value_info_ = GetPhaseBasedEEParameterization();
//  PrintOptIdxMap();
  SetNumberOfVariables(index_to_node_value_info_.size());
}

NodesVariablesEEForce::OptIndexMap
NodesVariablesEEMotion::GetPhaseBasedEEParameterization ()
{
  OptIndexMap index_map;

  int idx = 0; // index in variables set
  for (int node_id=0; node_id<nodes_.size(); ++node_id) {
    // swing node:
    if (!IsConstantNode(node_id)) {
      for (int dim=0; dim<GetDim(); ++dim) {
        // both position and velocity is optimized during swing phase
        index_map[idx++].push_back(NodeValueInfo(node_id, kPos, dim));
        index_map[idx++].push_back(NodeValueInfo(node_id, kVel, dim));
      }
    }
    // stance node
    else {
      // both position and velocity is optimized during stance phase (hybrid motion)
      int phase = GetPhase(node_id);
      bool is_first_node_in_phase = (node_id == GetNodeIDAtStartOfPhase(phase));

      for (int dim=0; dim<GetDim(); ++dim) {
        index_map[idx++].push_back(NodeValueInfo(node_id, kPos, dim));

		if (dim == Y) {
			if (y_lateral_constraint_) { // otherwise non-holonomic constraint will take care of this!
				nodes_.at(node_id).at(kVel).y() = 0.0;
			}
//			else
//				index_map[idx++].push_back(NodeValueInfo(node_id, kVel, dim));
		}
		else if (dim == Z) {
        	if (is_first_node_in_phase)
        		nodes_.at(node_id).at(kVel).z() = 0.0;
        }
        else {
          // wheel is allowed to move during stance!!
		  index_map[idx++].push_back(NodeValueInfo(node_id, kVel, dim));
        }

      }
    }
  }

  return index_map;
}

NodesVariablesEEForce::NodesVariablesEEForce(int phase_count,
                                              bool is_in_contact_at_start,
                                              const std::string& name,
											  std::vector<int> n_polys_per_phase)
    :NodesVariablesPhaseBased(phase_count,
                              !is_in_contact_at_start, // contact phase for force is non-constant
                              name, n_polys_per_phase)
{
  index_to_node_value_info_ = GetPhaseBasedEEParameterization();
//  PrintOptIdxMap ();
  SetNumberOfVariables(index_to_node_value_info_.size());
}

NodesVariablesEEForce::OptIndexMap
NodesVariablesEEForce::GetPhaseBasedEEParameterization ()
{
  OptIndexMap index_map;

  int idx = 0; // index in variables set
  for (int id=0; id<nodes_.size(); ++id) {
    // stance node: forces can be created during stance, so these nodes are optimized over.
    if (!IsConstantNode(id)) {
      for (int dim=0; dim<GetDim(); ++dim) {
        index_map[idx++].push_back(NodeValueInfo(id, kPos, dim));
        index_map[idx++].push_back(NodeValueInfo(id, kVel, dim));
      }
    }
    // swing node
    else {
      // forces can't exist during swing phase, so no need to be optimized
      // -> all node values simply set to zero.
      nodes_.at(id).at(kPos).setZero();
      nodes_.at(id).at(kVel).setZero();
    }
  }

  return index_map;
}

} /* namespace towr */
