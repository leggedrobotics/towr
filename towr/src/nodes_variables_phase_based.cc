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

#include <towr/initialization/quadruped_gait_generator.h>
#include <towr/initialization/gait_generator.h>

#include <towr/constraints/force_constraint.h>

#include <towr/variables/nodes_variables_phase_based.h>
#include <towr/variables/cartesian_dimensions.h>

#include <iostream>

namespace towr {


std::vector<NodesVariablesPhaseBased::PolyInfo>
BuildPolyInfos (int phase_count, bool first_phase_constant,
        int n_polys_in_first_drive_phase,
		 int n_polys_in_drift_phase,
		 int n_polys_in_second_drive_phase)
{
  using PolyInfo = NodesVariablesPhaseBased::PolyInfo;
  std::vector<PolyInfo> polynomial_info;

  bool phase_constant = first_phase_constant;	//first phase for all drive phase, so "constant"

  for (int i=0; i<phase_count; ++i) {

	  if (i == 0){
		  for (int j=0; j<n_polys_in_first_drive_phase; ++j)
		          polynomial_info.push_back(PolyInfo(i,j,n_polys_in_first_drive_phase, false));
	  }
	  if (i == 1){
		  for (int j=0; j<n_polys_in_drift_phase; ++j)
		  		  polynomial_info.push_back(PolyInfo(i,j,n_polys_in_drift_phase, false));
	  }
	  else {
		  for (int j=0; j<n_polys_in_second_drive_phase; ++j)
		  		  polynomial_info.push_back(PolyInfo(i,j,n_polys_in_second_drive_phase, false));
	  }

  }

  return polynomial_info;
}

NodesVariablesPhaseBased::NodesVariablesPhaseBased (int phase_count,
                                                    bool first_phase_constant,
                                                    const std::string& name,
													int n_polys_in_first_drive_phase,
													int n_polys_in_drift_phase,
													int n_polys_in_second_drive_phase,
													EE ee)
    :NodesVariables(name)
{
  polynomial_info_ = BuildPolyInfos(phase_count, first_phase_constant, n_polys_in_first_drive_phase, n_polys_in_drift_phase, n_polys_in_second_drive_phase);

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
NodesVariablesPhaseBased::GetNumberOfPrevPolynomialsInPhase(int poly_id) const
{
  return polynomial_info_.at(poly_id).poly_in_phase_;
}

bool
NodesVariablesPhaseBased::IsConstantNode (int node_id) const
{
	//no constant nodes if there is force and movement at the same time!
  bool is_constant = false;

  return is_constant;
}

bool
NodesVariablesPhaseBased::IsInConstantPhase(int poly_id) const
{
  return polynomial_info_.at(poly_id).is_constant_;
}

NodesVariablesPhaseBased::NodeIds
NodesVariablesPhaseBased::GetIndicesOfAllNodes() const
{
  NodeIds node_ids;

  for (int id=0; id<GetNodes().size(); ++id)
      node_ids.push_back(id);

  return node_ids;
}

NodesVariablesPhaseBased::NodeIds
NodesVariablesPhaseBased::GetIndicesNodesWOFirstAndLast() const
{
  NodeIds node_ids;

  for (int id=1; id<(GetNodes().size()-1); ++id)
      node_ids.push_back(id);

  return node_ids;
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
NodesVariablesPhaseBased::GetPhase (int node_id, EE ee) const
{

	int n_nodes_first_drive_phase = params_.polynomials_in_first_drive_phase_ + 1;
	int n_nodes_drift_phase = params_.polynomials_in_drift_phase_;
	int n_nodes_second_drive_phase = params_.polynomials_in_second_drive_phase_;


	if (params_.just_drive_){
		return 0;
	}
	else {
	if (node_id < n_nodes_first_drive_phase){
		if (ee == 0 or ee == 1)
			return 0;
		if (ee == 2 or ee == 3)
			return 3;
	}

	else if (node_id > (n_nodes_first_drive_phase + n_nodes_drift_phase - 1)){
		if (ee == 0 or ee == 1)
			return 2;
		if (ee == 2 or ee == 3)
			return 5;
	}
	else {
		if (ee == 0 or ee == 1)
					return 1;
		if (ee == 2 or ee == 3)
					return 4;
	}
	}
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
                               int num_polys_in_phase, bool is_constant)
    :phase_(phase),
     poly_in_phase_(poly_id_in_phase),
     n_polys_in_phase_(num_polys_in_phase),
     is_constant_(is_constant)
{
}

void
NodesVariablesPhaseBased::SetNumberOfVariables(int n_variables)
{
  bounds_ = VecBound(n_variables, ifopt::NoBound);
  SetRows(n_variables);
}

NodesVariablesEEMotion::NodesVariablesEEMotion(int phase_count,
                                               bool is_in_contact_at_start,
                                               const std::string& name,
											   int n_polys_in_first_drive_phase,
											   int n_polys_in_drift_phase,
											   int n_polys_in_second_drive_phase,
											   EE ee)
    :NodesVariablesPhaseBased(phase_count,
                              !is_in_contact_at_start, // contact phase for motion is !not! constant
                              name,
							  n_polys_in_first_drive_phase,
							  n_polys_in_drift_phase,
							  n_polys_in_second_drive_phase,
							  ee)
{
  index_to_node_value_info_ = GetPhaseBasedEEParameterization(ee);
  SetNumberOfVariables(index_to_node_value_info_.size());
}

NodesVariablesEEForce::OptIndexMap
NodesVariablesEEMotion::GetPhaseBasedEEParameterization (EE ee)
{
  OptIndexMap index_map;

  int idx = 0; // index in variables set
  for (int node_id=0; node_id<nodes_.size(); ++node_id) {
	 int phase = GetPhase(node_id, ee);

	 index_map[idx++].push_back(NodeValueInfo(node_id, kPos, X));
	 index_map[idx++].push_back(NodeValueInfo(node_id, kPos, Y));
	 index_map[idx++].push_back(NodeValueInfo(node_id, kVel, X));
	 index_map[idx++].push_back(NodeValueInfo(node_id, kVel, Y));
	 nodes_.at(node_id).at(kVel).z() = 0.0;
    }

  return index_map;
}

NodesVariablesEEForce::NodesVariablesEEForce(int phase_count,
                                              bool is_in_contact_at_start,
                                              const std::string& name,
											  int n_polys_in_first_drive_phase,
											  int n_polys_in_drift_phase,
											  int n_polys_in_second_drive_phase,
											  EE ee)
    :NodesVariablesPhaseBased(phase_count,
                              !is_in_contact_at_start, // contact phase for force is non-constant
                              name,
							  n_polys_in_first_drive_phase,
							  n_polys_in_drift_phase,
							  n_polys_in_second_drive_phase,
							  ee)
{
  index_to_node_value_info_ = GetPhaseBasedEEParameterization(ee);
  SetNumberOfVariables(index_to_node_value_info_.size());
}

NodesVariablesEEForce::OptIndexMap
NodesVariablesEEForce::GetPhaseBasedEEParameterization (EE ee)
{
  OptIndexMap index_map;

  int idx = 0; // index in variables set
  for (int id=0; id<nodes_.size(); ++id) {

      for (int dim=0; dim<GetDim(); ++dim) {
        index_map[idx++].push_back(NodeValueInfo(id, kPos, dim));
        index_map[idx++].push_back(NodeValueInfo(id, kVel, dim));

      }
  	}

  return index_map;
}

} /* namespace towr */
