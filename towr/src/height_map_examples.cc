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

#include <towr/terrain/examples/height_map_examples.h>

namespace towr {


FlatGround::FlatGround(double height)
{
  height_ = height;
}

// BLOCK
double
Block::GetHeight (double x, double y) const
{
  double h = 0.0;

  // very steep ramp leading up to block
  if (block_start <= x && x <=block_start+eps_)
    h = slope_*(x-block_start);

  if (block_start+eps_ <= x && x <= block_start+length_)
    h = height_;

  return h;
}

double
Block::GetHeightDerivWrtX (double x, double y) const
{
  double dhdx = 0.0;

  // very steep ramp leading up to block
  if (block_start <= x && x <=block_start+eps_)
    dhdx = slope_;

  return dhdx;
}

// BLOCK RIGHT
double
BlockRight::GetHeight (double x, double y) const
{
  double h = 0.0;

  if (y <= 0) {
    // very steep ramp leading up to block
    if (block_start <= x && x <=block_start+eps_)
	  h = slope_*(x-block_start);

    if (block_start+eps_ <= x && x <= block_start+length_)
	  h = height_;
  }

  return h;
}

double
BlockRight::GetHeightDerivWrtX (double x, double y) const
{
  double dhdx = 0.0;

  if (y <= 0) {
    // very steep ramp leading up to block
    if (block_start <= x && x <=block_start+eps_)
      dhdx = slope_;
  }

  return dhdx;
}

// STAIRS
double
Stairs::GetHeight (double x, double y) const
{
  double h = 0.0;

  if (x>=first_step_start_)
    h = height_first_step;

  if (x>=first_step_start_+first_step_width_)
    h = height_second_step;

  if (x>=first_step_start_+first_step_width_+width_top)
    h = 0.0;

  return h;
}


// GAP
double
Gap::GetHeight (double x, double y) const
{
  double h = height_offset;

  // modelled as parabola
  if (gap_start_ <= x && x <= gap_end_x)
    h = a*x*x + b*x + c + height_offset;

  return h;
}

double
Gap::GetHeightDerivWrtX (double x, double y) const
{
  double dhdx = 0.0;

  if (gap_start_ <= x && x <= gap_end_x)
    dhdx = 2*a*x + b;

  return dhdx;
}

double
Gap::GetHeightDerivWrtXX (double x, double y) const
{
  double dzdxx = 0.0;

  if (gap_start_ <= x && x <= gap_end_x)
    dzdxx = 2*a;

  return dzdxx;
}


// SLOPE
double
Slope::GetHeight (double x, double y) const
{
  double z = 0.0;
  if (x >= slope_start_)
    z = slope_*(x-slope_start_);

  // going back down
  if (x >= x_down_start_) {
    z = height_center - slope_*(x-x_down_start_);
  }

  // back on flat ground
  if (x >= x_flat_start_)
    z = 0.0;

  return z;
}

double
Slope::GetHeightDerivWrtX (double x, double y) const
{
  double dzdx = 0.0;
  if (x >= slope_start_)
    dzdx = slope_;

  if (x >= x_down_start_)
    dzdx = -slope_;

  if (x >= x_flat_start_)
    dzdx = 0.0;

  return dzdx;
}


// Chimney
double
Chimney::GetHeight (double x, double y) const
{
  double z = 0.0;

  if (x_start_<=x && x<=x_end_)
    z = slope_*(y-y_start_);

  return z;
}

double
Chimney::GetHeightDerivWrtY (double x, double y) const
{
  double dzdy = 0.0;

  if (x_start_<= x && x<= x_end_)
    dzdy = slope_;

  return dzdy;
}


// Chimney LR
double
ChimneyLR::GetHeight (double x, double y) const
{
  double z = 0.0;

  if (x_start_<=x && x<=x_end1_)
    z = slope_*(y-y_start_);

  if (x_end1_<=x && x<=x_end2_)
    z = -slope_*(y+y_start_);

  return z;
}

double
ChimneyLR::GetHeightDerivWrtY (double x, double y) const
{
  double dzdy = 0.0;

  if (x_start_ <= x && x <= x_end1_)
    dzdy = slope_;

  if (x_end1_<=x && x<=x_end2_)
    dzdy = -slope_;

  return dzdy;
}

// Steps

Steps::Steps(){
  level_slope_ends_.push_back(0);
  slopes_.push_back(0);
  for(int i=1;i<level_starts_.size();i++){
    level_slope_ends_.push_back(level_starts_.at(i)+slope_length_);
    slopes_.push_back((level_heights_.at(i)-level_heights_.at(i-1))/(level_slope_ends_.at(i)-level_starts_.at(i)));
  }
}
double
Steps::GetHeight(double x, double y) const {
  double h = 0.0;
  int index;
  for(index=0;index<level_starts_.size();index++){
    if(x < level_starts_.at(index)) break;
  }
  index-=1;

  if(index==0){
    h = level_heights_.at(index);
  }else{
    if (x>=level_starts_.at(index) && x<level_slope_ends_.at(index)) {
      h = slopes_.at(index)*(x-level_starts_.at(index))+level_heights_.at(index-1);
      if(h>level_heights_.at(index)) h=level_heights_.at(index);
      if(h<level_heights_.at(index-1)) h=level_heights_.at(index-1);
    }
    if (x>=level_slope_ends_.at(index)) {h = level_heights_.at(index);}
  }
  return h;
}

double
Steps::GetHeightDerivWrtX(double x, double y) const {
  double dhdx = 0.0;

  int index;
  for(index=0;index<level_starts_.size();index++){
    if(x < level_starts_.at(index)) break;
  }
  index-=1;

  if(index>0){
    if (x>=level_starts_.at(index) && x<level_slope_ends_.at(index)) {
      dhdx = slopes_.at(index);
    }
  }

  return dhdx;
}

} /* namespace towr */
