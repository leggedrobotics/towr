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

// StepFlat
double
StepFlat::GetHeight(double x, double y) const {
  double h = 0.0;

  // use continuous slope, uncomment this line
  //  if (x>=step_start_)
  //    h = slope_*(x-step_start_);

  // discrete step case
  if (x>=step_start_)
    h = step_height_;

  return h;
}

double
StepFlat::GetHeightDerivWrtX(double x, double y) const {
  double dhdx = 0.0;
  // use continuous slope, uncomment the following code
  //  if (x>=step_start_)
  //    dhdx = slope_;
  //
  //  if (x>=step_end_)
  //    dhdx = 0.0;

  return dhdx;
}
// GAP
double
Gap::GetHeight (double x, double y) const
{
  double h = 0.0;

  // modelled as parabola
  if (gap_start_ <= x && x <= gap_end_x)
    h = a*x*x + b*x + c;

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

double
Step::GetHeight(double x, double y) const
{
  double h = 0.0;

  if (step_start_ <= x && x <= step_end_)
//    h = coeff(0)*pow(x,3)+coeff(1)*pow(x,2)+coeff(2)*x+coeff(3);
    h = slope_*(x-step_start_);

  if (step_end_ <= x)
    h = height_;

  return h;
}

double
Step::GetHeightDerivWrtX(double x, double y) const
{
  double dhdx = 0.0;

  if (step_start_ <= x && x <= step_end_)
    dhdx = slope_;

  return dhdx;
}

double
Step::GetHeightDerivWrtXX(double x, double y) const
{
  double Dhdx = 0.0;

  if (step_start_ <= x && x <= step_end_)
    Dhdx = 0.0;

  return Dhdx;
}

// TWO SLOPE
double
TwoSlope::GetHeight(double x, double y) const
{
  double h = 0.0;

  if (y >= 0) {
	if (x >= step_up_start)
	  h = slope*(x-step_up_start);

	if (x >= step_up_end)
	  h = height;

	if (x >= step_down_start)
	  h = height - slope*(x-step_down_start);

	if (x >= step_down_end)
	  h = 0.0;
  }

  if (y < 0) {
	if (x >= (step_up_start + dist_steps))
	  h = slope*(x-(step_up_start+dist_steps));

	if (x >= (step_up_end + dist_steps))
	  h = height;

	if (x >= (step_down_start + dist_steps))
	  h = height - slope*(x-(step_down_start+dist_steps));

	if (x >= (step_down_end + dist_steps))
	  h = 0.0;
  }

  return h;
}

double
TwoSlope::GetHeightDerivWrtX(double x, double y) const
{
  double dhdx = 0.0;

  if (y >= 0) {
    if (x >= step_up_start)
      dhdx = slope;

    if (x >= step_up_end)
	  dhdx = 0.0;

    if (x >= step_down_start)
      dhdx = -slope;

    if (x >= step_down_end)
	  dhdx = 0.0;
  }

  if (y < 0) {
	if (x >= (step_up_start + dist_steps))
	  dhdx = slope;

	if (x >= (step_up_end + dist_steps))
	  dhdx = 0.0;

	if (x >= (step_down_start + dist_steps))
	  dhdx = -slope;

	if (x >= (step_down_end + dist_steps))
	  dhdx = 0.0;
  }

  return dhdx;
}

double
TwoSlope::GetHeightDerivWrtXX(double x, double y) const
{
  return 0.0;
}

// TWO STEP
double
TwoStep::GetHeight(double x, double y) const
{
  double h = 0.0;

  if (x >= step_start)
	h = slope*(x-step_start);

  if (x >= step_end)
	h = height;

  if (x >= (step_start+dist_steps))
	h = height + slope*(x-(step_start+dist_steps));

  if (x >= (step_end+dist_steps))
	h = 2*height;

  return h;
}

double
TwoStep::GetHeightDerivWrtX(double x, double y) const
{
  double dhdx = 0.0;

  if (x >= step_start)
	dhdx = slope;

  if (x >= step_end)
    dhdx = 0.0;

  if (x >= (step_start+dist_steps))
	dhdx = slope;

  if (x >= (step_end+dist_steps))
	dhdx = 0.0;

  return dhdx;
}

// FIVE STEPS
double
FiveSteps::GetHeight(double x, double y) const
{
  double h = 0.0;

  for (int i = 0; i < num_steps; ++i) {
	if (x >= (step_start + i*dist_steps))
	  h = i*step_height + slope*(x-(step_start + i*dist_steps));
	if (x >= (step_start + i*dist_steps + step_width))
	  h = (i+1)*step_height;
  }

  return h;
}

double
FiveSteps::GetHeightDerivWrtX(double x, double y) const
{
  double dhdx = 0.0;

  for (int i = 0; i < num_steps; ++i) {
	if (x >= (step_start + i*dist_steps))
	  dhdx = slope;
	if (x >= (step_start + i*dist_steps + step_width))
	  dhdx = 0.0;
  }

  return dhdx;
}

// SLOPE PLAT
double
SlopePlat::GetHeight(double x, double y) const
{
  double h = 0.0;

//  if (y > 0) {
  // going up
  if (x >= slope_start_)
    h = slope_up_*(x-slope_start_);

  // flat platform
  if (x >= x_plat_start_)
    h = height_center_;

  // going back down
  if (x >= x_down_start_) {
    h = height_center_ - slope_down_*(x-x_down_start_);
  }

  // back on flat ground
  if (x >= x_flat_start_)
    h = 0.0;
//  }

  return h;
}

double
SlopePlat::GetHeightDerivWrtX(double x, double y) const
{
  double dhdx = 0.0;

//  if (y > 0) {
  if (x >= slope_start_)
    dhdx = slope_up_;

  if (x >= x_plat_start_)
	dhdx = 0.0;

  if (x >= x_down_start_)
    dhdx = -slope_down_;

  if (x >= x_flat_start_)
    dhdx = 0.0;
//  }

  return dhdx;
}

// MULTIPLE SLOPES
double
MultipleSlopes::GetHeight(double x, double y) const
{
  double h = 0.0;

  if (y <= 0)
  {
	if (x >= slope_start_)
      h = slope_up_*(x-slope_start_);

	if (x >= x_plat_start_)
      h = height_center_;

	if (x >= x_down_start_)
      h = height_center_ - slope_down_*(x-x_down_start_);

	if (x >= x_flat_start_)
      h = 0.0;
  }

  if (y >= 0)
  {
	if (x >= (slope_start_+dist_slopes_) )
      h = slope_up_*(x-(slope_start_+dist_slopes_));

	if (x >= (x_plat_start_+dist_slopes_) )
      h = height_center_;

	if (x >= (x_down_start_+dist_slopes_) )
      h = height_center_ - slope_down_*(x-(x_down_start_+dist_slopes_));

	if (x >= (x_flat_start_+dist_slopes_) )
      h = 0.0;
  }

  return h;
}

double
MultipleSlopes::GetHeightDerivWrtX(double x, double y) const
{
  double dhdx = 0.0;

  if (y <= 0)
  {
	if (x >= slope_start_)
	  dhdx = slope_up_;

	if (x >= x_plat_start_)
	  dhdx = 0.0;

	if (x >= x_down_start_)
	  dhdx = -slope_down_;

	if (x >= x_flat_start_)
	  dhdx = 0.0;
  }

  if (y >= 0)
  {
	if (x >= (slope_start_+dist_slopes_) )
	  dhdx = slope_up_;

	if (x >= (x_plat_start_+dist_slopes_) )
	  dhdx = 0.0;

	if (x >= (x_down_start_+dist_slopes_) )
	  dhdx = -slope_down_;

	if (x >= (x_flat_start_+dist_slopes_) )
	  dhdx = 0.0;
  }

  return dhdx;
}

// LOW FREQUENCY SINE
double
SineLowFreq::GetHeight(double x, double y) const
{
  double h = h_offset_;
  if (x >= sine_start_ && x <= sine_end_)
	h = amp_*sin(freq_*(x-sine_start_)) + h_offset_;

  return h;
}

double
SineLowFreq::GetHeightDerivWrtX(double x, double y) const
{
  double dhdx = 0.0;
  if (x >= sine_start_ && x <= sine_end_)
	dhdx = freq_*amp_*cos(freq_*(x-sine_start_));

  return dhdx;
}

double
SineLowFreq::GetHeightDerivWrtXX(double x, double y) const
{
  double Dhdx = 0.0;
  if (x >= sine_start_ && x <= sine_end_)
	Dhdx = -freq_*amp_*freq_*sin(freq_*(x-sine_start_));

  return Dhdx;
}

// HIGH FREQUENCY SINE
double
SineHighFreq::GetHeight(double x, double y) const
{
  double h = h_offset_;
  if (x >= sine_start_ && x <= sine_end_)
	h = amp_*sin(freq_*(x-sine_start_)) + h_offset_;

  return h;
}

double
SineHighFreq::GetHeightDerivWrtX(double x, double y) const
{
  double dhdx = 0.0;
  if (x >= sine_start_ && x <= sine_end_)
	dhdx = freq_*amp_*cos(freq_*(x-sine_start_));

  return dhdx;
}

double
SineHighFreq::GetHeightDerivWrtXX(double x, double y) const
{
  double Dhdx = 0.0;
  if (x >= sine_start_ && x <= sine_end_)
	Dhdx = -freq_*amp_*freq_*sin(freq_*(x-sine_start_));

  return Dhdx;
}

// ROUGH
double
Rough::GetHeight(double x, double y) const
{
  double h = 0.0;
  if (x >= rough_start_)
	h = amp_*sin(freq_*(x-rough_start_))+slope_*(x-rough_start_);

  if (x > rough_end_)
	h = h_end_;

  return h;
}

double
Rough::GetHeightDerivWrtX(double x, double y) const
{
  double dhdx = 0.0;
  if (x >= rough_start_)
	dhdx = freq_*amp_*cos(freq_*(x-rough_start_))+slope_;

  if (x > rough_end_)
	dhdx = 0.0;

  return dhdx;
}

double
Rough::GetHeightDerivWrtXX(double x, double y) const
{
  double Dhdx = 0.0;
  if (x >= rough_start_)
	Dhdx = -freq_*amp_*freq_*sin(freq_*(x-rough_start_));

  if (x > rough_end_)
	Dhdx = 0.0;

  return Dhdx;
}

// ROUNDSTAIR
    double
    RoundStair::GetHeight(double x, double y) const
    {
        double h = 0.0;

        if (y >= 0) {
            if (x >= step_start)
                h = slope * (x - step_start);

            if (x >= step_end)
                h = height;

        } else {
            if (x >= step_start1)
                h = slope1 * (x - step_start1);

            if (x >= step_end1)
                h = height1;
        }
        return h;
    }

    double
    RoundStair::GetHeightDerivWrtX(double x, double y) const
    {
        double dhdx = 0.0;

        if (y >= 0) {
        if (x >= step_start)
            dhdx = slope;

        if (x >= step_end)
            dhdx = 0.0;

    } else {
            if (x >= step_start1)
                dhdx = slope1;

            if (x >= step_end1)
                dhdx = 0.0;

        }

        return dhdx;
    }



} /* namespace towr */
