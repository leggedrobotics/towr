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

#ifndef TOWR_TOWR_ROS_INCLUDE_TOWR_ROS_HEIGHT_MAP_EXAMPLES_H_
#define TOWR_TOWR_ROS_INCLUDE_TOWR_ROS_HEIGHT_MAP_EXAMPLES_H_

#include <cmath>
#include <iostream>
#include <towr/terrain/height_map.h>

namespace towr {

/**
 * @addtogroup Terrains
 * @{
 */

/**
 * @brief Sample terrain of even height.
 */
class FlatGround : public HeightMap {
public:
  FlatGround(double height = 0.0);
  double GetHeight(double x, double y)  const override { return height_; };

private:
  double height_; // [m]
};

/**
 * @brief Sample terrain with a step in height in x-direction.
 */
class Block : public HeightMap {
public:
  double GetHeight(double x, double y)  const override;
  double GetHeightDerivWrtX(double x, double y) const override;

private:
  double block_start = 0.7;
  double length_     = 3.5;
  double height_     = 0.5; // [m]

  double eps_ = 0.03; // approximate as slope
  const double slope_ = height_/eps_;
};

/**
 * @brief Sample terrain with a two-steps in height in x-direction.
 */
class Stairs : public HeightMap {
public:
  double GetHeight(double x, double y) const override;

private:
  double first_step_start_  = 1.0;
  double first_step_width_  = 0.4;
  double height_first_step  = 0.2;
  double height_second_step = 0.4;
  double width_top = 1.0;
};

/**
 * @brief Sample terrain with parabola-modeled gap in x-direction.
 */
class Gap : public HeightMap {
public:
  double GetHeight(double x, double y) const override;
  double GetHeightDerivWrtX(double x, double y) const override;
  double GetHeightDerivWrtXX(double x, double y) const override;

private:
  const double gap_start_ = 0.5; // 1.0;
  const double w = 2.0; //0.5;
  const double h = 0.5; //1.5;

  const double slope_ = h/w;
  const double dx = w/2.0;
  const double xc = gap_start_ + dx; // gap center
  const double gap_end_x = gap_start_ + w;

  const double h_offset_ = 0.5;

  // generated with matlab
  // see matlab/gap_height_map.m
  // coefficients of 2nd order polynomial
  // h = a*x^2 + b*x + c
  const double a = (4*h)/(w*w);
  const double b = -(8*h*xc)/(w*w);
  const double c = -(h*(w - 2*xc)*(w + 2*xc))/(w*w);
};

/**
 * @brief Sample terrain with an increasing and then decreasing slope in x-direction.
 */
class Slope : public HeightMap {
public:
  double GetHeight(double x, double y) const override;
  double GetHeightDerivWrtX(double x, double y) const override;

private:
  const double slope_start_  = 1.0;
  const double up_length_    = 1.0;
  const double down_length_  = 1.0;
  const double height_center = 0.7;

  const double x_down_start_ = slope_start_+up_length_;
  const double x_flat_start_ = x_down_start_ + down_length_;
  const double slope_ = height_center/up_length_;
};

/**
 * @brief Sample terrain with a tilted vertical wall to cross a gap.
 */
class Chimney : public HeightMap {
public:
  double GetHeight(double x, double y) const override;
  double GetHeightDerivWrtY(double x, double y) const override;

private:
  const double x_start_ = 1.0;
  const double length_  = 1.5;
  const double y_start_ = 0.5; // distance to start of slope from center at z=0
  const double slope_   = 3.0;

  const double x_end_ = x_start_+length_;
};

/**
 * @brief Sample terrain with two tilted vertical walls to cross a gap.
 */
class ChimneyLR : public HeightMap {
public:
  double GetHeight(double x, double y) const override;
  double GetHeightDerivWrtY(double x, double y) const override;

private:
  const double x_start_ = 0.5;
  const double length_  = 1.0;
  const double y_start_ = 0.5; // distance to start of slope from center at z=0
  const double slope_   = 2;

  const double x_end1_ = x_start_+length_;
  const double x_end2_ = x_start_+2*length_;
};

/** @}*/

/**
 * @brief Sample terrain with a step in height with a linear transition.
 */
class Step : public HeightMap {
public:
  double GetHeight(double x, double y) const override;
  double GetHeightDerivWrtX(double x, double y) const override;
  double GetHeightDerivWrtXX(double x, double y) const override;

private:
  double slope_		 = std::tan(63*M_PI/180.0);  // std::tan(45*M_PI/180.0);
  double height_     = 0.2;
  double step_start_ = 1.0;
  double step_end_   = step_start_ + height_/slope_;

//  Eigen::Vector4d coeff {-400, 1260, -1320, 460}; // block from 1.0 to 1.1
};

/**
 * @brief Sample terrain with two steps in height with a linear transition.
 */
class TwoStep : public HeightMap {
public:
  double GetHeight(double x, double y) const override;
  double GetHeightDerivWrtX(double x, double y) const override;

private:
  double step_start = 1.0;
  double step_end   = 1.1;
  double height     = 0.2;

  double dist_steps	= 0.5;
  double slope = 2;
};

/**
 * @brief Sample terrain with one small slope on the left and another on the right.
 */
class TwoSlope : public HeightMap {
public:
  double GetHeight(double x, double y) const override;
  double GetHeightDerivWrtX(double x, double y) const override;
  double GetHeightDerivWrtXX(double x, double y) const override;

private:
  double dx = 0.025; //0.0326352;
  double dh = 0.05;  //0.05;

  double step_up_start   = 1.0;
  double step_up_end     = 1.1 - dx;
  double step_down_start = 1.2 - dx;
  double step_down_end   = 1.3 - 2*dx;
  double height          = 0.2 - dh;
  double dist_steps		 = 1.0;

  double slope = 2;
//  Eigen::Vector4d coeff {-400, 1260, -1320, 460}; // block from 1.0 to 1.1
//  Eigen::Vector4d coeff {-400, 60, 0, 0};
};

/**
 * @brief Sample terrain with an increasing and then decreasing slope.
 */
class SlopePlat : public HeightMap {
public:
  double GetHeight(double x, double y) const override;
  double GetHeightDerivWrtX(double x, double y) const override;

private:
  const double slope_start_   = 1.0;
  const double up_length_     = 0.091651513899117; //0.101905089898886; // 0.2;
  const double down_length_   = 0.2;
  const double plat_length_   = 1.216; //1.2; //1.0;
  const double height_center_ = 0.2;
  const double slope_up_	  = height_center_/up_length_;
  const double slope_down_	  = height_center_/down_length_;

  const double x_plat_start_ = slope_start_ + up_length_;
  const double x_down_start_ = x_plat_start_ + plat_length_;
  const double x_flat_start_ = x_down_start_ + down_length_;
};

/**
 * @brief Sample terrain with an increasing and then decreasing slope.
 */
class MultipleSlopes : public HeightMap {
public:
  double GetHeight(double x, double y) const override;
  double GetHeightDerivWrtX(double x, double y) const override;

private:
  const double slope_start_   = 1.0;
  const double up_length_     = 0.091651513899117; //0.101905089898886;
  const double down_length_   = 0.2;
  const double plat_length_   = 1.216; //1.0;
  const double height_center_ = 0.2;
  const double slope_up_	  = height_center_/up_length_;
  const double slope_down_	  = height_center_/down_length_;

  const double x_plat_start_ = slope_start_ + up_length_;
  const double x_down_start_ = x_plat_start_ + plat_length_;
  const double x_flat_start_ = x_down_start_ + down_length_;

  const double dist_slopes_ = 6.0;
};

/**
 * @brief Sample terrain with a low frequency sine profile.
 */
class SineLowFreq : public HeightMap {
public:
  double GetHeight(double x, double y) const override;
  double GetHeightDerivWrtX(double x, double y) const override;
  double GetHeightDerivWrtXX(double x, double y) const override;

private:
  const double sine_start_ = 0.5;
  const double freq_ = 2.0;
  const double amp_  = 0.2;
  const double h_offset_ = amp_;
  const double n_cycles_ = 2.0;
  const double sine_end_ = n_cycles_*2*M_PI/freq_ + sine_start_;
};

/**
 * @brief Sample terrain with a low frequency sine profile.
 */
class SineHighFreq : public HeightMap {
public:
  double GetHeight(double x, double y) const override;
  double GetHeightDerivWrtX(double x, double y) const override;
  double GetHeightDerivWrtXX(double x, double y) const override;

private:
  const double sine_start_ = 0.5;
  const double freq_ = 2.0*M_PI/0.7;
  const double amp_  = 0.06;
  const double h_offset_ = amp_;
  const double n_cycles_ = 3.0;
  const double sine_end_ = n_cycles_*2*M_PI/freq_ + sine_start_;
};

/**
 * @brief Sample terrain with a slope and oscillations.
 */
class Rough : public HeightMap {
public:
  double GetHeight(double x, double y) const override;
  double GetHeightDerivWrtX(double x, double y) const override;
  double GetHeightDerivWrtXX(double x, double y) const override;

private:
  const double rough_start_ = 0.5;
  const double freq_  = 5.0;
  const double amp_   = 0.1;
  const double slope_ = 0.2;
  const double n_cycles_ = 2.0;
  const double rough_end_ = n_cycles_*2*M_PI/freq_ + rough_start_;
  const double h_end_ = amp_*sin(freq_*(rough_end_-rough_start_))+slope_*(rough_end_-rough_start_);
};

} /* namespace towr */

#endif /* TOWR_TOWR_ROS_INCLUDE_TOWR_ROS_HEIGHT_MAP_EXAMPLES_H_ */
