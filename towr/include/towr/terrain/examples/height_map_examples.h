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
  double height_     = 0.1; // [m]

  double eps_ = 0.03; // approximate as slope
  const double slope_ = height_/eps_;
};

/**
 * @brief Sample terrain with a step in height in the right side of the path.
 */
class BlockRight : public HeightMap {
public:
  double GetHeight(double x, double y)  const override;
  double GetHeightDerivWrtX(double x, double y) const override;

private:
  double block_start = 0.7;
  double length_     = 3.5;
  double height_     = 0.2; // [m]

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
  const double gap_start_ = 0.7;
  const double w = 0.25;
  const double h = 1.5;

  const double slope_ = h/w;
  const double dx = w/2.0;
  const double xc = gap_start_ + dx; // gap center
  const double gap_end_x = gap_start_ + w;

  const double height_offset = 0.0; //0.2;

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
  const double slope_start_ = 0.5;
  const double up_length_   = 2.0*cos(0.3);
  const double down_length_ = 2.0*cos(0.3);
  const double height_center = 2.0*sin(0.3);

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

/**
 * @brief Sample terrain with multiple steps in x-direction.
 */
class Steps : public HeightMap{
public:
  Steps();
  double GetHeight(double x, double y) const override;
  double GetHeightDerivWrtX(double x, double y) const override;


  /**
   * Each step is defined as the absolute height of the step and from where on the x-axis it starts, on the same place in the vector. First
   * step is always 0, as the ground step and it starts from -inf, which is here presented as -1000. This means that the number of steps will
   * be num_of_steps=level_heigths.size()-1. Adding numbers to these two vectors will change number of steps.
   */

  std::vector<double> level_heights_={0.0,0.1};
  // terrain step-up-and-down
  // 1
  //   std::vector<double> level_heights_={0.0,0.155,0.0};
  // 2
  //   std::vector<double> level_heights_={0.0,0.155,0.0};
  // 3
  //   std::vector<double> level_heights_={0.0,0.2,0.0};
  // 4
  //   std::vector<double> level_heights_={0.0,0.2,0.0};
  // terrain double-step
  // 1
  //   std::vector<double> level_heights_={0.0,0.2, 0.29};
  // 2
  //   std::vector<double> level_heights_={0.0,0.2, 0.29};
  // 3
  //   std::vector<double> level_heights_={0.0,0.2, 0.39};
  // 4
  //   std::vector<double> level_heights_={0.0,0.2, 0.39};
  // 5
  //   std::vector<double> level_heights_={0.0,0.2, 0.35};
  // 6
  //   std::vector<double> level_heights_={0.0,0.2, 0.35};




  std::vector<double> level_starts_={-1000.0,0.9};
  // terrain step-up-and-down
  // 1
  //   std::vector<double> level_starts_={-1000.0,0.9,2.415};
  // 2
  //   std::vector<double> level_starts_={-1000.0,0.9,1.915};
  // 3
  //   std::vector<double> level_starts_={-1000.0,0.9, 2.105};
  // 4
  //   std::vector<double> level_starts_={-1000.0,0.9, 1.7};
  // terrain double-step
  // 1
  // std::vector<double> level_starts_={-1000.0,0.9, 2.105};
  // 2
  // std::vector<double> level_starts_={-1000.0,0.9, 1.7};
  // 3
  // std::vector<double> level_starts_={-1000.0,0.9, 2.105};
  // 4
  // std::vector<double> level_starts_={-1000.0,0.9, 1.7};
  // 5
  // std::vector<double> level_starts_={-1000.0,0.9, 2.105};
  // 6
  // std::vector<double> level_starts_={-1000.0,0.9, 1.7};

private:

  /***
   * In order not to have discontinuities, steep slope is introduced between steps.
   */
  double slope_length_=0.02;
  std::vector<double> level_slope_ends_;
  std::vector<double> slopes_;
};
/** @}*/



} /* namespace towr */

#endif /* TOWR_TOWR_ROS_INCLUDE_TOWR_ROS_HEIGHT_MAP_EXAMPLES_H_ */
