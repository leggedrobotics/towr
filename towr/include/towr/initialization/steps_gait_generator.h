//
// Created by Vuk Pajovic on 30/11/2020.
//

#ifndef TOWR_INITIALIZATION_STEPS_GAIT_GENERATOR_H
#define TOWR_INITIALIZATION_STEPS_GAIT_GENERATOR_H

#include <memory>

#include <towr/terrain/examples/height_map_examples.h>
#include <towr/initialization/gait_generator.h>


namespace towr{

/**
 * @brief Produces the contact sequence for a case when there is a terrain made of steps
 */

class StepsGaitGenerator : public GaitGenerator {
public:
  StepsGaitGenerator(double des_v_x, std::shared_ptr<towr::HeightMap> terrain);
  virtual ~StepsGaitGenerator () = default;

  /**
   * This function makes duration of the phases based on the positions of the steps and desired velocity(goal/total_time)
   * @param T total time
   * @param ee end-effector
   * @return
   */
  GaitGenerator::VecTimes GetPhaseDurations(double T, EE ee) const override;

  /**
   * Dummy function just to override the superclass method
   * @param combo
   */
  void SetCombo(Combos combo) override;

  /**
   * Always returning true, because for the steps kind of terrain, there should be all contact start
   * @param ee
   * @return
   */
  bool IsInContactAtStart(EE ee) const;
private:
  /**
   * parameters for stepping
   */
   /// hard coded margins where should start the first swing phase based on the step in front related to the base position
  //std::array<double,4> step_margins_ = {0.3git stat7499999999999,0.17499999999,-0.3,-0.5};
  std::array<double,4> step_margins_ = {0.275,0.075,-0.3,-0.5};

  /// how long should be the swing phase over the step[s*1e2]
  int swing_phase_duration_=45;
  /**
   * simulation parameters used to generate durations
   */
  std::shared_ptr<towr::Steps> terrain_;
  double des_v_x_;

  ///only used for dummy purposes in GetGait() method
  ContactState BB_;

  ///only overriding the virtual method of the superclass
  GaitGenerator::GaitInfo GetGait(Gaits gait) const;
};

} /* namespace towr */
#endif // TOWR_INITIALIZATION_STEPS_GAIT_GENERATOR_H
