#include <towr/initialization/steps_gait_generator.h>
#include <iostream>

namespace towr{

StepsGaitGenerator::StepsGaitGenerator(double des_v_x, std::shared_ptr<towr::HeightMap> terrain){
  des_v_x_=des_v_x;
  terrain_=std::dynamic_pointer_cast<Steps>(terrain);
}

GaitGenerator::VecTimes
StepsGaitGenerator::GetPhaseDurations(double T, EE ee) const{


  // based on margin making touchdown_place;
  std::vector<double> touchdown_place;
  std::vector<double> steps_x_axis(terrain_->level_starts_.begin(), terrain_->level_starts_.end());

  for(int i=1;i<steps_x_axis.size();i++){
    touchdown_place.emplace_back(steps_x_axis.at(i)-step_margins_.at(ee));
  }

  std::vector<double> durations;


  // calculating the time of a touchdown, based on desired x-velocity of a base
  std::vector<int> touchdown_time;
  for(int i=0;i<touchdown_place.size();i++){
    touchdown_time.emplace_back(touchdown_place.at(i)/des_v_x_*100);
  }

  int number_of_steps=touchdown_time.size();


  // converting everything into simulation times
  int total_t_integer= (float)(T)*100;
  int last_step_touchdown_time=0;

  for(int i=0;i<number_of_steps;i++){
    durations.emplace_back(touchdown_time.at(i)-swing_phase_duration_-last_step_touchdown_time);
    if (durations.back() <= 0){
      std::cout<< "Error : Phase duration is zero or less than zero. Change the speed of the robot,step_margins or swing phase duration!!!";
    }
    durations.emplace_back(swing_phase_duration_);
    last_step_touchdown_time=touchdown_time.at(i);
  }
  durations.emplace_back(total_t_integer-last_step_touchdown_time);
  for(auto t=durations.begin(); t!=durations.end(); ++t){
    *t=*t/100;
  }
  return durations;
}

void StepsGaitGenerator::SetCombo(Combos combo) {}

bool StepsGaitGenerator::IsInContactAtStart(EE ee) const{
  return true;
}

StepsGaitGenerator::GaitInfo
StepsGaitGenerator::GetGait(Gaits gait) const {
  auto times =
      {
          0.6,
      };
  auto contacts =
      {
          BB_,
      };

  return std::make_pair(times, contacts);
}


} /* namespace towr */