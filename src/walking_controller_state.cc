/**
 @file    walking_controller_state.cc
 @author  Alexander W. Winkler (winklera@ethz.ch)
 @date    Jun 30, 2016
 @brief   Brief description
 */

#include <xpp/exe/walking_controller_state.h>
#include <xpp/exe/walking_controller.h>

namespace xpp {
namespace exe {

WalkingControllerState::WalkingControllerState ()
{
  // TODO Auto-generated constructor stub
}

WalkingControllerState::~WalkingControllerState ()
{
  // TODO Auto-generated destructor stub
}

WalkingControllerState::StatesMap
WalkingControllerState::BuildStates ()
{
  StatesMap states;

  states.emplace(kFirstPlanning,      std::make_shared<FirstPlanning>());
  states.emplace(kExecuting,          std::make_shared<Executing>());
  states.emplace(kUpdateAndExecuting, std::make_shared<UpdateAndExecuting>());
  states.emplace(kSleeping,           std::make_shared<Sleeping>());

  return states;
}

void FirstPlanning::Run(WalkingController* context) const
{
  context->EstimateCurrPose();
  context->PublishCurrentState();
  context->SetState(kSleeping); // waiting for nlp optimizer to finish
}

void Sleeping::Run(WalkingController* context) const
{
  if (context->optimal_trajectory_updated) {
    context->ResetTime();
    context->SetState(kUpdateAndExecuting);
    context->optimal_trajectory_updated = false;
  }
}

void UpdateAndExecuting::Run(WalkingController* context) const
{
  context->IntegrateOptimizedTrajectory();
  context->SetState(kExecuting);
  context->ExecuteLoop(); // to send values to the controller in every controller iteration
}

void Executing::Run(WalkingController* context) const
{
  if (context->EndCurrentExecution())
    context->SetState(kSleeping);

  if (context->IsTimeToSendOutState())
    context->PublishOptimizationStartState();

  context->ExecuteLoop();
}



} /* namespace exe */
} /* namespace xpp */
