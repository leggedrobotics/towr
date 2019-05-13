cmake_minimum_required(VERSION 2.8.3)
project(towr_test)

add_compile_options(-std=c++11)

find_package(catkin REQUIRED COMPONENTS
  roscpp
  roslib
  std_msgs
  towr_ros
  towr
  xpp_states
  xpp_msgs
  yaml_tools  
)

###################################
## catkin specific configuration ##
###################################
## The catkin_package macro generates cmake config files for your package
catkin_package(
   INCLUDE_DIRS
   LIBRARIES
   CATKIN_DEPENDS std_msgs roscpp towr_ros xpp_states xpp_msgs roslib yaml_tools
   DEPENDS towr
)

###########
## Build ##
###########
include_directories(
  include
  ${catkin_INCLUDE_DIRS}
)

add_executable(anymal_wheels_test 
  test/anymal_wheels_test.cc
)
target_link_libraries(anymal_wheels_test
  ${catkin_LIBRARIES}
)

