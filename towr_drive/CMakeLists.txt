cmake_minimum_required(VERSION 2.8.3)
project(towr_drive)

add_compile_options(-std=c++11)

find_package(catkin REQUIRED COMPONENTS
  roscpp
  roslib
  rosbag
  std_msgs
  towr_ros
  towr
  xpp_states
  xpp_msgs
  anymal_wheels_ctrl_track_msgs
  yaml_tools  
)

###################################
## catkin specific configuration ##
###################################
## The catkin_package macro generates cmake config files for your package
catkin_package(
   INCLUDE_DIRS
   LIBRARIES
   CATKIN_DEPENDS roscpp xpp_states
   DEPENDS towr
)

###########
## Build ##
###########
include_directories(
  include
  ${catkin_INCLUDE_DIRS}
)

# The interface to derive from for ros integration
add_library(${PROJECT_NAME} SHARED
  src/towr_drive_ros_interface.cc
)
add_dependencies(${PROJECT_NAME} 
  ${PROJECT_NAME}_gencpp
)
target_link_libraries(${PROJECT_NAME}
  ${catkin_LIBRARIES}
)

## An example executable that runs TOWR
add_executable(towr_drive_ros_app 
  src/towr_drive_ros_app.cc
)
add_dependencies(towr_drive_ros_app 
  ${PROJECT_NAME}_gencpp
)
target_link_libraries(towr_drive_ros_app
  ${catkin_LIBRARIES}
  ${PROJECT_NAME}
)

# Run optimization tests
add_executable(anymal_wheels_test 
  test/anymal_wheels_test.cc
)
target_link_libraries(anymal_wheels_test
  ${catkin_LIBRARIES}
)

# Send trajectory to controller from TOWR rosbag
add_executable(send_trajectory_from_bag 
  src/send_trajectory_from_bag.cc
)
target_link_libraries(send_trajectory_from_bag
  ${catkin_LIBRARIES}
)

# Extract geometric messges from TOWR rosbag (for matlab)
add_executable(extract_geom_msg_rosbag
  src/extract_geom_msg_rosbag.cc
)
target_link_libraries(extract_geom_msg_rosbag
  ${catkin_LIBRARIES}
)

#############
## Install ##
#############
# Mark library for installation
install(
  TARGETS towr_drive_ros_app
          ${PROJECT_NAME}
  ARCHIVE DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  LIBRARY DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  RUNTIME DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
)

# Mark header files for installation
install(
  DIRECTORY include/${PROJECT_NAME}/
  DESTINATION ${CATKIN_PACKAGE_INCLUDE_DESTINATION}
  FILES_MATCHING PATTERN "*.h"
)

# Mark other files for installation
install(
  DIRECTORY launch rviz bash
  DESTINATION ${CATKIN_PACKAGE_SHARE_DESTINATION}
)

