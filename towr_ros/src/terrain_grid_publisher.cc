//
// Created by chaoni on 8/25/20.
//

#include <ros/ros.h>
#include <rosbag/bag.h>
#include <grid_map_ros/grid_map_ros.hpp>
#include <grid_map_msgs/GridMap.h>

#include <Eigen/Dense>

#include <towr_ros/TowrCommand.h>
#include <towr_ros/topic_names.h>
#include <towr/terrain/height_map.h>


namespace towr {

//static ros::Publisher rviz_pub;
//
//void UserCommandCallback(const towr_ros::TowrCommand& msg_in)
//{
////  // get which terrain
////  // save the grid_map
//  grid_map::GridMap map({"elevation"});
//  map.setFrameId("map");
//  map.setGeometry(grid_map::Length(4.0,4.0),0.03);
////  ROS_INFO("Created map with size %f x %f m (%i x %i cells).",
////           map.getLength().x(), map.getLength().y(),
////           map.getSize()(0), map.getSize()(1));
//
//  ros::Rate rate(30.0);
//  while (towr::n.ok()){
//
//    ros::Time time = ros::Time::now();
//    for (grid_map::GridMapIterator it(map); !it.isPastEnd(); ++it){
//      grid_map::Position position;
//      map.getPosition(*it, position);
//      if (position.x()>0.8)
//        map.at("elevation", *it) = 0.05;
//      else
//        map.at("elevation", *it) = 0;
//    }
////
////    map.setTimestamp(time.toNSec());
////    grid_map_msgs::GridMap message;
////    grid_map::GridMapRosConverter::toMessage(map, message);
////    rviz_pub.publish(message);
////
////  }
////
//////  ::ros::Time t0(1e-6);
//////  rosbag::Bag bag;
//////  bag.open("test.bag", rosbag::bagmode::Write);
//////  bag.write("/towr/grid_info", t0, message);
//}
////
} // namespace towr

int main(int argc, char *argv[])
{
  ros::init(argc, argv, "terrain_grid_publisher");

  ros::NodeHandle n("~");
  ros::Publisher rviz_pub;

//  ros::Subscriber goal_sub;
//  goal_sub       = n.subscribe(towr_msgs::user_command, 1, towr::UserCommandCallback);
  rviz_pub = n.advertise<grid_map_msgs::GridMap>("grid_map", 1, true);

  grid_map::GridMap map({"elevation"});
  map.setFrameId("world");
  map.setGeometry(grid_map::Length(3.0, 1.6), 0.03, grid_map::Position(0.8,0));
  ROS_INFO("Created map with size %f x %f m (%i x %i cells).",
           map.getLength().x(), map.getLength().y(),
           map.getSize()(0), map.getSize()(1));

  ros::Rate rate(30.0);
  while (n.ok()){

    ros::Time time = ros::Time::now();
    for (grid_map::GridMapIterator it(map); !it.isPastEnd(); ++it){
      grid_map::Position position;
      map.getPosition(*it, position);
      if (position.x()>0.8)
        map.at("elevation", *it) = 0.05;
      else
        map.at("elevation", *it) = 0;
    }

    map.setTimestamp(time.toNSec());
    grid_map_msgs::GridMap message;
    grid_map::GridMapRosConverter::toMessage(map, message);
    rviz_pub.publish(message);
    ROS_INFO_THROTTLE(1.0, "Grid map (timestamp %f) published.", message.info.header.stamp.toSec());

    rate.sleep();

  }
//  ros::spin();
  return 0;
}
