cmake_minimum_required(VERSION 2.8.3)
project(ros_common_robin_msgs)

find_package(catkin REQUIRED COMPONENTS
  genmsg
  actionlib
  actionlib_msgs
  std_msgs
  geometry_msgs
)


add_service_files(
   FILES
   get_object_pose.srv
   add_marker_srv.srv
   remove_marker_srv.srv
   localization.srv
   execute.srv
 )
 


 generate_messages(
   DEPENDENCIES
   std_msgs
   actionlib_msgs
   geometry_msgs
 )


catkin_package(
   CATKIN_DEPENDS std_msgs message_runtime actionlib actionlib_msgs std_msgs geometry_msgs

)

include_directories(
  ${catkin_INCLUDE_DIRS}
)

