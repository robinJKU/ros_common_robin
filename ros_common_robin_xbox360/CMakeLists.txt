cmake_minimum_required(VERSION 2.8.3)
project(ros_common_robin_xbox360)

catkin_package(
#  INCLUDE_DIRS include
#  LIBRARIES compute_cdist
#  CATKIN_DEPENDS roscpp visualization_msgs robot_model_helper compute_cdist
#  DEPENDS assimp Eigen mlpack
)


find_package(catkin REQUIRED COMPONENTS
  roscpp  
  geometry_msgs  
  sensor_msgs
)

include_directories(
  include
  ${catkin_INCLUDE_DIRS}
)

add_executable(xbox360_teleop_joy src/teleop_joy.cpp)
target_link_libraries(xbox360_teleop_joy ${catkin_LIBRARIES})
