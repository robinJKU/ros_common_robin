cmake_minimum_required(VERSION 2.8.3)
project(ros_common_robin_clear_costmap)


find_package(catkin REQUIRED COMPONENTS
  roscpp
  std_msgs
  tf
)

catkin_package(
	CATKIN_DEPENDS roscpp std_msgs tf
)

###########
## Build ##
###########

include_directories(
  ${catkin_INCLUDE_DIRS}
)


add_executable(clear_costmap src/clear_costmap.cpp)
target_link_libraries(clear_costmap ${catkin_LIBRARIES}
)
