cmake_minimum_required(VERSION 2.8.3)
project(ros_common_robin_tools)

## Find catkin macros and libraries
## if COMPONENTS list like find_package(catkin REQUIRED COMPONENTS xyz)
## is used, also find other catkin packages
find_package(catkin REQUIRED 
	geometry_msgs 
	tf
)

catkin_package(
  INCLUDE_DIRS include
)

include_directories(
  include
  ${catkin_INCLUDE_DIRS}
)

install(DIRECTORY include/${PROJECT_NAME}/
 DESTINATION ${CATKIN_PACKAGE_INCLUDE_DESTINATION}
)

add_executable(ros_common_robin_testsuite src/ros_common_robin_testsuite.cpp)
target_link_libraries(ros_common_robin_testsuite ${catkin_LIBRARIES})

