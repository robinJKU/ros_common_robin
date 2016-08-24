cmake_minimum_required(VERSION 2.8.3)
project(ros_common_robin_driver)

find_package(catkin REQUIRED
	roscpp
	tf
	controller_manager
	hardware_interface
)

catkin_package(
	
)

add_executable(robin_odometry_publisher_simulation_only src/robin_odometry_publisher_simulation_only.cpp)
target_link_libraries(robin_odometry_publisher_simulation_only ${catkin_LIBRARIES})

add_executable(robin_base_interface src/robin_base_interface.cpp)
target_link_libraries(robin_base_interface ${catkin_LIBRARIES})
