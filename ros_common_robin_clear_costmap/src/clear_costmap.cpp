/*
  clears costmap everytime the robot reaches a sent goal

*/

#include <ros/ros.h>
#include <tf/transform_broadcaster.h>
#include <actionlib/client/simple_action_client.h>
#include <actionlib_msgs/GoalStatus.h>
#include <move_base_msgs/MoveBaseAction.h>

#include "sensor_msgs/LaserScan.h"
#include <std_srvs/Empty.h>

//  global
int last_status; 	// save last status, should only reset costmaps, when status changed

//  declare functions

void goalCallback(const actionlib_msgs::GoalStatusArray::ConstPtr& goal_msg);


int main(int argc, char **argv) {
  ros::init(argc,argv, "clear_costmap");
  ros::NodeHandle nh1;
  last_status = 0;
  
  ros::Subscriber goal_status = nh1.subscribe("move_base/status", 10, goalCallback);
  			
  ros::spin();
  return 0;
}

void goalCallback(const actionlib_msgs::GoalStatusArray::ConstPtr& goal_msg) {
	actionlib_msgs::GoalStatusArray latest_goal = (*goal_msg);
	
	if (latest_goal.status_list.size() != 0 && latest_goal.status_list[0].status != last_status) {
		switch (latest_goal.status_list[0].status) {
			case 2: case 3: case 4: case 5: 
				std_srvs::Empty srv;
				if(ros::service::call("move_base/clear_costmaps", srv)) {	//service call to reset costmaps
        				ROS_ERROR("goal_status = %d, clearing", latest_goal.status_list[0].status);
				} else {
					ROS_ERROR("goal_status = %d, clearing failed", latest_goal.status_list[0].status);
				} break;
		}
		last_status = latest_goal.status_list[0].status;		// save acual status
	}
	
	// Terminal states:
	// case 2: The goal received a cancel request after it started executing and has since completed its execution
	// case 3: The goal was achieved successfully
	// case 4: The goal was aborted during execution due to some failure
	// case 5: The goal was rejected because the goal was unattainable or invalid 
}
