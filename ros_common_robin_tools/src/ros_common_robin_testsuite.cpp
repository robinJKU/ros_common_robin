#include <ros/ros.h>
#include <geometry_msgs/Twist.h>
#include <std_srvs/Empty.h>
#include <ros_common_robin_tools/common_cpp_tools.h>


class TESTSUITE{
	private:
		geometry_msgs::Twist cmd_vel;
		
		// advertised services
		ros::ServiceServer test90_service;
				
		// published topics
		ros::Publisher cmd_vel_publisher;

	public:
		TESTSUITE( ros::NodeHandle handle, std::string group_name ){
			
			// advertised services
			test90_service = handle.advertiseService( "/test90",  &TESTSUITE::test90, this);
			
			// published topics
			cmd_vel_publisher = handle.advertise<geometry_msgs::Twist>( group_name + "/drives/control/cmd_vel", 10);
		}
		
	private:
		void reset_cmd_vel( void ){
			cmd_vel.linear.x = 0.0;
			cmd_vel.linear.y = 0.0;
			cmd_vel.angular.z = 0.0;
		}
		
		bool test90( std_srvs::Empty::Request &req, std_srvs::Empty::Response &res ){
			geometry_msgs::Twist cmd_vel;
			double accel = 0.5;
			double deccel = 2.0;
			double vmax = 1.0;
			double start = 0;
			double end = M_PI_2;
			double diff = end - start;
			
			double Ts = 0.05;
			
			double t, t1, t2, t3, tend;
			double yaw_real, yaw1, yaw2, yaw3;
			
			t = 0;
			t1 = vmax / accel;
			t3 = vmax / deccel;
			
			yaw_real = 0;			
			yaw1 = accel * t1 * t1 / 2.0;
			yaw3 = deccel * t3 * t3 / 2.0;
			yaw2 = diff - yaw1 - yaw3;
			
			t2 = yaw2 / vmax;
			tend = t1 + t2 + t3;
			
			ROS_INFO("yaw1+2+3 = %f", (yaw1 + yaw2 + yaw3));
			ROS_INFO("t1 %f t2 %f t3 %f", t1, t2, t3);
			ROS_INFO("tend %f", tend);
			
			ros::Duration rate(Ts);
			while(t < tend) {
				if (t < t1){
					cmd_vel.angular.z = accel*t;
				} else if (t < t1+t2){
					cmd_vel.angular.z = vmax;
				} else if (t < t1+t2+t3){
					cmd_vel.angular.z = vmax-deccel*(t-(t1+t2));
				}
				
				cmd_vel_publisher.publish(cmd_vel);				
				t += Ts;
				yaw_real += cmd_vel.angular.z * Ts;
				rate.sleep();
			}
			
			cmd_vel.angular.z = 0.0;
			cmd_vel_publisher.publish(cmd_vel);
			
			ROS_INFO("yaw_real = %f", yaw_real);
			
			return true;
		}
		
	public:
		/**
		 * This function has to be called repeatedly to continuously transmit motion commands to the simulation.
		 */
		
}; // class base

int main( int argc, char** argv) {

	// init ros node
	ros::init(argc, argv, "ros_common_robin_testsuite");
	ros::NodeHandle handle;
	
	TESTSUITE tester( handle, "/mobrob_robin/base");
	
	ros::spin();	
	
}// main
