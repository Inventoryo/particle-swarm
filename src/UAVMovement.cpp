//
// Created by bb on 2020/3/12.
//

#include "utility.h"
#include <vector>

#include <ros/ros.h>
#include <ros/console.h>
#include <sensor_msgs/PointCloud2.h>

using namespace std;
const double pi = 3.141592653;

vector<TARGET> target;
void init(){
    TARGET temp_target;
    temp_target.position.x = 20000;//PosMin.x + (PosMax.x - PosMin.x) * (double)rand() / RAND_MAX;10+i*100
    temp_target.position.y = 50000;//PosMin.y + (PosMax.y - PosMin.y) * (double)rand() / RAND_MAX;30+i*200
    temp_target.velocity.x = 150;//线速度
    temp_target.velocity.y = 2 * pi*rand() / RAND_MAX;//角度
    temp_target.accelerate.x = 0;//线加速度
    temp_target.accelerate.y = 0;//角速度
    target.push_back(temp_target);

    temp_target.position.x = 30000;//PosMin.x + (PosMax.x - PosMin.x) * (double)rand() / RAND_MAX;10+i*100
    temp_target.position.y = 60000;//PosMin.y + (PosMax.y - PosMin.y) * (double)rand() / RAND_MAX;30+i*200
    temp_target.velocity.x = 150;//线速度
    temp_target.velocity.y = 2 * pi*rand() / RAND_MAX;//角度
    temp_target.accelerate.x = 0;//线加速度
    temp_target.accelerate.y = 0;//角速度
    target.push_back(temp_target);

    temp_target.position.x = 50000;//PosMin.x + (PosMax.x - PosMin.x) * (double)rand() / RAND_MAX;10+i*100
    temp_target.position.y = 70000;//PosMin.y + (PosMax.y - PosMin.y) * (double)rand() / RAND_MAX;30+i*200
    temp_target.velocity.x = 150;//线速度
    temp_target.velocity.y = 2 * pi*rand() / RAND_MAX;//角度
    temp_target.accelerate.x = 0;//线加速度
    temp_target.accelerate.y = 0;//角速度*/
    target.push_back(temp_target);
}

int main(int argc, char** argv){
    ros::init(argc, argv, "UAV Movement node");
    ros::NodeHandle nh("~");

    init();



}