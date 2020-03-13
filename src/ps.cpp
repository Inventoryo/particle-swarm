//
// Created by bb on 2020/3/13.
//

#include <ros/ros.h>
#include <ros/console.h>
#include <sensor_msgs/PointCloud2.h>

#include "utility.h"
#include <cmath>
#include<ctime>
#include<cstdlib>
#include <iostream>
#include <fstream>
#include <set>
#include "dubins.h"
#include <ompl/base/spaces/ReedsSheppStateSpace.h>
#include <ompl/base/spaces/DubinsStateSpace.h>
#include <ompl/base/spaces/SE2StateSpace.h>
#include <ompl/base/State.h>


using namespace std;

const double pi = 3.141592653;
const int resolution = 500;//?????50
const int particle_num = 50;
const int UAV_num = 9;
const int target_num = 3;

const Eigen::Vector2d PosMin(0, 0);
const Eigen::Vector2d PosMax(100000, 100000);

const Eigen::Vector2d VelMax(200, pi / 6);
const Eigen::Vector2d VelMin(150, -pi / 6);
const double AccMax=0.6;

const int search_R = 8000;
const double dt = 1;
const double particle_R = VelMax.x * dt;
const double communication_R = 500000;//5000
const double g = 9.8;
const double tanTheta = sqrt(1.0 / 3);
const double forget_time = 1000;

const double tao = 1;

//粒子参数
const double w1 = 1;
const double w2 = 200;
const double weight = 1;
const double c1 = 2;
const double c2 = 2;
const double  w = 0.5;

ofstream output_uav;
ofstream output_target;
ofstream output_traj_Point;
const string uav_path = "uav.txt";
const string target_path = "target.txt";
const string traj_Point_path = "traj_Point.txt";

vector<UAV*> uav;
vector<TARGET*> target;
map<int, grid*> global_map;
int size_x = (PosMax.x - PosMin.x) / resolution + 1;
int size_y = (PosMax.y - PosMin.y) / resolution + 1;

vector<vector<double>> target_state(UAV_num, vector<double>(target_num, 0));
bool tracked[3] = { false };


double dist(Eigen::Vector2d& a, Eigen::Vector2d& b) {
    return sqrt((a.x - b.x)*(a.x - b.x) + (a.y - b.y)*(a.y - b.y));
}

void updateUAVStatesInDubinsState(UAV* uav){

    double theTi = acos((uav->traj_Point.x - uav->position.x) / dist(uav->traj_Point, uav->position));
    HybridAStar::DubinsPath* path = new HybridAStar::DubinsPath();
    double minR = uav->velocity.x * uav->velocity.x / (g * tanTheta);//最小转弯半径
    double q0[3]={uav->position.x,uav->position.y,uav->velocity.y};//start point state
    double q1[3]={uav->traj_Point.x,uav->traj_Point.y,theTi};//end point state;
    HybridAStar::dubins_init(q0,q1,1/minR,path);
    double velocity = uav->velocity.x;
    double t=dt*velocity,total_lenth = HybridAStar::dubins_path_length(path);

    while(!uav->path_.empty())
        uav->path_.pop();

    while(t<=total_lenth){
        double temp_point[3];
        HybridAStar::dubins_path_sample(path,t,temp_point);
        double temp[4];
        temp[0]=temp_point[0];
        temp[1]=temp_point[1];
        temp[2]=temp_point[2];
        temp[3]=velocity;
        uav->path_.push(temp);
        if(t>path->param[0] && t<path->param[1]&&velocity+AccMax*dt<VelMax.x){
            t+=velocity*dt+0.5*AccMax*dt*dt;
            velocity=min(VelMax.x,velocity+AccMax*dt);
        }else{
            t+=velocity*dt;
        }
    }

}

void spreadParticles(UAV * uav) {

    srand((unsigned)time(NULL));
    uav.swarm.clear();
    for (int j = 0; j < particle_num; j++) {
        particle temp_particle;
        double length = search_R + VelMax.x * dt * (rand() % 1000) / 1000.0;
        double theta = 2 * pi*((rand() % 1000) / 1000.0);
        temp_particle.position.x =  length * cos(theta) + uav.position.x;
        temp_particle.position.y =  length * sin(theta) + uav.position.y;
        temp_particle.velocity.x = VelMax.x + (VelMax.x - VelMin.x) * ((rand() % 1000) / 1000.0);//线速度
        temp_particle.velocity.y = 2 * pi*((rand() % 1000) / 1000.0);//角度
        temp_particle.accelerate.x = 0;//线加速度
        temp_particle.accelerate.y = 0;//角速度
        temp_particle.Pbest_position = temp_particle.position;//局部最优
        temp_particle.fitness = 0;

        //边界判断
        if (temp_particle.position.x < PosMin.x) {
            temp_particle.position.x = PosMin.x;
            temp_particle.velocity.y = pi*(0.5 -(rand() % 1000) / 1000.0);//角度
        }
        else if (temp_particle.position.x > PosMax.x) {
            temp_particle.position.x = PosMax.x;
            temp_particle.velocity.y = pi * (0.5 + (rand() % 1000) / 1000.0);//角度
        }

        if (temp_particle.position.y < PosMin.y) {
            temp_particle.position.y = PosMin.y;
            temp_particle.velocity.y = pi * ( (rand() % 1000) / 1000.0);//角度
        }
        else if (temp_particle.position.y > PosMax.y) {
            temp_particle.position.y = PosMax.y;
            temp_particle.velocity.y = pi * (1 + (rand() % 1000) / 1000.0);//角度
        }

        uav.swarm.push_back(temp_particle);
    }
    uav.Gbest_position = uav.swarm.front().Pbest_position;
    uav.Gbest_fitness = 0;
    uav.velocity.x = VelMin.x;
    return;
}

