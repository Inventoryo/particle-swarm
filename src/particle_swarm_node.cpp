//
// Created by bb on 2020/3/9.
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

//#include "steering_functions/hc_cc_state_space/cc_dubins_state_space.hpp"

using namespace std;

//common
const double pi = 3.141592653;
const double dt = 1;
const double g = 9.8;

//map
const int resolution = 100;//?????50
const utility::point2D PosMin(0, 0);
const utility::point2D PosMax(100000, 100000);
int size_x = ceil((PosMax.y - PosMin.y) / resolution) + 1;
int size_y = ceil((PosMax.y - PosMin.y) / resolution) + 1;
vector<vector<utility::grid*>> global_map;
const int rato = 200;
const int sparse_resolution = resolution*rato;
int sparse_size_x = floor(size_x/rato);
int sparse_size_y = floor(size_y/rato);
vector<vector<utility::grid*>> sparse_map;

//target
const int target_num = 1;
vector<utility::TARGET*> target;

//uav
const int UAV_num = 2;
const utility::point2D VelMax(200, pi / 6);
const utility::point2D VelMin(150, -pi / 6);
double AccMax = 0.6;
double AccMin = -0.6;
const int search_R = 8000;
const double communication_R = 500000;//5000
const double tanTheta = sqrt(1.0 / 3);
const double forget_time = 10000;
vector<utility::UAV*> uav;

//particle
const int particle_num = 100;
const double particle_R = resolution+2*VelMax.x * dt;
const double tao = 1;
const double w1 = 1;
const double w2 = 1000;
const double w3 = 1;
const double weight = 1;
const double c1 = 2;
const double c2 = 0.001;
const double  w = 0.5;

//output
ofstream output_uav;
ofstream output_target;
ofstream output_traj_Point;
const string uav_path = "uav.txt";
const string target_path = "target.txt";
const string traj_Point_path = "traj_Point.txt";

//others
utility::nion union_;//Disjoint Set
vector<vector<double>> target_state(UAV_num, vector<double>(target_num, 0));
vector<bool> tracked(target_num,false);
//CC_Dubins_State_Space  cc_dubins_state_space_ptr_(1/(VelMin.x * VelMin.x / (g * tanTheta)),2*g*AccMax/(VelMin.x * VelMin.x*VelMin.x ) ,10);
//steer::State s1,s2;


double dist(utility::point2D& a, utility::point2D& b) {
    return sqrt((a.x - b.x)*(a.x - b.x) + (a.y - b.y)*(a.y - b.y));
}

double dubinsDistance(utility::State* state1,utility::State* state2){
    double minR = state1->velocity.x * state1->velocity.x / (g * tanTheta);//最小转弯半径
    double l = dist(state1->position, state2->position);
    double theTa1 = acos((state2->position.x - state1->position.x) / l);
    theTa1 = (state2->position.y - state1->position.y)>0?theTa1:2*pi- theTa1;
    double theTa2 = state1->velocity.y-theTa1;
    theTa2 = theTa2<pi?theTa2:theTa2-2*pi;
    double theTi,d;
    if(theTa2<-pi/2){
        theTa2+=pi;
         d = sqrt((l-minR*sin(theTa2))*(l-minR*sin(theTa2))+minR*cos(theTa2)*minR*cos(theTa2));
        theTi = asin(minR*cos(theTa2)/d) + asin(minR/d) + theTa1;
    }else if(theTa2<0){
        theTa2*=-1;
         d = sqrt((l-minR*sin(theTa2))*(l-minR*sin(theTa2))+minR*cos(theTa2)*minR*cos(theTa2));
        theTi = asin(minR/d) - asin(minR*cos(theTa2)/d) + theTa1;
    }else if(theTa2>pi/2){
        theTa2 = pi-theTa2;
         d = sqrt((l-minR*sin(theTa2))*(l-minR*sin(theTa2))+minR*cos(theTa2)*minR*cos(theTa2));
        theTi = -asin(minR/d) - asin(minR*cos(theTa2)/d) + theTa1;
    } else{
         d = sqrt((l-minR*sin(theTa2))*(l-minR*sin(theTa2))+minR*cos(theTa2)*minR*cos(theTa2));
        theTi = -asin(minR/d) + asin(minR*cos(theTa2)/d) + theTa1;
    }

    HybridAStar::DubinsPath* path = new HybridAStar::DubinsPath();
    double q1[3]={state1->position.x,state1->position.y,state1->velocity.y};//start point state
    double q2[3]={state2->position.x,state2->position.y,theTi};//end point state;
    HybridAStar::dubins_init(q1,q2,minR,path);
    return HybridAStar::dubins_path_length(path);

    /*
    s1.x = state1->position.x;
    s1.y = state1->position.y;
    s1.theta = state1->velocity.y;
    s1.kappa = 1/minR;
    s2.x = state2->position.x;
    s2.y = state2->position.y;
    s2.theta = theTi;
    s2.kappa = 1;
    return cc_dubins_state_space_ptr_.get_distance(s1,s2);
     */

}

void updateSubMap(utility::UAV * uav){
    //更新地图信息
    /*
    uav->coverd_area_cnt=0;
    int id =uav->id;
    for(int i =0;i<size_y;i++){
        for(int j =0;j<size_x;j++){
            global_map[i][j]->search_count[id]++;
            point2D grid_pose((j+0.5)*resolution,(i+0.5)*resolution);
            if (dist(grid_pose, uav->state->position) < search_R)
                global_map[i][j]->search_count[id] = 0;

            if(global_map[i][j]->search_count[id]>=forget_time)
                global_map[i][j]->search_count[id]=forget_time;
            else{
                uav->coverd_area_cnt++;
                //cout<<"global_map["<<i<<"]["<<j<<"]->search_count["<<id<<"]="<<global_map[i][j]->search_count[id]<<endl;
            }
        }
    }
    */
    for (int i = max(0,int((uav->state->position.y-search_R)/resolution)); i < min(size_y,int((uav->state->position.y+search_R)/resolution)); i++) {//j是在无人机搜索半径划分的正方形区域内所在的行数
        for (int j = max(0,int((uav->state->position.x-search_R)/resolution)); j < min(size_x,int((uav->state->position.x+search_R)/resolution)); j++) {//列
            utility::grid* temp_grid = global_map[i][j];//通过ID在全局地图中找出该栅格，并拿出
            utility::point2D grid_pose((j+0.5)*resolution,(i+0.5)*resolution);
            if (dist(grid_pose, uav->state->position) < search_R) //如果栅格在无人机探测半径内
                temp_grid->search_count[uav->id] = 0;
        }
    }

    return;
}

void updateUAVStatesInDubinsState(utility::UAV* uav){
    double minR = uav->state->velocity.x * uav->state->velocity.x / (g * tanTheta);//最小转弯半径
    double l = dist(uav->traj_Point, uav->state->position);
    double theTa1 = acos((uav->traj_Point.x - uav->state->position.x) / l);
    theTa1 = (uav->traj_Point.y - uav->state->position.y)>0?theTa1:2*pi- theTa1;
    double theTa2 = uav->state->velocity.y-theTa1;
    theTa2 = theTa2<pi?theTa2:theTa2-2*pi;
    double theTi,d;
    if(theTa2<-pi/2){
        theTa2+=pi;
        d = sqrt((l-minR*sin(theTa2))*(l-minR*sin(theTa2))+minR*cos(theTa2)*minR*cos(theTa2));
        theTi = asin(minR*cos(theTa2)/d) + asin(minR/d) + theTa1;
    }else if(theTa2<0){
        theTa2*=-1;
        d = sqrt((l-minR*sin(theTa2))*(l-minR*sin(theTa2))+minR*cos(theTa2)*minR*cos(theTa2));
        theTi = asin(minR/d) - asin(minR*cos(theTa2)/d) + theTa1;
    }else if(theTa2>pi/2){
        theTa2 = pi-theTa2;
        d = sqrt((l-minR*sin(theTa2))*(l-minR*sin(theTa2))+minR*cos(theTa2)*minR*cos(theTa2));
        theTi = -asin(minR/d) - asin(minR*cos(theTa2)/d) + theTa1;
    } else{
        d = sqrt((l-minR*sin(theTa2))*(l-minR*sin(theTa2))+minR*cos(theTa2)*minR*cos(theTa2));
        theTi = -asin(minR/d) + asin(minR*cos(theTa2)/d) + theTa1;
    }


    /*
    s1.x = uav->state->position.x;
    s1.y = uav->state->position.y;
    s1.theta = uav->state->velocity.y;
    s1.kappa = minR;
    s2.x = uav->traj_Point.x;
    s2.y = uav->traj_Point.y;
    s2.theta = theTi;
    s2.kappa = 1;
    double velocity = uav->state->velocity.x;
    double t=dt*velocity,total_lenth = cc_dubins_state_space_ptr_.get_distance(s1,s2);
    while(!uav->path_.empty())
        uav->path_.pop();
    vector<State> temp_path = cc_dubins_state_space_ptr_.get_path(s1,s2);

    while(t<=total_lenth){
        int index = t/10;
        vector<double> temp;
        temp.push_back(temp_path[index].x);
        temp.push_back(temp_path[index].y);
        temp.push_back(velocity);
        temp.push_back(temp_path[index].theta);
        uav->path_.push(temp);
        t+=velocity*dt;
    }
    */

    HybridAStar::DubinsPath* path = new HybridAStar::DubinsPath();
    double q0[3]={uav->state->position.x,uav->state->position.y,uav->state->velocity.y};//start point state
    double q1[3]={uav->traj_Point.x,uav->traj_Point.y,theTi};//end point state;
    HybridAStar::dubins_init(q0,q1,minR,path);
    double velocity = uav->state->velocity.x;
    double t=dt*velocity,total_lenth = HybridAStar::dubins_path_length(path);

    while(!uav->path_.empty())
        uav->path_.pop();

    while(t<=total_lenth){
        double temp_point[3];
        HybridAStar::dubins_path_sample(path,t,temp_point);
        vector<double> temp;
        temp.push_back(temp_point[0]);
        temp.push_back(temp_point[1]);
        temp.push_back(velocity);
        temp.push_back(temp_point[2]);
        uav->path_.push(temp);
        if(t>path->param[0] && t<path->param[1]&&velocity+AccMax*dt<VelMax.x){
            t+=velocity*dt+0.5*AccMax*dt*dt;
            velocity=min(VelMax.x,velocity+AccMax*dt);
        }else{
            t+=velocity*dt;
        }
    }

    return;
}

void spreadParticles(utility::UAV * uav) {
    srand((unsigned)time(NULL));
    for (int j = 0; j < particle_num; j++) {
        utility::particle* temp_particle =uav->swarm[j] ;
        double length = search_R + resolution + particle_R * (rand() % 1000) / 1000.0;
        double theta = 2 * pi*((rand() % 1000) / 1000.0);
        temp_particle->state->position.x =  length * cos(theta) + uav->state->position.x;
        temp_particle->state->position.y =  length * sin(theta) + uav->state->position.y;
        temp_particle->state->velocity.x = VelMin.x + (VelMax.x - VelMin.x) * ((rand() % 1000) / 1000.0);//线速度
        temp_particle->state->velocity.y = 2 * pi*((rand() % 1000) / 1000.0);//角度
        temp_particle->state->accelerate.x = 0;//线加速度
        temp_particle->state->accelerate.y = 0;//角速度
        //边界判断
        if (temp_particle->state->position.x < PosMin.x) {
            temp_particle->state->position.x = PosMin.x;
            temp_particle->state->velocity.y = pi*(0.5 -(rand() % 1000) / 1000.0);//角度
        }
        else if (temp_particle->state->position.x > PosMax.x) {
            temp_particle->state->position.x = PosMax.x;
            temp_particle->state->velocity.y = pi * (0.5 + (rand() % 1000) / 1000.0);//角度
        }

        if (temp_particle->state->position.y < PosMin.y) {
            temp_particle->state->position.y = PosMin.y;
            temp_particle->state->velocity.y = pi * ( (rand() % 1000) / 1000.0);//角度
        }
        else if (temp_particle->state->position.y > PosMax.y) {
            temp_particle->state->position.y = PosMax.y;
            temp_particle->state->velocity.y = pi * (1 + (rand() % 1000) / 1000.0);//角度
        }
        temp_particle->Pbest_position = temp_particle->state->position;//局部最优
        temp_particle->fitness = 0;
    }
    uav->Gbest_position = uav->state->position;
    uav->Gbest_fitness = 0;
    return;
}

void updateParticleStates(utility::UAV * uav){//更新粒子
    uav->Gbest_fitness = 0;
    int best_x=0,best_y=0;
    double best_grid_fit =0,temp_grid_fit=0;
    utility::point2D best_grid(0,0);
    for(int i =0;i<sparse_size_y;i++){
        for(int j=0;j<sparse_size_x;j++){
            utility::State* temp = new utility::State();
            temp->position.x=(i+0.5)*sparse_resolution;
            temp->position.y=(j+0.5)*sparse_resolution;
            temp_grid_fit = 1/max(dubinsDistance(uav->state,temp),dist(uav->state->position,temp->position)) + sparse_map[i][j]->search_count[uav->id];
            if(temp_grid_fit>best_grid_fit){
                best_grid_fit = temp_grid_fit;
                best_grid.x=(j+0.5)*sparse_resolution;
                best_grid.y=(i+0.5)*sparse_resolution;
            }
        }
    }
    for(int i = 0;i<50;i++){
        uav->last_Gbest_position = uav->Gbest_position;
        for (int j = 0; j < particle_num; j++) {
            utility::particle & temp_part = *(uav->swarm[j]);
            //速度更新
            double r1 = (double)rand() / RAND_MAX;
            double r2 = (double)rand() / RAND_MAX;
            double line_X = temp_part.state->velocity.x * cos(temp_part.state->velocity.y);//速度在XY向的分量
            double line_Y = temp_part.state->velocity.x * sin(temp_part.state->velocity.y);
            line_X = weight * line_X + c1 * r1 *(temp_part.Pbest_position.x - temp_part.state->position.x) + c2 * r2 *(uav->Gbest_position.x - temp_part.state->position.x);
            line_Y = weight * line_Y + c1 * r1 *(temp_part.Pbest_position.y - temp_part.state->position.y) + c2 * r2 *(uav->Gbest_position.y - temp_part.state->position.y);
            temp_part.state->velocity.x = sqrt(line_X * line_X + line_Y * line_Y);//粒子的速度更新公式
            if (line_Y >= 0)
                temp_part.state->velocity.y = acos(line_X / temp_part.state->velocity.x);
            else
                temp_part.state->velocity.y = 2*pi - acos(line_X / temp_part.state->velocity.x);

            //速度限制
            if (temp_part.state->velocity.x > VelMax.x)
                temp_part.state->velocity.x = VelMax.x;
            else if(temp_part.state->velocity.x < VelMin.x)
                temp_part.state->velocity.x = VelMin.x;
            //粒子的位置更新公式
            temp_part.state->position.x += temp_part.state->velocity.x * dt*cos(temp_part.state->velocity.y);
            temp_part.state->position.y += temp_part.state->velocity.x * dt*sin(temp_part.state->velocity.y);
            //边界判断
            if (temp_part.state->position.x < PosMin.x) {
                temp_part.state->position.x = PosMin.x;
                temp_part.state->velocity.y = pi * (0.5 - (rand() % 1000) / 1000.0);//角度
            }
            else if (temp_part.state->position.x > PosMax.x) {
                temp_part.state->position.x = PosMax.x;
                temp_part.state->velocity.y = pi * (0.5 + (rand() % 1000) / 1000.0);//角度
            }
            if (temp_part.state->position.y < PosMin.y) {
                temp_part.state->position.y = PosMin.y;
                temp_part.state->velocity.y = pi * ((rand() % 1000) / 1000.0);//角度
            }
            else if (temp_part.state->position.y > PosMax.y) {
                temp_part.state->position.y = PosMax.y;
                temp_part.state->velocity.y = pi * (1 + (rand() % 1000) / 1000.0);//角度
            }

            //概率更新
            int x = temp_part.state->position.x/resolution;
            int y = temp_part.state->position.y/resolution;
            int t = global_map[y][x]->search_count[uav->id];
            //cout<<"size_x*size_y-uav->coverd_area_cnt="<<size_x*size_y-uav->coverd_area_cnt<<endl;
            //cout<<"1 - exp(-tao *t )="<<1 - exp(-tao *t )<<endl;
            for (int k = 0; k < target_num; k++){
                temp_part.p[k] = 1000000*(1 - exp(-tao *t )) / (size_x*size_y-uav->coverd_area_cnt);
            }

            //适应值更新
            //fitness1
            double temp_fitness_1 = 0, temp_fitness_2 = 0, temp_fitness_3 = 0;
            for (int k = 0; k < target_num; k++) {
                double Dik =   max(dubinsDistance(uav->state, temp_part.state),dist(uav->state->position, temp_part.state->position));//无人机k与粒子i的距离
                double Pij = temp_part.p[k];
                double dik = exp(-(Dik/search_R));
                int Tj = uav->Tj[k];//目标k是否未被与之更匹配的无人机跟踪
                //temp_fitness_1 += Tj * (Pij  + Dik +Cjk);
                temp_fitness_1 += (Pij * dik)*  Tj ;
                //cout<<"1"<<endl;
            }
            //fitness2
            if(t>=forget_time)
                temp_fitness_2 = 10000*w * 1.0 / (size_x*size_y - uav->coverd_area_cnt);
            else
                temp_fitness_2 = 10000*(1 - w)*t / (forget_time * (size_x*size_y - uav->coverd_area_cnt));
            //cout << "temp_fitness_1 = " << temp_fitness_1 <<";temp_fitness_2 = "<<temp_fitness_2<<endl;

            //fitness3
            double angle1 = acos((best_grid.x-uav->state->position.x)/dist(uav->state->position,best_grid));
            angle1 = (best_grid.y-uav->state->position.y)>0?angle1:2*pi-angle1;//0~2*pi
            double angle2 = acos((temp_part.state->position.x-uav->state->position.x)/dist(uav->state->position,temp_part.state->position));
            angle2 = (temp_part.state->position.y-uav->state->position.y)>0?angle2:2*pi-angle2;//0~2*pi
            temp_fitness_3 = 1/abs(angle2-angle1);

            double temp_fitness = w1 * temp_fitness_1 + w2 * temp_fitness_2 + w3* temp_fitness_3;//适应值
            //cout << "temp_fitness = " << temp_fitness <<endl;
            if (temp_fitness >= temp_part.fitness) {//局部最优位置判断
                temp_part.fitness = temp_fitness;
                temp_part.Pbest_position = temp_part.state->position;
            }
            //cout<<j<<endl;

            if (temp_fitness >= uav->Gbest_fitness) {//全局最优位置判断
                //cout << "更新temp_fitness" << temp_fitness << " current particle" << j << endl;
                uav->Gbest_fitness = temp_fitness;
                uav->Gbest_position = temp_part.state->position;
            }
        }
        //if (uav->last_Gbest_position.x() == uav->Gbest_position.x() && uav->last_Gbest_position.y() == uav->Gbest_position.y())
        //    break;
    }
    //cout<<"uav->Gbest_fitness"<<uav->Gbest_fitness<<endl;
    uav->traj_Point = uav->Gbest_position;
    return;
}

bool init() {
    //target
    for(int i =0;i<target_num;i++){
        utility::TARGET* temp_target = new utility::TARGET();
        temp_target->state->position.x = 30000 +i*20000;//PosMin.x + (PosMax.x - PosMin.x) * (double)rand() / RAND_MAX;10+i*100
        temp_target->state->position.y = 50000+i*20000;//PosMin.y + (PosMax.y - PosMin.y) * (double)rand() / RAND_MAX;30+i*200
        temp_target->state->velocity.x = 40;//线速度
        temp_target->state->velocity.y = 2 * pi*rand() / RAND_MAX;//角度
        temp_target->state->accelerate.x = 0;//线加速度
        temp_target->state->accelerate.y = 0;//角速度
        target.push_back(temp_target);
    }
    //map
    for (int i = 0; i < size_y ; i++) {//初始化地图，i=0代表第一行，j=0代表第一列
        vector<utility::grid*> temp;
        for (int j = 0; j < size_x ; j++) {
            utility::grid* tempGridPtr = new utility::grid();
            for(int k=0;k<UAV_num;k++)
                tempGridPtr->search_count.push_back(forget_time);
            temp.push_back(tempGridPtr);
        }
        global_map.push_back(temp);
    }
    //sparse_map
    for (int i = 0; i < sparse_size_y ; i++) {//初始化地图，i=0代表第一行，j=0代表第一列
        vector<utility::grid*> temp;
        for (int j = 0; j < sparse_size_x ; j++) {
            utility::grid* tempGridPtr = new utility::grid();
            for(int k=0;k<UAV_num;k++)
                tempGridPtr->search_count.push_back(forget_time);
            temp.push_back(tempGridPtr);
        }
        sparse_map.push_back(temp);
    }
    //UAV
    for (int i = 0; i < UAV_num; i++) {//初始化无人机
        utility::UAV* uav_temp = new utility::UAV();
        srand((int(time(NULL)) + i));
        uav_temp->id =i;
        uav_temp->state->position.x = 1000;//位置
        uav_temp->state->position.y = 10000*i;
        uav_temp->state->velocity.x = VelMin.x + (VelMax.x - VelMin.x) * (rand() / double(RAND_MAX));//线速度
        uav_temp->state->velocity.y = 0;//角度
        uav_temp->state->accelerate.x = 0;//线加速度
        uav_temp->state->accelerate.y = 0;//角速度
        uav_temp->search_r = search_R;
        uav_temp->particle_r = particle_R;
        uav_temp->coverd_area_cnt =5000;
        for (int j = 0; j < target_num; j++) {
            utility::point2D target_pose(-1,-1);
            uav_temp->target_position.push_back(make_pair(target_pose,forget_time));
            uav_temp->Tj[j] = 1;
        }
        updateSubMap(uav_temp);
        for (int j = 0; j < particle_num; j++) {
            utility::particle* temp = new utility::particle();
            uav_temp->swarm.push_back(temp);
        }
        cout<<"1"<<endl;
        spreadParticles(uav_temp);
        cout<<"2"<<endl;
        updateParticleStates(uav_temp);
        cout<<"3"<<endl;
        updateUAVStatesInDubinsState(uav_temp);
        cout<<"4"<<endl;
        while(uav_temp->path_.size()>1){//only keep one state
            uav_temp->path_.pop();
        }

        uav_temp->track_target_num = -1;
        uav.push_back(uav_temp);
        union_.parent.push_back(i);
    }
    return 1;
};

void updateTargetStates() {
    //目标状态更新，这里是等速运动
    output_target.open(target_path.c_str(), ios::app | ios::binary);
    for (int target_tag = 0; target_tag < target_num; target_tag++) {
        double line_X = target[target_tag]->state->velocity.x * cos(target[target_tag]->state->velocity.y);//速度在XY的分量
        double line_Y = target[target_tag]->state->velocity.x * sin(target[target_tag]->state->velocity.y);
        target[target_tag]->state->position.x += line_X * dt;//位置更新
        target[target_tag]->state->position.y += line_Y * dt;

        //边界判断
        if (target[target_tag]->state->position.x < PosMin.x) {
            target[target_tag]->state->position.x = PosMin.x;
            target[target_tag]->state->velocity.y = pi * (0.5 - (rand() % 1000) / 1000.0);//角度
        }
        else if (target[target_tag]->state->position.x > PosMax.x) {
            target[target_tag]->state->position.x = PosMax.x;
            target[target_tag]->state->velocity.y = pi * (0.5 + (rand() % 1000) / 1000.0);//角度
        }

        if (target[target_tag]->state->position.y < PosMin.y) {
            target[target_tag]->state->position.y = PosMin.y;
            target[target_tag]->state->velocity.y = pi * ((rand() % 1000) / 1000.0);//角度
        }
        else if (target[target_tag]->state->position.y > PosMax.y) {
            target[target_tag]->state->position.y = PosMax.y;
            target[target_tag]->state->velocity.y = pi * (1 + (rand() % 1000) / 1000.0);//角度
        }
        output_target << target[target_tag]->state->position.x << " " << target[target_tag]->state->position.y << " ";
    }
    output_target << endl;
    output_target.close();
    return;
};

void updateUAVStates1(){
    output_uav.open(uav_path.c_str(), ios::app | ios::binary);
    output_traj_Point.open(traj_Point_path.c_str(), ios::app | ios::binary);
    for(int i =0;i<UAV_num;i++){
        if(uav[i]->path_.empty()){
            uav[i]->state->velocity.x=max(uav[i]->state->velocity.x+AccMin*dt,VelMin.x);
            uav[i]->state->position.x += uav[i]->state->velocity.x*cos(uav[i]->state->velocity.y)*dt;
            uav[i]->state->position.y += uav[i]->state->velocity.x*sin(uav[i]->state->velocity.y)*dt;
            output_uav << uav[i]->state->position.x << " " << uav[i]->state->position.y << " " << uav[i]->track_target_num << " ";
            output_traj_Point << uav[i]->traj_Point.x << " " << uav[i]->traj_Point.y << " ";
            cout<<"error"<<endl;
            continue;
        }
        uav[i]->state->position.x = uav[i]->path_.front()[0];
        uav[i]->state->position.y = uav[i]->path_.front()[1];
        uav[i]->state->velocity.x = uav[i]->path_.front()[2];
        uav[i]->state->velocity.y = uav[i]->path_.front()[3];
        uav[i]->path_.pop();
        output_uav << uav[i]->state->position.x << " " << uav[i]->state->position.y << " " << uav[i]->track_target_num << " ";
        output_traj_Point << uav[i]->traj_Point.x << " " << uav[i]->traj_Point.y << " ";
        updateSubMap(uav[i]);
        for(int j=0;j<target_num;j++){//target inf0
            if(dist(target[j]->state->position,uav[i]->state->position)<search_R) {
                uav[i]->target_position[j].first = target[j]->state->position;
                uav[i]->target_position[j].second =0;
            } else
                uav[i]->target_position[j].second++;//search_cunt++
        }
    }
    output_uav << endl;
    output_uav.close();
    output_traj_Point << endl;
    output_traj_Point.close();
    return;
}

void informationShare() {
    for (int i = 0; i < UAV_num; i++) {//reset
        union_.parent[i]=i;
        uav[i]->coverd_area_cnt=0;
    }
    for (int i = 0; i < UAV_num; i++) {
        for (int j = i; j < UAV_num; j++) {
            if(dist(uav[i]->state->position, uav[j]->state->position) <= communication_R){//union
                union_.join(i,j);
            }
            for(int k=0;k<target_num;k++){//exchange target info
                if(uav[i]->target_position[k].second>uav[j]->target_position[k].second)
                    uav[i]->target_position[k] = uav[j]->target_position[k];
                else
                    uav[j]->target_position[k] = uav[i]->target_position[k];
            }
        }
    }

    union_.setup_barrel();//setup barrel
    int sum[9] ={0};
    for(int i =0;i<size_y;i++){
        for (int j = 0; j <size_y ; j++) {
            for(int k = 0;k<union_.barrel.size();k++){
                int min_count = forget_time;
                for(int l =0;l<union_.barrel[k].size();l++){
                    min_count =std::min(min_count,++global_map[i][j]->search_count[union_.barrel[k][l]]);//find the lowest search time
                }
                for(int l =0;l<union_.barrel[k].size();l++){
                    global_map[i][j]->search_count[union_.barrel[k][l]] = min_count;//set the search time
                    if(min_count<forget_time)
                        uav[union_.barrel[k][l]]->coverd_area_cnt++;//count the coverd area
                    //sparse_map update
                    if(i==0||j==0)
                        continue;
                    if(i%rato==0&&j%rato==0){
                        sparse_map[j/rato-1][i/rato-1]->search_count[union_.barrel[k][l]] = sum[union_.barrel[k][l]]/40000;
                        sum[union_.barrel[k][l]]=0;
                    } else{
                        sum[union_.barrel[k][l]]+=min_count;
                    }

                }
            }

        }
    }


    return;
}

void updateMission() {
    for (int i = 0; i < UAV_num; i++) {
        for (int j = 0; j < target_num; j++) {//如果无人机i发现的目标j所在的栅格的搜索次数大于1000次
            //cout << "uav[i]->target_position[j].second = " << uav[i]->target_position[j].second  << endl;
            if (uav[i]->target_position[j].second >= forget_time)
                target_state[i][j] = 0;//不让无人机i去追目标j
            else
                target_state[i][j] = 1000.0 / dist(uav[i]->state->position, target[j]->state->position);
        }
    }
    set<int> row, col;
    for (int l = 0; l < target_num; l++) {
        double maxTargetState = 0;
        int uav_cnt, target_cnt;
        for (int i = 0; i < UAV_num; i++) {//寻找满足条件的最大值，并行与列
            for (int j = 0; j < target_num; j++) {
                if (row.find(i) == row.end() && col.find(j) == col.end() && target_state[i][j] > maxTargetState) {
                    maxTargetState = target_state[i][j];
                    uav_cnt = i;
                    target_cnt = j;
                }
            }
        }
        if (maxTargetState > 0.001) {//当能够找到有效的任务分配时，执行任务分配
            cout << "uav_cnt = " << uav_cnt << " target_cnt = " << target_cnt << endl;
            if (dist(uav[uav_cnt]->state->position, target[target_cnt]->state->position) < 500) {
                tracked[target_cnt] = true;
                uav[uav_cnt]->state->velocity.x = target[target_cnt]->state->velocity.x;
                cout << "target " << target_cnt << " is successfully tracked by uav " << uav_cnt << endl;
            }

            for (int i = 0; i < UAV_num; i++) {
                if (i == uav_cnt) {
                    uav[i]->Tj[target_cnt] = 1;//分配该无人机追踪该目标
                }
                else
                    uav[i]->Tj[target_cnt] = 0;//否则分配该无人机不追踪该目标
            }
            row.insert(uav_cnt); col.insert(target_cnt);//将行与列放入set中，下次循环跳过该行与列
        }
    }
    for (int i = 0; i < UAV_num; i++) {
        int j = 0;
        for (; j < target_num; j++) {
            //cout << "uav[" << i << "].TJ[" << j << "]= " << uav[i].Tj[j] << " "<< uav[i].covered_target_id[j].second<< endl;
            if (uav[i]->Tj[j] == 1 && uav[i]->target_position[j].second < forget_time) {//分配该无人机追踪该目标
                uav[i]->track_target_num = j;
                uav[i]->traj_Point = target[j]->state->position;
                uav[i]->state->velocity.x = 2*target[j]->state->velocity.x;
                updateUAVStatesInDubinsState(uav[i]);
                cout << "uav " << i << " is tracking target " << j << endl;
                break;
            }
        }
        if (j >= target_num) {
            uav[i]->track_target_num = -1;
            cout << "uav_cnt = " << i << " is free" << endl;
        }//该无人机尚未分配任何目标
    }
}

int main(int argc, char** argv)
{
    //ros::init(argc, argv, "particle_swarm_node");
    //ros::NodeHandle nh("~");

    //初始化
    init();
    int cunt = 0;
    while (cunt < 10000) {
        //如果无人机到达全局最优，则重新洒粒子
        for (int i = 0; i < UAV_num; i++) {
            if (uav[i]->path_.empty()) {
                cout<<endl<<endl<<endl<<"replanning uav"<<i<<endl<<endl<<endl;
                spreadParticles(uav[i]);
                updateParticleStates(uav[i]);
                updateUAVStatesInDubinsState(uav[i]);
            }
        }

        //更新目标位置
        updateTargetStates();

        //更新无人机状态
        updateUAVStates1();

        //机间通信
        if(UAV_num>1)
            informationShare();

        //任务分配
        updateMission();

        cout << "cunt" << cunt << endl;

        //判断终止
        int i =0;
        for(;i<target_num;i++){
            if(!tracked[i])
                break;
        }
        if (i>=target_num)
            break;

        cunt++;
    }
    //std::cout << "target_state:" << endl << target_state << endl;
    std::cout << "All target has been tracked!The total number of step is " << cunt << std::endl;
    char aa = getchar();
    return 0;
}