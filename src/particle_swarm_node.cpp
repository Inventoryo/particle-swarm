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
double AccMax = 0.6;
double AccMin = -0.6;
const int search_R = 8000;
const double dt = 1;
const double particle_R = VelMax.x() * dt;
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
vector<vector<grid*>> global_map;
nion union_;//Disjoint Set
int size_x = ceil((PosMax.y() - PosMin.y()) / resolution) + 1;
int size_y = ceil((PosMax.y() - PosMin.y()) / resolution) + 1;

vector<vector<double>> target_state(UAV_num, vector<double>(target_num, 0));
bool tracked[3] = { false };

double dist(Eigen::Vector2d& a, Eigen::Vector2d& b) {
    return sqrt((a.x() - b.x())*(a.x() - b.x()) + (a.y() - b.y())*(a.y() - b.y()));
}

void updateUAVStatesInDubinsState(UAV* uav){

    double theTi = acos((uav->traj_Point.x() - uav->position.x()) / dist(uav->traj_Point, uav->position));
    HybridAStar::DubinsPath* path = new HybridAStar::DubinsPath();
    double minR = uav->velocity.x() * uav->velocity.x() / (g * tanTheta);//最小转弯半径
    double q0[3]={uav->position.x(),uav->position.y(),uav->velocity.y()};//start point state
    double q1[3]={uav->traj_Point.x(),uav->traj_Point.y(),theTi};//end point state;
    HybridAStar::dubins_init(q0,q1,minR,path);
    double velocity = uav->velocity.x();
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
        if(t>path->param[0] && t<path->param[1]&&velocity+AccMax*dt<VelMax.x()){
            t+=velocity*dt+0.5*AccMax*dt*dt;
            velocity=min(VelMax.x(),velocity+AccMax*dt);
        }else{
            t+=velocity*dt;
        }
    }
}

void spreadParticles(UAV * uav) {
    srand((unsigned)time(NULL));
    for (int j = 0; j < particle_num; j++) {
        particle* temp_particle =uav->swarm[j] ;
        double length = search_R + VelMax.x() * dt * (rand() % 1000) / 1000.0;
        double theta = 2 * pi*((rand() % 1000) / 1000.0);
        temp_particle->position.x() =  length * cos(theta) + uav->position.x();
        temp_particle->position.y() =  length * sin(theta) + uav->position.y();
        temp_particle->velocity.x() = VelMax.x() + (VelMax.x() - VelMin.x()) * ((rand() % 1000) / 1000.0);//线速度
        temp_particle->velocity.y() = 2 * pi*((rand() % 1000) / 1000.0);//角度
        temp_particle->accelerate.x() = 0;//线加速度
        temp_particle->accelerate.y() = 0;//角速度
        temp_particle->Pbest_position = temp_particle->position;//局部最优
        temp_particle->fitness = 0;

        //边界判断
        if (temp_particle->position.x() < PosMin.x()) {
            temp_particle->position.x() = PosMin.x();
            temp_particle->velocity.y() = pi*(0.5 -(rand() % 1000) / 1000.0);//角度
        }
        else if (temp_particle->position.x() > PosMax.x()) {
            temp_particle->position.x() = PosMax.x();
            temp_particle->velocity.y() = pi * (0.5 + (rand() % 1000) / 1000.0);//角度
        }

        if (temp_particle->position.y() < PosMin.y()) {
            temp_particle->position.y() = PosMin.y();
            temp_particle->velocity.y() = pi * ( (rand() % 1000) / 1000.0);//角度
        }
        else if (temp_particle->position.y() > PosMax.y()) {
            temp_particle->position.y() = PosMax.y();
            temp_particle->velocity.y() = pi * (1 + (rand() % 1000) / 1000.0);//角度
        }
    }
    uav->Gbest_position = uav->swarm.front()->Pbest_position;
    uav->Gbest_fitness = 0;
    return;
}

void updateSubMap(UAV * uav){
    //更新地图信息
    for (int i = max(0,int((uav->position.y()-search_R)/resolution)); i < min(size_y,int((uav->position.y()+search_R)/resolution)); i++) {//j是在无人机搜索半径划分的正方形区域内所在的行数
        for (int j = max(0,int((uav->position.x()-search_R)/resolution)); j < min(size_x,int((uav->position.x()+search_R)/resolution)); j++) {//列
            grid* temp_grid = global_map[i][j];//通过ID在全局地图中找出该栅格，并拿出
            Eigen::Vector2d grid_pose((j+0.5)*resolution,(i+0.5)*resolution);
            if (dist(grid_pose, uav->position) < search_R) //如果栅格在无人机探测半径内
                temp_grid->search_count[uav->id] = 0;
        }
    }
}

void updateParticleStates(UAV * uav){//更新粒子
    //uav[i].Gbest_fitness = 0;
    if (uav->track_target_num != -1 && tracked[uav->track_target_num])//如果无人机已经追踪到目标，则粒子群更新停止
        return;
    //if (uav[i].isConvergenced)//如果粒子群已收敛，则不更新
    uav->last_Gbest_position = uav->Gbest_position;
    for (int j = 0; j < particle_num; j++) {
        particle & temp_part = *(uav->swarm[j]);
        //速度更新
        double r1 = (double)rand() / RAND_MAX;
        double r2 = (double)rand() / RAND_MAX;
        double line_X = temp_part.velocity.x() * cos(temp_part.velocity.y());//速度在XY向的分量
        double line_Y = temp_part.velocity.x() * sin(temp_part.velocity.y());
        line_X = weight * line_X + c1 * r1 *(temp_part.Pbest_position.x() - temp_part.position.x()) + c2 * r2 *(uav->Gbest_position.x() - temp_part.position.x());
        line_Y = weight * line_Y + c1 * r1 *(temp_part.Pbest_position.y() - temp_part.position.y()) + c2 * r2 *(uav->Gbest_position.y() - temp_part.position.y());
        temp_part.velocity.x() = sqrt(line_X * line_X + line_Y * line_Y);//粒子的速度更新公式
        if (line_Y >= 0)
            temp_part.velocity.y() = acos(line_X / temp_part.velocity.x());
        else
            temp_part.velocity.y() = 2*pi - acos(line_X / temp_part.velocity.x());

        //速度限制
        if (temp_part.velocity.x() > VelMax.x())
            temp_part.velocity.x() = VelMax.x();
        else if(temp_part.velocity.x() < VelMin.x())
            temp_part.velocity.x() = VelMin.x();
        //粒子的位置更新公式
        temp_part.position.x() += line_X * dt;
        temp_part.position.y() += line_Y * dt;
        //边界判断
        if (temp_part.position.x() < PosMin.x()) {
            temp_part.position.x() = PosMin.x();
            temp_part.velocity.y() = pi * (0.5 - (rand() % 1000) / 1000.0);//角度
        }
        else if (temp_part.position.x() > PosMax.x()) {
            temp_part.position.x() = PosMax.x();
            temp_part.velocity.y() = pi * (0.5 + (rand() % 1000) / 1000.0);//角度
        }

        if (temp_part.position.y() < PosMin.y()) {
            temp_part.position.y() = PosMin.y();
            temp_part.velocity.y() = pi * ((rand() % 1000) / 1000.0);//角度
        }
        else if (temp_part.position.y() > PosMax.y()) {
            temp_part.position.y() = PosMax.y();
            temp_part.velocity.y() = pi * (1 + (rand() % 1000) / 1000.0);//角度
        }

        //概率更新
        int t = global_map[int(temp_part.position.y()/resolution)][int(temp_part.position.y()/resolution)]->search_count[uav->id];
        for (int k = 0; k < target_num; k++){
            temp_part.p[k] = (1 - exp(-tao *t )) / (global_map.size()-uav->coverd_area_cnt);
        }
        //适应值更新
        double temp_fitness_1 = 0, temp_fitness_2 = 0;
        for (int k = 0; k < target_num; k++) {
            double Cjk = 0;
            if (dist(target[k]->position, uav->position) < search_R)
                Cjk =  (dist(target[k]->position, temp_part.position));//无人机k与存在目标j的匹配程度
            else
                Cjk = 1;
            double Dik =   dist(uav->position, temp_part.position);//无人机k与粒子i的距离
            double Pij = temp_part.p[k];

            int Tj = uav->Tj[k];//目标k是否未被与之更匹配的无人机跟踪
            //temp_fitness_1 += Tj * (Pij  + Dik +Cjk);
            temp_fitness_1 += Pij *  Tj / (Dik* Cjk);
        }

        if(t>=forget_time)
            temp_fitness_2 = w * 1.0 / (global_map.size() - uav->coverd_area_cnt);
        else
            temp_fitness_2 = (1 - w)*t / (forget_time * (global_map.size() - uav->coverd_area_cnt));

        double temp_fitness = w1 *1000* temp_fitness_1 + w2 *1000* temp_fitness_2;//适应值
        //cout << "temp_fitness = " << temp_fitness <<endl;
        if (temp_fitness >= temp_part.fitness) {//局部最优位置判断
            temp_part.fitness = temp_fitness;
            temp_part.Pbest_position = temp_part.position;
        }

        if (temp_fitness >= uav->Gbest_fitness) {//全局最优位置判断
            //cout << "更新temp_fitness" << temp_fitness << " current particle" << j << endl;
            uav->Gbest_fitness = temp_fitness;
            uav->Gbest_position = temp_part.position;
        }
    }

    if (uav->last_Gbest_position.x() == uav->Gbest_position.x() && uav->last_Gbest_position.y() == uav->Gbest_position.y())
        uav->isConvergenced = true;
    return;
}

bool init() {
    /*for (int i = 0; i < target_num; i++) {//目标初始化
        TARGET temp_target;
        temp_target.position.x = 100+i*i*100;//PosMin.x + (PosMax.x - PosMin.x) * (double)rand() / RAND_MAX;10+i*100
        temp_target.position.y = 50+i*i*300;//PosMin.y + (PosMax.y - PosMin.y) * (double)rand() / RAND_MAX;30+i*200
        temp_target.velocity.x = 0;//线速度
        temp_target.velocity.y = 2 * pi*rand() / RAND_MAX;//角度
        temp_target.accelerate.x = 0;//线加速度
        temp_target.accelerate.y = 0;//角速度
        target.push_back(temp_target);

    }*/

    TARGET* temp_target = new TARGET();
    temp_target->position.x() = 20000;//PosMin.x + (PosMax.x - PosMin.x) * (double)rand() / RAND_MAX;10+i*100
    temp_target->position.y() = 50000;//PosMin.y + (PosMax.y - PosMin.y) * (double)rand() / RAND_MAX;30+i*200
    temp_target->velocity.x() = 150;//线速度
    temp_target->velocity.y() = 2 * pi*rand() / RAND_MAX;//角度
    temp_target->accelerate.x() = 0;//线加速度
    temp_target->accelerate.y() = 0;//角速度
    target.push_back(temp_target);

    temp_target->position.x() = 30000;//PosMin.x + (PosMax.x - PosMin.x) * (double)rand() / RAND_MAX;10+i*100
    temp_target->position.y() = 60000;//PosMin.y + (PosMax.y - PosMin.y) * (double)rand() / RAND_MAX;30+i*200
    temp_target->velocity.x() = 150;//线速度
    temp_target->velocity.y() = 2 * pi*rand() / RAND_MAX;//角度
    temp_target->accelerate.x() = 0;//线加速度
    temp_target->accelerate.y() = 0;//角速度
    target.push_back(temp_target);

    temp_target->position.x() = 50000;//PosMin.x + (PosMax.x - PosMin.x) * (double)rand() / RAND_MAX;10+i*100
    temp_target->position.y() = 70000;//PosMin.y + (PosMax.y - PosMin.y) * (double)rand() / RAND_MAX;30+i*200
    temp_target->velocity.x() = 150;//线速度
    temp_target->velocity.y() = 2 * pi*rand() / RAND_MAX;//角度
    temp_target->accelerate.x() = 0;//线加速度
    temp_target->accelerate.y() = 0;//角速度*/
    target.push_back(temp_target);

    for (int i = 0; i < size_y ; i++) {//初始化地图，i=0代表第一行，j=0代表第一列
        vector<grid*> temp;
        for (int j = 0; j < size_x ; j++) {
            grid* tempGridPtr = new grid();
            for(int k=0;k<UAV_num;k++)
                tempGridPtr->search_count.push_back(forget_time);
            temp.push_back(tempGridPtr);
        }
        global_map.push_back(temp);
    }

    for (int i = 0; i < UAV_num; i++) {//初始化无人机
        UAV* uav_temp = new UAV();
        srand((int(time(NULL)) + i));
        uav_temp->id =1;
        uav_temp->position.x() = 0;//位置
        uav_temp->position.y() = 10000*i;
        uav_temp->velocity.x() = VelMin.x() + (VelMax.x() - VelMin.x()) * (rand() / double(RAND_MAX));//线速度
        uav_temp->velocity.y() = 0;//角度
        uav_temp->accelerate.x() = 0;//线加速度
        uav_temp->accelerate.y() = 0;//角速度
        uav_temp->search_r = search_R;
        uav_temp->particle_r = particle_R;
        updateSubMap(uav_temp);
        for (int j = 0; j < particle_num; j++) {
            particle* temp = new particle();
            uav_temp->swarm.push_back(temp);
        }
        spreadParticles(uav_temp);
        uav.push_back(uav_temp);
    }

    for(auto uav_ptr:uav){
        uav_ptr->traj_Point = uav_ptr->Gbest_position;
        updateUAVStatesInDubinsState(uav_ptr);
        uav_ptr->isConvergenced = false;
        for (int i = 0; i < target_num; i++) {
            Eigen::Vector2d target_pose(-1,-1);
            uav_ptr->target_position.push_back(make_pair(target_pose,forget_time));
            uav_ptr->Tj[i] = 1;
        }
        uav_ptr->track_target_num = -1;
    }
    for (int i = 0; i < UAV_num; i++) {//reset
        union_.parent.push_back(i);
    }
    return 1;
};

void updateTargetStates() {
    //目标状态更新，这里是等速运动
    output_target.open(target_path.c_str(), ios::app | ios::binary);
    for (int target_tag = 0; target_tag < target_num; target_tag++) {
        double line_X = target[target_tag]->velocity.x() * cos(target[target_tag]->velocity.y());//速度在XY的分量
        double line_Y = target[target_tag]->velocity.x() * sin(target[target_tag]->velocity.y());
        target[target_tag]->position.x() += line_X * dt;//位置更新
        target[target_tag]->position.y() += line_Y * dt;

        //边界判断
        if (target[target_tag]->position.x() < PosMin.x()) {
            line_X = -line_X;
            target[target_tag]->position.x() = PosMin.x();
        }
        else if (target[target_tag]->position.x() > PosMax.x()) {
            line_X = -line_X;
            target[target_tag]->position.x() = PosMax.x();
        }

        if (target[target_tag]->position.y() < PosMin.y()) {
            line_Y = -line_Y;
            target[target_tag]->position.y() = PosMin.y();
        }
        else if (target[target_tag]->position.y() > PosMax.y()) {
            line_Y = -line_Y;
            target[target_tag]->position.y() = PosMax.y();
        }
        output_target << target[target_tag]->position.x() << " " << target[target_tag]->position.y() << " ";
    }
    output_target << endl;
    output_target.close();
    return;
};

void updateUAVStates1(){
    output_uav.open(uav_path.c_str(), ios::app | ios::binary);
    output_traj_Point.open(traj_Point_path.c_str(), ios::app | ios::binary);
    for(int i =0;i<UAV_num;i++){

        uav[i]->position.x() = uav[i]->path_.front()[0];
        uav[i]->position.y() = uav[i]->path_.front()[1];
        uav[i]->velocity.x() = uav[i]->path_.front()[2];
        uav[i]->velocity.y() = uav[i]->path_.front()[3];
        uav[i]->path_.pop();
        output_uav << uav[i]->position.x() << " " << uav[i]->position.y() << " " << uav[i]->track_target_num << " ";
        output_traj_Point << uav[i]->traj_Point.x() << " " << uav[i]->traj_Point.y() << " ";
        updateSubMap(uav[i]);
        for(int j=0;j<target_num;j++){//target inf0
            if(dist(target[j]->position,uav[i]->position)<search_R) {
                uav[i]->target_position[j].first = target[j]->position;
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
    for (int i = 0; i < UAV_num-1; i++) {
        for (int j = i+1; j < UAV_num; j++) {
            if(dist(uav[i]->position, uav[j]->position) <= communication_R){//union
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
                }
            }
        }
    }
    return;
}

void updateMission() {
    for (int i = 0; i < UAV_num; i++) {
        for (int j = 0; j < target_num; j++) {//如果无人机i发现的目标j所在的栅格的搜索次数大于1000次
            //cout << "uav_temp.covered_target_id[" << i << "].first = " << uav[i].covered_target_id[i].first << ";second = " << uav[i].covered_target_id[i].second << endl;
            if (uav[i]->target_position[j].second >= forget_time)
                target_state[i][j] = 0;//不让无人机i去追目标j
            else
                target_state[i][j] = 1000.0 / dist(uav[i]->position, target[j]->position);
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
        if (maxTargetState > 0) {//当能够找到有效的任务分配时，执行任务分配
            cout << "uav_cnt = " << uav_cnt << " target_cnt = " << target_cnt << endl;
            if (dist(uav[uav_cnt]->position, target[target_cnt]->position) < 1000) {
                tracked[target_cnt] = true;
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
                uav[i]->traj_Point = target[j]->position;
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
                spreadParticles(uav[i]);
                updateParticleStates(uav[i]);
                uav[i]->traj_Point = uav[i]->Gbest_position;
                updateUAVStatesInDubinsState(uav[i]);
            }
        }
        //更新目标位置
        updateTargetStates();
        //更新无人机状态
        updateUAVStates1();
        //机间通信
        informationShare();
        //任务分配
        updateMission();
        cout << "cunt" << cunt << endl;
        //判断终止
        if (tracked[0] && tracked[1] && tracked[2])
            break;
        cunt++;
    }
    //std::cout << "target_state:" << endl << target_state << endl;
    std::cout << "All target has been tracked!The total number of step is " << cunt << std::endl;
    char aa = getchar();
    return 0;
}