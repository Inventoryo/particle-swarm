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

const point2D PosMin(0, 0);
const point2D PosMax(100000, 100000);

const point2D VelMax(200, pi / 6);
const point2D VelMin(150, -pi / 6);
double AccMax = 0.6;
double AccMin = -0.6;
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

vector<UAV> uav;
vector<TARGET> target;
map<int, grid> global_map;
int size_x = (PosMax.x - PosMin.x) / resolution + 1;
int size_y = (PosMax.y - PosMin.y) / resolution + 1;

vector<vector<double>> target_state(UAV_num, vector<double>(target_num, 0));
bool tracked[3] = { false };



double dist(point2D a, point2D b) {
    return sqrt((a.x - b.x)*(a.x - b.x) + (a.y - b.y)*(a.y - b.y));
}


void updateUAVStatesInDubinsState(UAV& uav){

    double theTi = acos((uav.traj_Point.x - uav.position.x) / dist(uav.traj_Point, uav.position));
    HybridAStar::DubinsPath* path = new HybridAStar::DubinsPath();
    double minR = uav.velocity.x * uav.velocity.x / (g * tanTheta);//最小转弯半径
    double q0[3]={uav.position.x,uav.position.y,uav.velocity.y};//start point state
    double q1[3]={uav.traj_Point.x,uav.traj_Point.y,theTi};//end point state;
    HybridAStar::dubins_init(q0,q1,minR,path);
    double velocity = uav.velocity.x;
    double t=dt*velocity,total_lenth = HybridAStar::dubins_path_length(path);

    while(!uav.path_.empty())
        uav.path_.pop();

    while(t<=total_lenth){
        double temp_point[3];
        HybridAStar::dubins_path_sample(path,t,temp_point);
        vector<double> temp;
        temp.push_back(temp_point[0]);
        temp.push_back(temp_point[1]);
        temp.push_back(velocity);
        temp.push_back(temp_point[2]);
        uav.path_.push(temp);
        if(t>path->param[0] && t<path->param[1]&&velocity+AccMax*dt<VelMax.x){
            t+=velocity*dt+0.5*AccMax*dt*dt;
            velocity=min(VelMax.x,velocity+AccMax*dt);
        }else{
            t+=velocity*dt;
        }
    }

}


/* 根据初始点、速度矢量和圆半径，计算圆心坐标，但对应圆心坐标有2个！！！！！！！！ */
void calCircleCenter_right(const point2D &beginPoint, const double radius, const point2D &velocity, point2D &circleCenter) {
    //double vel_x = beginPoint.x + velocity.x * cos(velocity.y);
    //double vel_y = beginPoint.y + velocity.x * sin(velocity.y);
    //const point2D vel_point(vel_x, vel_y);
    // 设力与速度右手定则90°！！！！设了其中一个！！！！
    double direction_center_sin = -cos(velocity.y);
    double direction_center_cos = sin(velocity.y);
    double center_x = beginPoint.x + direction_center_cos * radius;
    double center_y = beginPoint.y + direction_center_sin * radius;
    circleCenter.x = center_x;
    circleCenter.y = center_y;
}

/* 根据初始点、速度矢量和圆半径，计算圆心坐标，但对应圆心坐标有2个！！！！！！！！ */
void calCircleCenter_left(const point2D& beginPoint, const double radius, const point2D& velocity, point2D& circleCenter) {
    //double vel_x = beginPoint.x + velocity.x * cos(velocity.y);
    //double vel_y = beginPoint.y + velocity.x * sin(velocity.y);
    //const point2D vel_point(vel_x, vel_y);
    // 设力与速度左手定则90°！！！！设了其中一个！！！！
    double direction_center_sin = cos(velocity.y);
    double direction_center_cos = -sin(velocity.y);
    double center_x = beginPoint.x + direction_center_cos * radius;
    double center_y = beginPoint.y + direction_center_sin * radius;
    circleCenter.x = center_x;
    circleCenter.y = center_y;
}

// 根据目标点goalPoint，确定飞行模式：顺时针+直飞=0，顺时针=1，直飞=2，逆时针+直飞=3，逆时针=4
enum flyPattern {
    //	clockwise_straight = 0,
            clockwise = 1,
    straight = 2,
    counterclockwise = 3,
    //	counterclockwise_straight = 4,
};

void chooseFlyPattern(const point2D& beginPoint, const point2D& velocity, const point2D& goalPoint, int& fly_pattern) {
    double thetaTI;//无人机与目标的方位角
    if ((goalPoint.y - beginPoint.y) >= 0)
        thetaTI = acos((goalPoint.x - beginPoint.x) / dist(goalPoint, beginPoint));
    else
        thetaTI = 2 * pi - acos((goalPoint.x - beginPoint.x) / dist(goalPoint, beginPoint));
    double thetaVT = velocity.y;//无人机的速度方向角
    double thetaVTI = thetaTI - thetaVT;
    if (thetaVTI > pi) {
        thetaVTI -= 2 * pi;
    }

    else if (thetaVTI < -pi) {
        thetaVTI = thetaVTI + 2 * pi;
    }
    if (thetaVTI > pi / 180) {
        fly_pattern = counterclockwise;
    }
    else if (thetaVTI < -pi / 180) {
        fly_pattern = clockwise;
    }
    else {
        fly_pattern = straight;
    }
}


// 计算切点坐标
void calCutOffPoint(const point2D& ptCenter, const point2D& ptOutside, const double dbRadious, const point2D& beginPoint, const point2D& velocity, point2D& cutOffPoint) {
    // 计算两个符合条件的切点
    point2D E, F, G, cutOffPoint_1, cutOffPoint_2;
    double r = dbRadious;
    //1. 坐标平移到圆心ptCenter处,求园外点的新坐标E
    E.x = ptOutside.x - ptCenter.x;
    E.y = ptOutside.y - ptCenter.y; //平移变换到E

    //2. 求园与OE的交点坐标F, 相当于E的缩放变换
    double t = r / sqrt(E.x * E.x + E.y * E.y);  //得到缩放比例
    F.x = E.x * t;   F.y = E.y * t;   //缩放变换到F

    //3. 将E旋转变换角度a到切点G，其中cos(a)=r/OF=t, 所以a=arccos(t);
    if (t > 0.9999999) t = 1;
    else if (t < -0.9999999) t = -1;
    double a = acos(t);   //得到旋转角度
    G.x = F.x * cos(a) - F.y * sin(a);
    G.y = F.x * sin(a) + F.y * cos(a);    //旋转变换到G

    //4. 将G平移到原来的坐标下得到新坐标H
    cutOffPoint_1.x = G.x + ptCenter.x;
    cutOffPoint_1.y = G.y + ptCenter.y;             //平移变换到H
    //cout << "  cutoff point_1: (" << cutOffPoint_1.x << ", " << cutOffPoint_1.y << ") " << endl;

    //3. 将E旋转变换角度a到切点G，其中cos(a)=r/OF=t, 所以a=-arccos(t);
    a = -acos(t);   //得到旋转角度
    G.x = F.x * cos(a) - F.y * sin(a);
    G.y = F.x * sin(a) + F.y * cos(a);    //旋转变换到G

    //4. 将G平移到原来的坐标下得到新坐标H
    cutOffPoint_2.x = G.x + ptCenter.x;
    cutOffPoint_2.y = G.y + ptCenter.y;             //平移变换到H
    cout << "  cutoff point_2: (" << cutOffPoint_2.x << ", " << cutOffPoint_2.y << ") " << endl;

    //计算两个切点中哪个切点符合要求
    point2D vel_vector(velocity.x * cos(velocity.y), velocity.x * sin(velocity.y));
    point2D cutoff_vec_1(cutOffPoint_1.x - beginPoint.x, cutOffPoint_1.y - beginPoint.y);
    point2D cutoff_vec_2(cutOffPoint_2.x - beginPoint.x, cutOffPoint_2.y - beginPoint.y);

    double a1 = vel_vector.x * cutoff_vec_1.x + vel_vector.y * cutoff_vec_1.y;
    double b1 = sqrt(vel_vector.x * vel_vector.x + vel_vector.y * vel_vector.y);
    double c1 = sqrt(cutoff_vec_1.x * cutoff_vec_1.x + cutoff_vec_1.y * cutoff_vec_1.y);
    double cos_1 = a1 / b1 / c1;

    double a2 = vel_vector.x * cutoff_vec_2.x + vel_vector.y * cutoff_vec_2.y;
    double b2 = sqrt(vel_vector.x * vel_vector.x + vel_vector.y * vel_vector.y);
    double c2 = sqrt(cutoff_vec_2.x * cutoff_vec_2.x + cutoff_vec_2.y * cutoff_vec_2.y);
    double cos_2 = a2 / b2 / c2;

    cutOffPoint = cos_1 > cos_2 ? cutOffPoint_1 : cutOffPoint_2;

}


void spreadParticles(UAV & uav) {

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

    for (int i = -1; i <= (PosMax.y - PosMin.y) / resolution + 1; i++) {//初始化地图，i=0代表第一行，j=0代表第一列
        for (int j = -1; j <= (PosMax.x - PosMin.x) / resolution + 1; j++) {
            grid temp_grid;
            int id = i * size_x + j;
            temp_grid.position.x = j * resolution + resolution / 2.0;
            temp_grid.position.y = i * resolution + resolution / 2.0;
            global_map[id] = temp_grid;
        }
    }
    for (int i = 0; i < UAV_num; i++) {//初始化无人机
        UAV uav_temp;
        srand((int(time(NULL)) + i));
        uav_temp.position.x = 0;//位置
        uav_temp.position.y = 10000*i;
        uav_temp.velocity.x = VelMin.x + (VelMax.x - VelMin.x) * (rand() / double(RAND_MAX));//线速度
        uav_temp.velocity.y = 0;//角度
        uav_temp.accelerate.x = 0;//线加速度
        uav_temp.accelerate.y = 0;//角速度
        uav_temp.search_r = search_R;
        uav_temp.particle_r = particle_R;
        for (int j = 0; j < particle_num; j++) {//初始化粒子
            particle  temp_particle;
            double length = 2*search_R + VelMax.x * dt * (double)rand() / RAND_MAX;
            double theta = 2 * pi*rand() / RAND_MAX;
            temp_particle.position.x = 2*length * cos(theta) + uav_temp.position.x;
            temp_particle.position.y =2*length * sin(theta) + uav_temp.position.y;
            temp_particle.velocity.x = VelMax.x + (VelMax.x - VelMin.x) * (double)rand() / RAND_MAX;//线速度
            temp_particle.velocity.y = 2 * pi*rand() / RAND_MAX;//角度
            temp_particle.accelerate.x = 0;//线加速度
            temp_particle.accelerate.y = 0;//角速度

            temp_particle.fitness = 0;
            //边界判断
            if (temp_particle.position.x < PosMin.x) {
                temp_particle.position.x = PosMin.x;
                temp_particle.velocity.y = pi * (0.5 - (rand() % 1000) / 1000.0);//角度
            }
            else if (temp_particle.position.x > PosMax.x) {
                temp_particle.position.x = PosMax.x;
                temp_particle.velocity.y = pi * (0.5 + (rand() % 1000) / 1000.0);//角度
            }

            if (temp_particle.position.y < PosMin.y) {
                temp_particle.position.y = PosMin.y;
                temp_particle.velocity.y = pi * ((rand() % 1000) / 1000.0);//角度
            }
            else if (temp_particle.position.y > PosMax.y) {
                temp_particle.position.y = PosMax.y;
                temp_particle.velocity.y = pi * (1 + (rand() % 1000) / 1000.0);//角度
            }
            temp_particle.Pbest_position = temp_particle.position;//局部最优
            uav_temp.swarm.push_back(temp_particle);
        }
        uav_temp.Gbest_position = uav_temp.swarm.front().Pbest_position;//初始化全局最优
        uav_temp.Gbest_fitness = 0;

        uav_temp.traj_Point = uav_temp.Gbest_position;
        updateUAVStatesInDubinsState(uav_temp);
        uav_temp.isConvergenced = false;
        for (int i = 0; i < target_num; i++) {
            uav_temp.covered_target_id[i] = make_pair(0, forget_time);//初始化无人机认为目标所在位置
            cout << "uav_temp.covered_target_id[" << i << "].first = " << uav_temp.covered_target_id[i].first << ";second = " << uav_temp.covered_target_id[i].second << endl;
            uav_temp.Tj[i] = 1;
        }
        uav_temp.track_target_num = -1;
        //更新地图信息
        for (int j = 0; j < (2 * search_R) / resolution; j++) {//j是在无人机搜索半径划分的正方形区域内所在的行数
            for (int k = 0; k < (2 * search_R) / resolution; k++) {//列

                int id = (static_cast<int>(uav_temp.position.y - search_R) / resolution + j) *size_x + static_cast<int>(uav_temp.position.x - search_R) / resolution + k;//栅格在全局地图中的ID
                if (global_map.count(id)) {//如果此栅格在全局地图里
                    grid& temp_grid = global_map[id];//通过ID在全局地图中找出该栅格，并拿出
                    if (dist(temp_grid.position, uav_temp.position) < search_R) //如果栅格在无人机探测半径内
                    {
                        for (int l = 0; l < target_num; l++) {//栅格内目标状态更新
                            if (dist(temp_grid.position, target[l].position) < resolution*sqrt(0.5)) {
                                uav_temp.covered_target_id[l] = make_pair(id, 0);
                                cout << "!" << endl;
                            } //如果当前栅格有目标

                        }
                        uav_temp.covered_grid_id[id] = 0;
                    }
                }
            }
        }
        //spreadParticles(uav_temp);
        uav.push_back(uav_temp);
    }
    return 0;
};


void updateParticleStates(){//更新粒子

    for (int i = 0; i < UAV_num;i++) {
        //uav[i].Gbest_fitness = 0;
        if (uav[i].track_target_num != -1 && tracked[uav[i].track_target_num])//如果无人机已经追踪到目标，则粒子群更新停止
            continue;
        //if (uav[i].isConvergenced)//如果粒子群已收敛，则不更新
        //continue;
        uav[i].last_Gbest_position = uav[i].Gbest_position;
        for (int j = 0; j < particle_num; j++) {
            particle & temp_part = uav[i].swarm[j];
            //速度更新
            double r1 = (double)rand() / RAND_MAX;
            double r2 = (double)rand() / RAND_MAX;
            double line_X = temp_part.velocity.x * cos(temp_part.velocity.y);//速度在XY向的分量
            double line_Y = temp_part.velocity.x * sin(temp_part.velocity.y);
            line_X = weight * line_X + c1 * r1 *(temp_part.Pbest_position.x - temp_part.position.x) + c2 * r2 *(uav[i].Gbest_position.x - temp_part.position.x);
            line_Y = weight * line_Y + c1 * r1 *(temp_part.Pbest_position.y - temp_part.position.y) + c2 * r2 *(uav[i].Gbest_position.y - temp_part.position.y);
            temp_part.velocity.x = sqrt(line_X * line_X + line_Y * line_Y);//粒子的速度更新公式
            if (line_Y >= 0)
                temp_part.velocity.y = acos(line_X / temp_part.velocity.x);
            else
                temp_part.velocity.y = 2*pi - acos(line_X / temp_part.velocity.x);

            //速度限制
            if (temp_part.velocity.x > VelMax.x)
                temp_part.velocity.x = VelMax.x;
            else if(temp_part.velocity.x < VelMin.x)
                temp_part.velocity.x = VelMin.x;
            //粒子的位置更新公式
            temp_part.position.x += line_X * dt;
            temp_part.position.y += line_Y * dt;
            //边界判断
            if (temp_part.position.x < PosMin.x) {
                temp_part.position.x = PosMin.x;
                temp_part.velocity.y = pi * (0.5 - (rand() % 1000) / 1000.0);//角度
            }
            else if (temp_part.position.x > PosMax.x) {
                temp_part.position.x = PosMax.x;
                temp_part.velocity.y = pi * (0.5 + (rand() % 1000) / 1000.0);//角度
            }

            if (temp_part.position.y < PosMin.y) {
                temp_part.position.y = PosMin.y;
                temp_part.velocity.y = pi * ((rand() % 1000) / 1000.0);//角度
            }
            else if (temp_part.position.y > PosMax.y) {
                temp_part.position.y = PosMax.y;
                temp_part.velocity.y = pi * (1 + (rand() % 1000) / 1000.0);//角度
            }


            //概率更新
            int id = static_cast<int>(temp_part.position.y / resolution) * size_x + static_cast<int>(temp_part.position.x / resolution);
            if (global_map.count(id)) {//如果全局地图中有此id
                for (int k = 0; k < target_num; k++){
                    if (uav[i].covered_grid_id.find(id) != uav[i].covered_grid_id.end())//如果该栅格被搜索过
                        temp_part.p[k] = (1 - exp(-tao * uav[i].covered_grid_id[id]))*(1.0 / (global_map.size() - uav[i].covered_grid_id.size()));
                    else//如果该栅格没有被搜索过
                        temp_part.p[k] =  (1 - exp(-tao*forget_time )) / (global_map.size() - uav[i].covered_grid_id.size());
                }
            }
            else
                for (int k = 0; k < target_num; k++) {
                    temp_part.p[k] = 0;
                }
            //适应值更新
            double temp_fitness_1 = 0, temp_fitness_2 = 0;
            for (int k = 0; k < target_num; k++) {
                double Cjk = 0;
                if (dist(target[k].position, uav[i].position) < search_R)
                    Cjk =  (dist(target[k].position, temp_part.position));//无人机k与存在目标j的匹配程度
                else
                    Cjk = 1;
                double Dik =   dist(uav[i].position, temp_part.position);//无人机k与粒子i的距离
                double Pij = temp_part.p[k];

                int Tj = uav[i].Tj[k];//目标k是否未被与之更匹配的无人机跟踪
                //temp_fitness_1 += Tj * (Pij  + Dik +Cjk);
                temp_fitness_1 += Pij *  Tj / (Dik* Cjk);
            }
            if (uav[i].covered_grid_id.find(id) == uav[i].covered_grid_id.end()) {//如果该栅格被搜索过
                temp_fitness_2 = w * 1.0 / (global_map.size() - uav[i].covered_grid_id.size());
            }
            else {
                temp_fitness_2 = (1 - w)*uav[i].covered_grid_id[id] / (forget_time * (global_map.size() - uav[i].covered_grid_id.size()));
            }
            double temp_fitness = w1 *1000* temp_fitness_1 + w2 *1000* temp_fitness_2;//适应值
            //cout << "temp_fitness = " << temp_fitness <<endl;
            if (temp_fitness >= temp_part.fitness) {//局部最优位置判断
                temp_part.fitness = temp_fitness;
                temp_part.Pbest_position = temp_part.position;
            }

            if (temp_fitness >= uav[i].Gbest_fitness) {//全局最优位置判断
                //cout << "更新temp_fitness" << temp_fitness << " current particle" << j << endl;
                uav[i].Gbest_fitness = temp_fitness;
                uav[i].Gbest_position = temp_part.position;
            }
        }

        if (uav[i].last_Gbest_position.x == uav[i].Gbest_position.x&&uav[i].last_Gbest_position.y == uav[i].Gbest_position.y)
            uav[i].isConvergenced = true;
    }
    return;
}


void updateTargetStates() {
    //目标状态更新，这里是等速运动
    output_target.open(target_path.c_str(), ios::app | ios::binary);
    for (int target_tag = 0; target_tag < target_num; target_tag++) {
        double line_X = target[target_tag].velocity.x * cos(target[target_tag].velocity.y);//速度在XY的分量
        double line_Y = target[target_tag].velocity.x * sin(target[target_tag].velocity.y);
        target[target_tag].position.x += line_X * dt;//位置更新
        target[target_tag].position.y += line_Y * dt;

        //边界判断
        if (target[target_tag].position.x < PosMin.x) {
            line_X = -line_X;
            target[target_tag].position.x = PosMin.x;
        }
        else if (target[target_tag].position.x > PosMax.x) {
            line_X = -line_X;
            target[target_tag].position.x = PosMax.x;
        }

        if (target[target_tag].position.y < PosMin.y) {
            line_Y = -line_Y;
            target[target_tag].position.y = PosMin.y;
        }
        else if (target[target_tag].position.y > PosMax.y) {
            line_Y = -line_Y;
            target[target_tag].position.y = PosMax.y;
        }
        output_target << target[target_tag].position.x << " " << target[target_tag].position.y << " ";
    }
    output_target << endl;
    output_target.close();
    return;
};

void updateUAVStates1(){

    output_uav.open(uav_path.c_str(), ios::app | ios::binary);
    output_traj_Point.open(traj_Point_path.c_str(), ios::app | ios::binary);
    for(int i =0;i<UAV_num;i++){

        uav[i].position.x = uav[i].path_.front()[0];
        uav[i].position.y = uav[i].path_.front()[1];
        uav[i].velocity.x = uav[i].path_.front()[2];
        uav[i].velocity.y = uav[i].path_.front()[3];
        uav[i].path_.pop();
        cout<<"updating UAV "<<i<<endl;
        output_uav << uav[i].position.x << " " << uav[i].position.y << " " << uav[i].track_target_num << " ";
        output_traj_Point << uav[i].traj_Point.x << " " << uav[i].traj_Point.y << " ";
    }
    output_uav << endl;
    output_uav.close();
    output_traj_Point << endl;
    output_traj_Point.close();

    //更新地图状态
    for (int i = 0; i < UAV_num; i++) {//无人机
        //旧的栅格点++

        for (auto iter = uav[i].covered_grid_id.begin(); iter != uav[i].covered_grid_id.end(); ++iter) {
            iter->second++;
            if (iter->second >= forget_time) //如果跟踪次数超过阈值，则视为没有被跟踪过
                iter->second = forget_time;
        }

        for (int l = 0; l < target_num; l++) {
            uav[i].covered_target_id[l].second++;
            if (uav[i].covered_target_id[l].second >= forget_time) //如果跟踪次数超过阈值，则视为没有被跟踪过
                uav[i].covered_target_id[l].second = forget_time;
            //cout << "uav_temp.covered_target_id[" << i << "].first = " << uav[i].covered_target_id[i].first << ";second = " << uav[i].covered_target_id[i].second << endl;
        }//目标跟踪信息更新

        //当前视角下的栅格点处理
        for (int j = 0; j < (2 * search_R) / resolution; j++) {//Y向分割
            for (int k = 0; k < (2 * search_R) / resolution; k++) {//X向分割
                if ((uav[i].position.y - search_R + j * resolution) <= PosMax.y && (uav[i].position.y - search_R + j * resolution) >= PosMin.y &&//地图内点
                    (uav[i].position.x - search_R + k * resolution) <= PosMax.x && (uav[i].position.x - search_R + k * resolution) >= PosMin.x) {
                    int id = static_cast<int>((uav[i].position.y - search_R) / resolution + j) *size_x + static_cast<int>((uav[i].position.x - search_R) / resolution) + k;
                    grid& temp_grid = global_map[id];
                    if (dist(temp_grid.position, uav[i].position) < search_R)//当前栅格在无人机视野内
                    {
                        for (int l = 0; l < target_num; l++)
                            if (dist(temp_grid.position, target[l].position) <= resolution * sqrt(0.5)) {
                                uav[i].covered_target_id[l] = make_pair(id, 0);
                                cout << " target " << l << " found in grid " << id << endl;
                            } //当前栅格含目标


                        uav[i].covered_grid_id[id] = 0;//保存该栅格的ID,跟踪次数
                    }
                }
            }
        }
    }
    return;
}

void updateUAVStates() {//更新UAV信息

    //按dubins更新无人机运动状态
    output_uav.open(uav_path.c_str(), ios::app | ios::binary);
    output_traj_Point.open(traj_Point_path.c_str(), ios::app | ios::binary);
    for (int i = 0; i < UAV_num; i++) {
        cout << "uav" << i << ".traj_Point" << uav[i].traj_Point.x << " " << uav[i].traj_Point.y << endl;

        double thetaTI;//无人机与目标的方位角
        if ((uav[i].traj_Point.y - uav[i].position.y) >= 0)
            thetaTI = acos((uav[i].traj_Point.x - uav[i].position.x) / dist(uav[i].traj_Point, uav[i].position));
        else
            thetaTI = 2 * pi - acos((uav[i].traj_Point.x - uav[i].position.x) / dist(uav[i].traj_Point, uav[i].position));
        double minR = uav[i].velocity.x * uav[i].velocity.x / (g * tanTheta);//最小转弯半径
        double omega = g * tanTheta / uav[i].velocity.x;//转弯角速率
        double thetaVT = uav[i].velocity.y;//无人机的速度方向角
        double thetaAI;
        double thetaIncre;
        int fly_pattern = -1;
        double R = sqrt(dist(uav[i].position, uav[i].traj_Point));
        double lanmude = atan((uav[i].traj_Point.y - uav[i].position.y) / (uav[i].traj_Point.x - uav[i].position.x));
        double xPoint = uav[i].velocity.x *cos(uav[i].velocity.y);
        double yPoint = uav[i].velocity.x*sin(uav[i].velocity.y);
        double xTPoint = 0;
        double yTPoint = 0;
        double rPoint = ((uav[i].traj_Point.x - uav[i].position.x)*(xTPoint - xPoint) + (uav[i].traj_Point.y - uav[i].position.y)*(yTPoint - yPoint)) / R;
        double lanPoint = ((uav[i].traj_Point.x - uav[i].position.x)*(xTPoint - xPoint) + (uav[i].traj_Point.y - uav[i].position.y)*(yTPoint - yPoint)) / R;
        point2D circle_center(0, 0);
        point2D cutoff_point;
        point2D T(uav[i].position.x, uav[i].position.y);
        point2D I(uav[i].traj_Point.x, uav[i].traj_Point.y);
        point2D poleT;//极坐标系下T点的坐标
        point2D poleA;//极坐标系下A点的坐标
        chooseFlyPattern(T, uav[i].velocity, I, fly_pattern);
        if (fly_pattern == clockwise) calCircleCenter_right(T, minR, uav[i].velocity, circle_center);
        else if (fly_pattern == counterclockwise) calCircleCenter_left(T, minR, uav[i].velocity, circle_center);
        else  cout << "fly pattern: " << fly_pattern << "  straight!!" << endl;
        switch (fly_pattern) {
            /*case 0:
            {
                calCutOffPoint(circle_center, I, minR, T, uav.velocity, cutoff_point);
                uav.position.x = cutoff_point.x;
                uav.position.y = cutoff_point.y;
                if ((I.y - cutoff_point.y) >= 0)
                    thetaAI = acos((I.x - cutoff_point.x) / dist(I, cutoff_point));

                else
                    thetaAI = 2 * pi - acos((I.x - cutoff_point.x) / dist(I, cutoff_point));

                thetaIncre =  thetaAI- thetaVT;
                if (thetaIncre > 0)
                    thetaIncre = thetaIncre - 2 * pi;

                uav.velocity.y += thetaIncre;//角度改变
                if (uav.velocity.y > 2 * pi) {
                    uav.velocity.y -= 2 * pi;
                }
                else if (uav.velocity.y < 0)
                    uav.velocity.y += 2 * pi;

                double t_remain = dt - thetaIncre*dt / omega;//直飞时间
                //加速直线
                if (uav.velocity.x >= VelMax.x) {//达到最大速度
                    uav.velocity.x = VelMax.x;
                    uav.position.x = uav.position.x + uav.velocity.x * t_remain *cos(uav.velocity.y);
                    uav.position.y = uav.position.y + uav.velocity.x* t_remain *sin(uav.velocity.y);

                }
                else {//尚未达到最大速度，则加速
                    uav.position.x = uav.position.x + (uav.velocity.x *t_remain + 0.5 * AccMax *t_remain *t_remain)*cos(uav.velocity.y);
                    uav.position.y= uav.position.y + (uav.velocity.x *t_remain + 0.5 * AccMax *t_remain *t_remain)*sin(uav.velocity.y);
                    uav.velocity.x = uav.velocity.x + AccMax * t_remain;
                }
                break;
            }*/
            case 1://顺时针
            {
                double direction_center_sin = -cos(uav[i].velocity.y);
                double direction_center_cos = sin(uav[i].velocity.y);
                double center_x = uav[i].position.x + direction_center_cos * minR;
                double center_y = uav[i].position.y + direction_center_sin * minR;
                double circleCenter_x = center_x;
                double circleCenter_y = center_y;
                thetaIncre = omega * dt;
                double polethetaTO = 0.0;
                double TO = uav[i].position.y - circleCenter_y;


                if (TO > 0) {
                    polethetaTO = acos((uav[i].position.x - circleCenter_x) / minR);
                }
                else if (TO < 0) {
                    polethetaTO = 2 * pi - acos((uav[i].position.x - circleCenter_x) / minR);
                }
                point2D poleT(minR*cos(polethetaTO), minR*sin(polethetaTO));
                point2D poleA(minR*cos(polethetaTO - thetaIncre), minR*sin(polethetaTO - thetaIncre));
                uav[i].position.x = poleA.x + circleCenter_x;
                uav[i].position.y = poleA.y + circleCenter_y;
                uav[i].velocity.y -= thetaIncre;//角度改变

                if (uav[i].velocity.y > 2 * pi) {
                    uav[i].velocity.y -= 2 * pi;
                }
                else if (uav[i].velocity.y < 0)
                    uav[i].velocity.y += 2 * pi;
                for (int j = 0; j < target_num; j++) {
                    if (dist(uav[i].position, target[j].position) < uav[i].velocity.x * uav[i].velocity.x / (g * tanTheta))
                        uav[i].velocity.x = VelMin.x;
                }
                break;
            }
            case 2://直飞
            {
                if (uav[i].velocity.x >= VelMax.x) {//达到最大速度
                    uav[i].velocity.x = VelMax.x;
                    uav[i].position.x = uav[i].position.x + uav[i].velocity.x * dt *cos(uav[i].velocity.y);
                    uav[i].position.y = uav[i].position.y + uav[i].velocity.x * dt *sin(uav[i].velocity.y);
                }
                else {//尚未达到最大速度，则加速
                    uav[i].position.x = uav[i].position.x + (uav[i].velocity.x *dt + 0.5 * AccMax *dt *dt)*cos(uav[i].velocity.y);
                    uav[i].position.y = uav[i].position.y + (uav[i].velocity.x *dt + 0.5 * AccMax *dt *dt)*sin(uav[i].velocity.y);
                    uav[i].velocity.x = uav[i].velocity.x + AccMax * dt;
                }
                for (int j = 0; j < target_num; j++) {
                    if (dist(uav[i].position, target[j].position) < uav[i].velocity.x * uav[i].velocity.x / (g * tanTheta))
                        uav[i].velocity.x = VelMin.x;
                }
                break;
            }
                /*case 3:
                {
                    calCutOffPoint(circle_center, I, minR, T, uav.velocity, cutoff_point);
                    uav.position.x = cutoff_point.x;
                    uav.position.y = cutoff_point.y;
                    if ((I.y - cutoff_point.y) >= 0)
                        thetaAI = acos((I.x - cutoff_point.x) / dist(I, cutoff_point));
                    else
                        thetaAI = 2 * pi - acos((I.x - cutoff_point.x) / dist(I, cutoff_point));
                    thetaIncre = thetaAI - thetaVT;
                    if (thetaIncre < 0)
                        thetaIncre = thetaIncre + 2 * pi;
                    uav.position.x = cutoff_point.x;
                    uav.position.y = cutoff_point.y;

                    uav.velocity.y += thetaIncre;//角度改变
                    if (uav.velocity.y > 2 * pi) {
                        uav.velocity.y -= 2 * pi;
                    }
                    else if (uav.velocity.y < 0)
                        uav.velocity.y += 2 * pi;

                    double t_remain = dt - thetaIncre*dt / omega;
                    //加速直线
                    if (uav.velocity.x >= VelMax.x) {//达到最大速度
                        uav.velocity.x = VelMax.x;
                        uav.position.x = uav.position.x + uav.velocity.x * t_remain *cos(uav.velocity.y);
                        uav.position.y = uav.position.y + uav.velocity.x * t_remain *sin(uav.velocity.y);

                    }
                    else {//尚未达到最大速度，则加速
                        uav.position.x = uav.position.x + (uav.velocity.x *t_remain + 0.5 * AccMax *t_remain *t_remain)*cos(uav.velocity.y);
                        uav.position.y = uav.position.y + (uav.velocity.x *t_remain + 0.5 * AccMax *t_remain *t_remain)*sin(uav.velocity.y);
                        uav.velocity.x = uav.velocity.x + AccMax * t_remain;
                    }
                    break;
                }*/
            case 3://逆时针
            {
                double direction_center_sin = cos(uav[i].velocity.y);
                double direction_center_cos = -sin(uav[i].velocity.y);
                double center_x = uav[i].position.x + direction_center_cos * minR;
                double center_y = uav[i].position.y + direction_center_sin * minR;
                double circleCenter_x = center_x;
                double circleCenter_y = center_y;
                thetaIncre = omega * dt;
                double polethetaTO = 0;
                double TO = uav[i].position.y - circleCenter_y;
                if (TO > 0) {
                    polethetaTO = acos((uav[i].position.x - circleCenter_x) / minR);
                }
                else if (TO < 0) {
                    polethetaTO = 2 * pi - acos((uav[i].position.x - circleCenter_x) / minR);
                }

                point2D poleT(minR*cos(polethetaTO), minR*sin(polethetaTO));
                point2D poleA(minR*cos(polethetaTO + thetaIncre), minR*sin(polethetaTO + thetaIncre));
                uav[i].position.x = poleA.x + circleCenter_x;
                uav[i].position.y = poleA.y + circleCenter_y;
                uav[i].velocity.y += thetaIncre;//角度改变
                if (uav[i].velocity.y > 2 * pi) {
                    uav[i].velocity.y -= 2 * pi;
                }
                else if (uav[i].velocity.y < 0)
                    uav[i].velocity.y += 2 * pi;

                for (int j = 0; j < target_num; j++) {
                    if (dist(uav[i].position, target[j].position) < uav[i].velocity.x * uav[i].velocity.x / (g * tanTheta))
                        uav[i].velocity.x = VelMin.x;
                }
                break;
            }
            default:
                break;
        }
        cout << uav[i].position.x << " " << uav[i].position.y << " " << "fly_pattern" << fly_pattern << endl;
        output_uav << uav[i].position.x << " " << uav[i].position.y << " " << uav[i].track_target_num << " ";
        output_traj_Point << uav[i].traj_Point.x << " " << uav[i].traj_Point.y << " ";
    }
    output_uav << endl;
    output_uav.close();
    output_traj_Point << endl;
    output_traj_Point.close();



    //更新地图状态
    for (int i = 0; i < UAV_num; i++) {//无人机
        //旧的栅格点++

        for (auto iter = uav[i].covered_grid_id.begin(); iter != uav[i].covered_grid_id.end(); ++iter) {
            iter->second++;
            if (iter->second >= forget_time) //如果跟踪次数超过阈值，则视为没有被跟踪过
                iter->second = forget_time;
        }

        for (int l = 0; l < target_num; l++) {
            uav[i].covered_target_id[l].second++;
            if (uav[i].covered_target_id[l].second >= forget_time) //如果跟踪次数超过阈值，则视为没有被跟踪过
                uav[i].covered_target_id[l].second = forget_time;
            //cout << "uav_temp.covered_target_id[" << i << "].first = " << uav[i].covered_target_id[i].first << ";second = " << uav[i].covered_target_id[i].second << endl;
        }//目标跟踪信息更新


        //当前视角下的栅格点处理
        for (int j = 0; j < (2 * search_R) / resolution; j++) {//Y向分割
            for (int k = 0; k < (2 * search_R) / resolution; k++) {//X向分割
                if ((uav[i].position.y - search_R + j * resolution) <= PosMax.y && (uav[i].position.y - search_R + j * resolution) >= PosMin.y &&//地图内点
                    (uav[i].position.x - search_R + k * resolution) <= PosMax.x && (uav[i].position.x - search_R + k * resolution) >= PosMin.x) {
                    int id = static_cast<int>((uav[i].position.y - search_R) / resolution + j) *size_x + static_cast<int>((uav[i].position.x - search_R) / resolution) + k;
                    grid& temp_grid = global_map[id];
                    if (dist(temp_grid.position, uav[i].position) < search_R)//当前栅格在无人机视野内
                    {
                        for (int l = 0; l < target_num; l++)
                            if (dist(temp_grid.position, target[l].position) <= resolution * sqrt(0.5)) {
                                uav[i].covered_target_id[l] = make_pair(id, 0);
                                cout << " target " << l << " found in grid " << id << endl;
                            } //当前栅格含目标


                        uav[i].covered_grid_id[id] = 0;//保存该栅格的ID,跟踪次数
                    }
                }
            }
        }
    }
    return;
}


void informationShare() {
    for (int i = 0; i < UAV_num; i++) {
        for (int j = 0; j < UAV_num; j++) {
            if (dist(uav[i].position, uav[j].position) <= communication_R) {//两架无人机在通信半径内
                //交换栅格信息
                for (auto c : uav[j].covered_grid_id)
                    if ((uav[i].covered_grid_id.find(c.first) != uav[i].covered_grid_id.end() && uav[i].covered_grid_id[c.first] > c.second)
                        || (uav[i].covered_grid_id.find(c.first) == uav[i].covered_grid_id.end()))
                        uav[i].covered_grid_id[c.first] = c.second;

                //交换目标信息
                for (int l = 0; l < target_num; l++) //交换目标信息，如果当前i搜索目标的次数大于j，则i的信息过时了，把J的信息给i；
                    if (uav[i].covered_target_id[l].second > uav[j].covered_target_id[l].second)
                        uav[i].covered_target_id[l] = uav[j].covered_target_id[l];
            }
        }
    }
    return;
}

void updateMission() {
    for (int i = 0; i < UAV_num; i++) {
        for (int j = 0; j < target_num; j++) {//如果无人机i发现的目标j所在的栅格的搜索次数大于1000次
            //cout << "uav_temp.covered_target_id[" << i << "].first = " << uav[i].covered_target_id[i].first << ";second = " << uav[i].covered_target_id[i].second << endl;
            if (uav[i].covered_target_id[j].second >= forget_time)
                target_state[i][j] = 0;//不让无人机i去追目标j
            else
                target_state[i][j] = 1000.0 / dist(uav[i].position, target[j].position);
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
            if (dist(uav[uav_cnt].position, target[target_cnt].position) < 1000) {
                tracked[target_cnt] = true;
                cout << "target " << target_cnt << " is successfully tracked by uav " << uav_cnt << endl;
            }

            for (int i = 0; i < UAV_num; i++) {
                if (i == uav_cnt) {
                    uav[i].Tj[target_cnt] = 1;//分配该无人机追踪该目标
                }
                else
                    uav[i].Tj[target_cnt] = 0;//否则分配该无人机不追踪该目标
            }
            row.insert(uav_cnt); col.insert(target_cnt);//将行与列放入set中，下次循环跳过该行与列
        }
    }

    for (int i = 0; i < UAV_num; i++) {
        int j = 0;
        for (; j < target_num; j++) {
            //cout << "uav[" << i << "].TJ[" << j << "]= " << uav[i].Tj[j] << " "<< uav[i].covered_target_id[j].second<< endl;
            if (uav[i].Tj[j] == 1 && uav[i].covered_target_id[j].second < forget_time) {//分配该无人机追踪该目标
                uav[i].track_target_num = j;

                uav[i].traj_Point = target[j].position;
                updateUAVStatesInDubinsState(uav[i]);
                cout << "uav " << i << " is tracking target " << j << endl;
                break;
            }
        }
        if (j >= target_num) {
            uav[i].track_target_num = -1;
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
    while (cunt < 500) {
        for (double i = 0; i < 20; i++) {
            //更新粒子状态
            updateParticleStates();
        }
        //更新粒子状态
        //updateParticleStates();

        //更新目标位置
        updateTargetStates();

        //更新无人机状态
        updateUAVStates1();

        //机间通信
        informationShare();

        //如果无人机到达全局最优，则重新洒粒子
        for (int i = 0; i < UAV_num; i++) {
            if (uav[i].path_.empty()) {
                cout << "spreadParticles" << i << endl;
                uav[i].isConvergenced = false;
                uav[i].traj_Point = uav[i].Gbest_position;
                updateUAVStatesInDubinsState(uav[i]);
                spreadParticles((uav[i]));
            }
        }
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