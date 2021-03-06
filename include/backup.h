//
// Created by bb on 2020/3/9.
//

#ifndef SRC_UTILITY_H
#define SRC_UTILITY_H



#include <vector>
#include <map>
#include <Eigen/Eigen>
#include <queue>
#include <ompl/base/spaces/ReedsSheppStateSpace.h>
#include <ompl/base/spaces/DubinsStateSpace.h>
#include <ompl/base/spaces/SE2StateSpace.h>
#include <ompl/base/State.h>

class point2D {
public:
    point2D() {

    };
    point2D(double a, double b) {
        x = a; y = b;
    };

    ~point2D() {};

    double x;
    double y;
};



class TARGET
{
public:
    TARGET() {};
    ~TARGET() {};

    point2D position;
    point2D velocity;
    point2D accelerate;
};

class particle
{
public:
    particle() {};
    ~particle() {};
    point2D position;
    point2D velocity;
    point2D accelerate;

    double p[3];

    double fitness;
    point2D Pbest_position;

};

class grid
{
public:
    grid() {};
    ~grid() {};
    point2D position;
    int search_cunt[9];
    bool has_target[9][3];
};

class UAV
{
public:
    UAV() {};
    ~UAV() {};

    point2D position;

    point2D velocity;

    point2D accelerate;

    double search_r;

    double particle_r;

    point2D Gbest_position;
    point2D last_Gbest_position;
    point2D traj_Point;
    double Gbest_fitness;

    std::queue<std::vector<double>> path_;//states

    std::vector<particle> swarm;
    std::map<int,int> covered_grid_id;////第一项栅格ID,第二项被跟踪次数
    std::map<int,std::pair<int,int> > covered_target_id;//第一项被跟踪的目标ID，第二项(栅格ID,该栅格已被搜索次数)
    std::map<int, int> Tj;//第一项被跟踪的目标ID，第二项该目标是否被更好的无人机追踪，0,1
    int track_target_num;
    bool isConvergenced;//

};

#endif //SRC_UTILITY_H
