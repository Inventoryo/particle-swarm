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

    Eigen::Vector2d position;
    Eigen::Vector2d velocity;
    Eigen::Vector2d accelerate;
};

class particle
{
public:
    particle() {};
    ~particle() {};
    Eigen::Vector2d position;
    Eigen::Vector2d velocity;
    Eigen::Vector2d accelerate;

    double p[3];

    double fitness;
    Eigen::Vector2d Pbest_position;

};

class grid
{
public:
    grid() {};
    ~grid() {};
    std::vector<int> search_count;
};

class UAV
{
public:
    UAV() {};
    ~UAV() {};

    int id ;
    Eigen::Vector2d position;

    Eigen::Vector2d velocity;

    Eigen::Vector2d accelerate;

    double search_r;

    double particle_r;
    int coverd_area_cnt;

    Eigen::Vector2d Gbest_position;
    Eigen::Vector2d last_Gbest_position;
    Eigen::Vector2d traj_Point;
    double Gbest_fitness;

    std::queue<std::vector<double>> path_;//states
    std::vector<particle*> swarm;
    std::vector<std::pair<Eigen::Vector2d,int>> target_position;//
    std::map<int, int> Tj;//第一项被跟踪的目标ID，第二项该目标是否被更好的无人机追踪，0,1
    int track_target_num;
    bool isConvergenced;//

};

class nion{
public:
    nion() {};
    ~nion() {};
    std::vector<int> parent;
    std::vector<std::vector<int>> barrel;
    int unionsearch(int root) //查找根结点
    {
        int son, tmp;
        son = root;
        while(root != parent[root]) //寻找根结点
            root = parent[root];
        while(son != root) //路径压缩
        {
            tmp = parent[son];
            parent[son] = root;
            son = tmp;
        }
        return root;
    }

    void join(int root1, int root2) //判断是否连通，不连通就合并
    {
        int x, y;
        x = unionsearch(root1);
        y = unionsearch(root2);
        if(x != y) //如果不连通，就把它们所在的连通分支合并
            parent[y] = x;
    }
    void setup_barrel(){
        barrel.clear();
        for(int i =0;i<parent.size();i++){
            if(i==parent[i]){
                std::vector<int> temp;
                temp.push_back(i);
                barrel.push_back(temp);
            }else{
                for(int j =0;j<barrel.size();j++){
                    if(barrel[j].front()==parent[i]){
                        barrel[j].push_back(i);
                        break;
                    }
                }
            }
        }
    };
};

#endif //SRC_UTILITY_H
