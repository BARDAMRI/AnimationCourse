#pragma once

#include "Scene.h"
#include <utility>

#include <igl/collapse_edge.h>

#include "igl/opengl/glfw/Viewer.h"
#include "igl/aabb.h"
#include <igl/circulation.h>
#include <igl/collapse_edge.h>
#include <igl/edge_flaps.h>
#include <igl/decimate.h>
#include <igl/shortest_edge_and_midpoint.h>
#include <igl/parallel_for.h>
#include <igl/read_triangle_mesh.h>
#include <Eigen/Core>
#include <iostream>
#include <set>
#include <igl/vertex_triangle_adjacency.h>

class BasicScene : public cg3d::Scene
{
public:
    explicit BasicScene(std::string name, cg3d::Display* display) : Scene(std::move(name), display) {};
    void Init(float fov, int width, int height, float near, float far);
    void Update(const cg3d::Program& program, const Eigen::Matrix4f& proj, const Eigen::Matrix4f& view, const Eigen::Matrix4f& model) override;
    void KeyCallback(cg3d::Viewport* _viewport, int x, int y, int key, int scancode, int action, int mods) override;
    void Simp();
    void reset();
    void caculateQMatrix();
    void caculateCostAndPlacment(int edge);

private:
    std::vector <Eigen::Matrix4d> Qmatrix;
    std::shared_ptr<Movable> root;
    std::shared_ptr<cg3d::Model> cyl, sphere1 ,cube, curr_mesh;
    std::string curr_mesh_name;
    Eigen::VectorXi EMAP;
    Eigen::MatrixXi F,IF,E,EF,EI;
    Eigen::VectorXi EQ;
    // If an edge were collapsed, we'd collapse it to these points:
    Eigen::MatrixXd V,IV, C;
    igl::min_heap< std::tuple<double,int,int> > Q;
    typedef std::set<std::pair<double, int> > PriorityQueue;
    std::vector<PriorityQueue::iterator > Qit;
    int num_collapsed;
    PriorityQueue Qt;
    std::vector < Eigen::MatrixXd> prev_V;
    std::vector < Eigen::MatrixXi> prev_F;
    int curr_collaps_index = 0;
    int max = 20;
    std::vector <bool> done_collaps;
    bool collapse_edge();
    void set_mesh( Eigen::MatrixXd &OV,  Eigen::MatrixXi &OF);
    void Simp2();


};
