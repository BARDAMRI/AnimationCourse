#include "BasicScene.h"
#include <read_triangle_mesh.h>
#include <utility>
#include "ObjLoader.h"
#include "IglMeshLoader.h"
#include "igl/read_triangle_mesh.cpp"
#include "igl/edge_flaps.h"
#include <igl/vertex_triangle_adjacency.h>
#include "AutoMorphingModel.h"

#include <igl/circulation.h>
#include <igl/collapse_edge.h>
#include <igl/edge_flaps.h>
#include <igl/decimate.h>
#include <igl/shortest_edge_and_midpoint.h>
#include <igl/parallel_for.h>
#include <igl/read_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include <Eigen/Core>
#include <iostream>
#include <set>

#include "igl/edge_flaps.h"
#include "igl/collapse_edge.h"
#include <igl/vertex_triangle_adjacency.h>
#include "Eigen/dense"
#include <igl/circulation.h>
#include <functional>


// #include "AutoMorphingModel.h"

using namespace cg3d;
using namespace std;
using namespace Eigen;
using namespace igl;

void BasicScene::Init(float fov, int width, int height, float near, float far) {
    camera = Camera::Create("camera", fov, float(width) / height, near, far);

    AddChild(root = Movable::Create("root")); // a common (invisible) parent object for all the shapes
    auto daylight{std::make_shared<Material>("daylight", "shaders/cubemapShader")};
    daylight->AddTexture(0, "textures/cubemaps/Daylight Box_", 3);
    auto background{Model::Create("background", Mesh::Cube(), daylight)};
    AddChild(background);
    background->Scale(120, Axis::XYZ);
    background->SetPickable(false);
    background->SetStatic();


    auto program = std::make_shared<Program>("shaders/basicShader");
    auto material{std::make_shared<Material>("material", program)}; // empty material
//    SetNamedObject(cube, Model::Create, Mesh::Cube(), material, shared_from_this());

    material->AddTexture(0, "textures/box0.bmp", 2);


    /*
    auto sphereMesh{IglLoader::MeshFromFiles("sphere_igl", "data/sphere.obj")};
    auto sphere{Model::Create("sphere", sphereMesh, material)};
    auto morphFunc = [](Model* model, cg3d::Visitor* visitor) {
        return model->meshIndex;
    };

    sphere1 = AutoMorphingModel::Create(*sphere, morphFunc);
    sphere1->Scale(7);
    sphere1->showWireframe = true;
    sphere1->Translate({-3, 0, 0});
    root->AddChild(sphere1);
     */

    auto cylMesh{IglLoader::MeshFromFiles("cyl_igl","data/camel_b.obj")};
    cyl = Model::Create( "cyl", cylMesh, material);
    cyl->Translate({3,0,0});
    cyl->Scale(0.12f);
    cyl->showWireframe = true;
    root->AddChild(cyl);

    /*
    auto cubeMesh{IglLoader::MeshFromFiles("cube_igl","data/cube.off")};
    cube = Model::Create( "cube", cubeMesh, material);
    cube->showWireframe = true;
    root->AddChild(cube);
     */

    curr_mesh = cyl;
    curr_mesh_name = "cyl";


    auto mesh = curr_mesh->GetMeshList();


    // Function to reset original mesh and data structures
    IV = mesh[0]->data[0].vertices;
    IF = mesh[0]->data[0].faces;


    camera->Translate(20, Axis::Z);
    core().is_animating = false;

    reset();


    /*
    std::cout<< "faces: \n" << F <<std::endl;
    
    std::cout<< "edges: \n" << E.transpose() <<std::endl;
    std::cout<< "edges to faces: \n" << EF.transpose() <<std::endl;
    std::cout<< "faces to edges: \n "<< EMAP.transpose()<<std::endl;
    std::cout<< "edges indices: \n" << EI.transpose() <<std::endl;
     */

}

void BasicScene::Update(const Program& program, const Eigen::Matrix4f& proj, const Eigen::Matrix4f& view, const Eigen::Matrix4f& model)
{
    Scene::Update(program, proj, view, model);
    program.SetUniform4f("lightColor", 1.0f, 1.0f, 1.0f, 0.5f);
    program.SetUniform4f("Kai", 1.0f, 1.0f, 1.0f, 1.0f);
    if(core().is_animating && !Qt.empty()){
        Simp();
        core().is_animating = false;
    }
    //sphere1->Rotate(0.01f, Axis::XYZ);
}
void BasicScene::KeyCallback(Viewport* _viewport, int x, int y, int key, int scancode, int action, int mods)
{
    if (action == GLFW_PRESS || action == GLFW_REPEAT)
    {
        if (key == GLFW_KEY_LEFT){
            if(curr_collaps_index < max-1){
                curr_mesh->meshIndex++;
                curr_collaps_index++;
            }
            core().is_animating =true ;
        }
        if (key == GLFW_KEY_RIGHT){
            if(curr_collaps_index > 0){
                curr_mesh->meshIndex--;
                curr_collaps_index--;
            }

            core().is_animating =true ;
        }
        if(key == GLFW_KEY_1){
            max+=10;
        }
        if (key == GLFW_KEY_0){
            reset();
        }

    }
}


void BasicScene::reset(){
    done_collaps.resize(max);
    prev_V.resize(max);
    prev_F.resize(max);


    done_collaps[0] = true;
    V = IV;
    F = IF;

    prev_V[0] = V;
    prev_F[0] = F;

    igl::edge_flaps(F,E,EMAP,EF,EI);
    C.resize(E  .rows(),V.cols());
    VectorXd costs(E.rows());

    EQ = Eigen::VectorXi::Zero(E.rows());
    num_collapsed = 0;
    curr_collaps_index = 0;


    set_mesh(V,F);
    data().set_face_based(true);

    //add this functionality
    Qit.resize(E.rows()); //number of edges
    caculateQMatrix();
    Qt.clear();
    for (int j = 0; j < E.rows(); j++) {
        caculateCostAndPlacment(j);
    }
}


void BasicScene::Simp(){

    // If animating then collapse 10% of edges
    if(core().is_animating && !Qt.empty() && (curr_collaps_index < max && curr_collaps_index >= 0))
    {
        if(!done_collaps[curr_collaps_index]) {
            bool something_collapsed = false;
            // collapse edge
            const int max_iter = std::ceil(0.01 * Qt.size());
            for (int j = 0; j < max_iter; j++) {
                if (!collapse_edge()) {
                    break;
                }
                something_collapsed = true;
                num_collapsed++;
            }

            if (something_collapsed) {
                //data().clear();
                set_mesh(V, F);
                data().set_face_based(true);
                //std::cout<< "ver: \n" << F <<std::endl;
                core().is_animating = false;
                done_collaps[curr_collaps_index] = true;
                prev_V[curr_collaps_index] = V;
                prev_F[curr_collaps_index] = F;

            }
        } else{
            //data().clear();
            set_mesh(prev_V[curr_collaps_index],prev_F[curr_collaps_index]);
            data().set_face_based(true);
            //std::cout<< "ver: \n" << F <<std::endl;
            core().is_animating = false;
        }
    }
    std::cout<< "index isssss: " << curr_collaps_index << "\n";
}

void BasicScene::caculateQMatrix() {
    std::vector<std::vector<int> > VF;// vertex to faces
    std::vector<std::vector<int> > VFi;//not used
    int n = V.rows();
    Qmatrix.resize(n);
    vertex_triangle_adjacency(n, F, VF, VFi);

    //for viewer
    Eigen::MatrixXd F_normals = data().F_normals;

    for (int i = 0; i < n; i++) {
        //initialize
        Qmatrix[i] = Eigen::Matrix4d::Zero();

        //caculate vertex  Q matrix
        for (int j = 0; j < VF[i].size(); j++) {
            Eigen::Vector3d normal = F_normals.row(VF[i][j]).normalized();//get face normal
            // the equation is ax+by+cz+d=0
            Eigen::Matrix4d curr;
            double a = normal[0];
            double b = normal[1];
            double c = normal[2];
            double d = V.row(i) * normal;
            d *= -1;
            curr.row(0) = Eigen::Vector4d(a*a, a*b, a*c, a*d);
            curr.row(1) = Eigen::Vector4d(a*b, b*b, b*c, b*d);
            curr.row(2) = Eigen::Vector4d(a*c, b*c, c*c, c*d);
            curr.row(3) = Eigen::Vector4d(a*d,b*d, c*d, d*d);
            Qmatrix[i] += curr;
        }

    }

}
void BasicScene:: caculateCostAndPlacment(int edge){
    //vertexes of the edge
    int v1 = E(edge, 0);
    int v2 = E(edge, 1);

    Eigen::Matrix4d Qedge= Qmatrix[v1] + Qmatrix[v2];

    Eigen::Matrix4d Qposition = Qedge; //we will use this to find v` position
    Qposition.row(3) = Eigen::Vector4d(0, 0, 0, 1);
    Eigen::Vector4d vposition;
    double cost;
    bool isInversable;
    Qposition.computeInverseWithCheck(Qposition, isInversable);
    if (isInversable) {
        vposition = Qposition * (Eigen::Vector4d(0, 0, 0, 1));
        cost = vposition.transpose() * Qedge * vposition;
    }
    else {
        //find min error from v1 v2 v1+v2/2
        Eigen::Vector4d v1p;
        v1p<< V.row(v1), 1;;
        double cost1 = v1p.transpose() * Qedge * v1p;

        Eigen::Vector4d v2p;
        v1p << V.row(v2), 1;;
        double cost2 = v2p.transpose() * Qedge * v2p;

        Eigen::Vector4d v12p;
        v1p << ((V.row(v1)+ V.row(v2))/2), 1;;
        double cost3 = v12p.transpose() * Qedge * v12p;
        if (cost1 < cost2 && cost1 < cost3) {
            vposition = v1p;
            cost = cost1;
        }
        else if (cost2 < cost1 && cost2 < cost3) {
            vposition = v2p;
            cost = cost2;
        }
        else {
            vposition = v12p;
            cost = cost3;

        }
    }
    Eigen::Vector3d pos;
    pos[0] = vposition[0];
    pos[1] = vposition[1];
    pos[2] = vposition[2];
    C.row(edge) = pos;
    Qit[edge] = Qt.insert(std::pair<double, int>(cost, edge)).first;
}
bool BasicScene::collapse_edge(){
    PriorityQueue&  curr_Q = Qt;
    std::vector<PriorityQueue::iterator >& curr_Qit = Qit;
    int e1, e2, f1, f2; //be used in the igl collapse_edge function
    if (curr_Q.empty())
    {
        // no edges to collapse
        return false;
    }
    std::pair<double, int> pair = *(curr_Q.begin());
    if (pair.first == std::numeric_limits<double>::infinity())
    {
        // min cost edge is infinite cost
        return false;
    }
    curr_Q.erase(curr_Q.begin()); //delete from the queue
    int e = pair.second; //the lowest cost edge in the queue
    //the 2 vertix of the edge
    int v1 = E.row(e)[0];
    int v2 = E.row(e)[1];

    curr_Qit[e] = curr_Q.end();

    //get the  list of faces around the end point the edge
    std::vector<int> N = igl::circulation(e, true, EMAP, EF, EI);
    std::vector<int> Nd = igl::circulation(e, false, EMAP, EF, EI);
    N.insert(N.begin(), Nd.begin(), Nd.end());

    //collapse the edage
    bool is_collapsed = igl::collapse_edge(e, C.row(e), V, F, E, EMAP, EF, EI, e1, e2, f1, f2);
    if(is_collapsed){


        // Erase the two, other collapsed edges
        curr_Q.erase(curr_Qit[e1]);
        curr_Qit[e1] = curr_Q.end();
        curr_Q.erase(curr_Qit[e2]);
        curr_Qit[e2] = curr_Q.end();

        //update the Q matrix for the 2 veterixes we collapsed
        Qmatrix[v1] = Qmatrix[v1] + Qmatrix[v2];
        Qmatrix[v2] = Qmatrix[v1] + Qmatrix[v2];

        Eigen::VectorXd newPosition;
        // update local neighbors
        // loop over original face neighbors
        for (auto n : N)
        {
            if (F(n, 0) != IGL_COLLAPSE_EDGE_NULL ||
                F(n, 1) != IGL_COLLAPSE_EDGE_NULL ||
                F(n, 2) != IGL_COLLAPSE_EDGE_NULL)
            {
                for (int v = 0; v < 3; v++)
                {
                    // get edge id
                    const  int ei = EMAP(v * F.rows() + n);
                    // erase old entry
                    curr_Q.erase(curr_Qit[ei]);
                    // compute cost and potential placement and place in queue
                    caculateCostAndPlacment(ei);
                    newPosition = C.row(ei);
                }
            }
        }
        std::cout << "edge " << e << ",cost " << pair.first << ",new position (" << newPosition[0] << ","
                  << newPosition[1] << "," << newPosition[2] << ")" << std::endl;
    }
    else
    {
        // reinsert with infinite weight (the provided cost function must **not**
        // have given this un-collapsable edge inf cost already)
        pair.first = std::numeric_limits<double>::infinity();
        curr_Qit[e] = curr_Q.insert(pair).first;
    }
    return is_collapsed;

}
void BasicScene::set_mesh(Eigen::MatrixXd &OV,Eigen::MatrixXi &OF) {
    data().set_mesh(OV,OF);
    Mesh *temp = new Mesh(curr_mesh_name,OV,OF,data().V_normals,curr_mesh->GetMeshList()[0]->data[0].textureCoords);
    shared_ptr<Mesh> m = (std::make_shared<Mesh>(*temp));
    std::vector <shared_ptr<Mesh >> new_mesh_list;
    new_mesh_list.resize(1);
    new_mesh_list[0] = m;
    curr_mesh->SetMeshList(new_mesh_list);
}
void BasicScene::Simp2(){
    Eigen::MatrixXi FI = curr_mesh->GetMeshList()[0]->data[0].faces;
    Eigen::MatrixXd VI = curr_mesh->GetMeshList()[0]->data[0].vertices;
    Mesh *temp1 = new Mesh(curr_mesh_name,VI,FI,curr_mesh->GetMeshList()[0]->data[0].vertexNormals,curr_mesh->GetMeshList()[0]->data[0].textureCoords);
    shared_ptr<Mesh> m1 = (std::make_shared<Mesh>(*temp1));
    std::vector <shared_ptr<Mesh >> new_mesh_list;
    new_mesh_list.resize(max);
    new_mesh_list[0] = m1;

    for(int i = 1; i <max; i++) {
        // If animating then collapse 10% of edges
        if (!Qt.empty() && (curr_collaps_index < max && curr_collaps_index >= 0)) {

            bool something_collapsed = false;
            // collapse edge
            const int max_iter = std::ceil(0.01 * Qt.size());
            for (int j = 0; j < max_iter; j++) {
                if (!collapse_edge()) {
                    break;
                }
                something_collapsed = true;
                num_collapsed++;
            }

            if (something_collapsed) {
                //data().clear();
                FI = F;
                VI = V;

                set_mesh(VI, FI);
                Mesh *temp = new Mesh(curr_mesh_name, V, F, data().V_normals,
                                      curr_mesh->GetMeshList()[0]->data[0].textureCoords);
                shared_ptr<Mesh> m = (std::make_shared<Mesh>(*temp));
                new_mesh_list[i] = m;
            }

        }
    }
    curr_mesh->SetMeshList(new_mesh_list);
    std::cout<< "index isssss: " << curr_collaps_index << "\n";
}

