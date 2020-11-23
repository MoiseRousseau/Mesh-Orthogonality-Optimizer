#include "OrthOpt.h"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <vector>
#include <map>
#include <array>
#include <iostream>
#include <fstream>

#include "Mesh.h"
#include "Point.h"
#include "Connection.h"
#include <Eigen/Core>


void OrthOpt::build_connection(std::array<unsigned int, 4> key, \
                               Element* elem_key, unsigned int opposite, \
                               std::map<std::array<unsigned int, 4>, Connection*> &unique_id_map) {

    std::map<std::array<unsigned int, 4>, Connection*>::iterator it;
    Connection* con;

    //if (key[0] == 0) {} //triangle
    //find connection in map
    it = unique_id_map.find(key);
    if (it != unique_id_map.end()) { //key already in
        //terminate connection connection
        it->second->element_id_up = elem_key;
        it->second->vertice_up = mesh->vertices[opposite-1];
        unique_id_map.erase(it);
        it->second->check_orientation();
        n_connections += 1;
    }
    else {
        con = new Connection();
        con->vertices.push_back(mesh->vertices[key[1]-1]);
        con->vertices.push_back(mesh->vertices[key[2]-1]);
        con->vertices.push_back(mesh->vertices[key[3]-1]);
        con->element_id_dn = elem_key;
        con->vertice_dn = mesh->vertices[opposite-1];
        con->type = 3;
        unique_id_map[key] = con;
        connections.push_back(con);
    }
}


void OrthOpt::populate_connections() {
    //first populate elem with pointer to vertice
    //from the input we need to identify all the connection in the mesh
    //therefore, for each face, we create a connection object
    // with vertices uniquely oriented
    mesh->decompose_mesh();
    std::array<unsigned int, 4> key;
    unsigned int opposite;
    std::map<std::array<unsigned int, 4>, Connection*> unique_id_map;
    for (Element* elem : mesh->elements) {
        if (elem->type == 4) { //tetrahedron
            //first face 0 1 2, opposite 3
            key[0] = 0;
            key[1] = elem->vertice_ids[0];
            key[2] = elem->vertice_ids[1];
            key[3] = elem->vertice_ids[2];
            opposite = elem->vertice_ids[3];
            std::sort(key.begin(),key.end());
            build_connection(key, elem, opposite, unique_id_map);

            //first face 1 2 3, opposite 0
            key[0] = 0;
            key[1] = elem->vertice_ids[1];
            key[2] = elem->vertice_ids[2];
            key[3] = elem->vertice_ids[3];
            opposite = elem->vertice_ids[0];
            std::sort(key.begin(),key.end());
            build_connection(key, elem, opposite, unique_id_map);

            //first face 0 1 3, opposite 2
            key[0] = 0;
            key[1] = elem->vertice_ids[0];
            key[2] = elem->vertice_ids[1];
            key[3] = elem->vertice_ids[3];
            opposite = elem->vertice_ids[2];
            std::sort(key.begin(),key.end());
            build_connection(key, elem, opposite, unique_id_map);

            //first face 0 2 3, opposite 1
            key[0] = 0;
            key[1] = elem->vertice_ids[0];
            key[2] = elem->vertice_ids[2];
            key[3] = elem->vertice_ids[3];
            opposite = elem->vertice_ids[1];
            std::sort(key.begin(),key.end());
            build_connection(key, elem, opposite, unique_id_map);
        }
        else if (elem->type == 5) { //pyramid
            // TODO: mesh decompose other element
        }
        else if (elem->type == 6) { //prisms

        }
        else if (elem->type == 8) { //hex

        }
        else {
            exit(1);
        }
    }
    //for (auto con = connections.begin(); con != connections.end(); con++) {
    for (auto con = connections.begin(); con != connections.end();) {
        if ((*con)->vertice_up == nullptr) {
            for (Vertice* v : (*con)->vertices) {
                v->fix_vertice();
                //TODO undo
            }
            connections.erase(con);
            boundary_connections.push_back(*con);
        }
        else {++con;}
    }
    for (Vertice* v : mesh->vertices) {
        if (v->fixed) {n_vertices_to_opt -= 1;}
    }
}



void OrthOpt::computeCostFunction()
{
    unsigned int count = 0;
    for (Connection* icon : connections) {
        if (icon->element_id_dn == nullptr) {continue;} //boundary connection
        face_error[count] = icon->weight * pow(icon->compute_error(), penalizing_power);
        count += 1;
    }
    cost_function_value = 0.;
    for (double err : face_error) {cost_function_value += err;}
}


void OrthOpt::computeCostDerivative()
{
    computeCostFunction();
    unsigned int index;
    Point deriv;
    for (index=0; index != face_error_derivative.size(); index++) {
        face_error_derivative[index] = 0.;
    }
    for (Connection* con : connections) {
        if (con->element_id_dn->type == 4 and
            con->element_id_up->type == 4) { //two tet case 1
            for (Vertice* p : con->vertices) { //A position
                if (p->fixed) {continue;}
                index = derivative_vertice_ids[p->natural_id-1];
                face_error_derivative[index-1] += derivative_A_position(con, p);
            }
            deriv = derivative_E_position(con);
            //E position
            if (con->vertice_dn->fixed == false) {
                index = derivative_vertice_ids[con->vertice_dn->natural_id-1];
                face_error_derivative[index-1] += -deriv;
            }
            //F position
            if (con->vertice_up->fixed == false) {
                index = derivative_vertice_ids[con->vertice_up->natural_id-1];
                face_error_derivative[index-1] += deriv;
            }
        }
    }
}

void OrthOpt::update_vertices_position(const Eigen::VectorXd &x) {
    unsigned int count=0;
    unsigned int index = 0;
    for (Vertice* v : mesh->vertices) {
        if (v->fixed == false) {
            //vertice is not fixed, update position
            v->coor->x = x[index];
            v->coor->y = x[index+1];
            v->coor->z = x[index+2];
            index += 3;
            count++;
        }
    }
}


void OrthOpt::save_face_non_orthogonality_angle(std::string f_out) {
    std::ofstream src(f_out);
    for (Connection* con : connections) {
        src << con->element_id_dn->natural_id << " ";
        src << con->element_id_up->natural_id << " ";
        src << std::acos(std::abs(1-con->error)) * 57.29583 << "\n";
    }
}
