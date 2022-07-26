#include "OrthOpt.h"
#include <cmath>
#include <iostream>
#include <vector>
#include <iostream>
#include <fstream>
#include <omp.h>

#include "Mesh.h"
#include "Connection.h"
#include <Eigen/Core>

typedef std::array<unsigned int, 3> face_id;
typedef std::map<face_id, Connection> map_faces;


void OrthOpt::computeCostFunction() {
    cost_function_value = 0;
    #pragma omp parallel for reduction (+:cost_function_value)
    for (unsigned int count=0; count != mesh->n_connections_internal; count++) {
    	mesh->connections_internal[count]->compute_orthogonality();
        cost_function_value += face_weight[count] * \
                       Ef->get_value(mesh->connections_internal[count]->orthogonality);
    }
}


void OrthOpt::computeCostDerivative(Eigen::VectorXd& grad)
{
    computeCostFunction(); //update normal, cell center vector and other
    unsigned int index;
    double prefactor;
    Eigen::Vector3d deriv;
    Connection* con;
    for (size_t i=0; i!=grad.size(); i++) {
        grad[i] = 0;
    }
    #pragma omp parallel for private(con, deriv, index, prefactor) shared(grad)
    for (unsigned int count=0; count != mesh->n_connections_internal; count++) {
        con = mesh->connections_internal[count];
        prefactor =  face_weight[count] * Ef->get_derivative(con->orthogonality);
        /*if (con->element_id_dn->type + con->element_id_up->type == -8) { 
            //2D case, segment between two quad
        }
        else if (con->element_id_dn->type + con->element_id_up->type == -6) { 
            //2D case, segment between two triangle
        }*/
        if (con->element_id_dn->type + con->element_id_up->type == 8) { 
            //two tet case a
            for (Vertice* p : con->vertices) { //A position
                if (p->fixed) {continue;}
                index = mesh->dim * derivative_vertice_ids[p->natural_id-1];
                deriv = derivative_A_position(con, p) * prefactor;
                for (size_t i=0; i<mesh->dim; i++) {
                    #pragma omp atomic update
                    grad[index+i] += deriv[i];
                }
            }
            deriv = derivative_E_position(con) * prefactor;
            //E position
            if (con->vertice_dn->fixed == false) {
                index = mesh->dim * derivative_vertice_ids[con->vertice_dn->natural_id-1];
                for (size_t i=0; i<mesh->dim; i++) {
                    #pragma omp atomic update
                    grad[index+i] += deriv[i];
                }
            }
            //F position
            if (con->vertice_up->fixed == false) {
                index = mesh->dim * derivative_vertice_ids[con->vertice_up->natural_id-1];
                for (size_t i=0; i<mesh->dim; i++) {
                    #pragma omp atomic update
                    grad[index+i] += deriv[i];
                }
            }
        }
        else if (con->element_id_dn->type + con->element_id_up->type == 10 and
            con->vertices.size() == 3) { 
            //triangle between two pyr
            //TODO
        }
        else if (con->element_id_dn->type + con->element_id_up->type == 10 and
            con->vertices.size() == 4) { 
            //quad between two pyr
            //TODO
        }
    }
}


void OrthOpt::update_vertices_position(const Eigen::VectorXd &x) {
    unsigned int index = 0;
    for (Vertice* v : mesh->vertices) {
        if (v->fixed == false) {
            //vertice is not fixed, update position
            for (size_t i=0; i<mesh->dim; i++) {
                (*(v->coor))[i] += x[index];
            }
            index += mesh->dim;
        }
    }
}


void OrthOpt::save_face_error(std::string f_out) {
    std::ofstream out(f_out);
    for (Connection* con : mesh->connections_internal) {
        out << con->element_id_dn->natural_id << " ";
        out << con->element_id_up->natural_id << " ";
        //out << std::acos(std::abs(con->orthogonality)) * 57.29583 << std::endl; //57.2583 = 180 / pi
        out << 1-con->orthogonality << std::endl;
    }
    out.close();
}

/*void OrthOpt::save_face_error_derivative(std::string f_out) {
    std::ofstream out(f_out);
    unsigned int index;
    for (Vertice* v : mesh->vertices) {
        if (v->fixed) {continue;}
        index = derivative_vertice_ids[v->natural_id-1];
        out << v->natural_id << " ";
        out << face_error_derivative[index-1].x << " ";
        out << face_error_derivative[index-1].y << " ";
        out << face_error_derivative[index-1].z << std::endl;
    }
    out.close();
}*/

/*
void OrthOpt::decompose_mesh() {
    std::set<face_id> faces_set;
    for (Element* elem : mesh->elements) {
        if (elem->type == 4) { //tetrahedron
            //face 0 1 2
            face_id key = build_face_key(elem, 0, 1, 2);
            //face 1 2 3
            face_id key = build_face_key(elem, 1, 2, 3);
            //face 0 1 3, opposite 2
            face_id key = build_face_key(elem, 0, 1, 3);
            //face 0 2 3, opposite 1
            face_id key = build_face_key(elem, 0, 2, 3);
        }
        else if (elem->type == 5) { //pyramid
            face_id key = build_face_key(elem, 0, 1, 2, 3);
            //other are treated by tet above
        }
        else if (elem->type == 6) { //prisms
            //not optimized, so do not figure
            //TODO: I must fix the point in this case, else, they would move and degenerate the prisms. Verify
        }
        else if (elem->type == 8) { //hex
            //not optimized, so do not figure
        }
        else {
            std::cerr << std::endl;
            std::cerr << "Mesh decompose error" << std::endl;
            std::cerr << "element type not recognized: " << elem->type;
            std::cerr << std::endl;
            std::cerr << std::endl;
            exit(1);
        }
        faces_set.insert(key);
    }
    for (auto x : faces_set) {
        Connection* face = new Connection();
        faces.push_back(face);
    }
}

face_id OrthOpt::build_face_key(Element* elem, unsigned int i, unsigned int j, unsigned int k, unsigned int h) {
    //build key
    face_id key;
    key[0] = elem->vertice_ids[i];
    key[1] = elem->vertice_ids[j];
    key[2] = elem->vertice_ids[k];
    std::sort(key.begin(),key.end());
    if (h>0) {
        if (h > key[2]) {
            key[2} = elem->vertice_ids[h];
            std::sort(key.begin(),key.end());
        }
    }
    return key
}*/
