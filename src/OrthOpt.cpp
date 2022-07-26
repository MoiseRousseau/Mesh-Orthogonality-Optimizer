#include "OrthOpt.h"
#include <cmath>
#include <iostream>
#include <vector>
#include <iostream>
#include <fstream>
#include <omp.h>

#include "Mesh.h"
#include "Point.h"
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
    Point deriv, temp;
    Connection* con;
    for (size_t i=0; i!=grad.size(); i++) {
        grad[i] = 0;
    }
    #pragma omp parallel for private(con, deriv, index, prefactor) shared(grad)
    for (unsigned int count=0; count != mesh->n_connections_internal; count++) {
        con = mesh->connections_internal[count];
        prefactor =  face_weight[count] * Ef->get_derivative(con->orthogonality);
        temp = con->normal - con->cell_center_vector * (con->orthogonality);
        temp /= con->cell_center_vector_norm;
        for (Vertice* p : con->vertices) {
            //normal derivative
            deriv = derivative_A_position(con, p);
            if (con->element_id_dn->type + con->element_id_up->type != 8 or
                con->element_id_dn->type + con->element_id_up->type != -6) {
                deriv += temp * (con->element_id_up->center_derivative(p) - con->element_id_dn->center_derivative(p) )
            }
            index = mesh->dim * derivative_vertice_ids[p->natural_id-1];
            #pragma omp atomic update
            grad[index] -= prefactor * deriv.x;
            #pragma omp atomic update
            grad[index+1] -= prefactor * deriv.y;
            #pragma omp atomic update
            grad[index+2] -= prefactor * deriv.z;
        }
        // add contribution of vertices of element dn
        for (Vertice* p : con->element_id_dn->vertices) {
            if (p->fixed) {continue;}
            if (con->vertices->contains(p)) {
                continue; //already treated
            }
            deriv = temp * con->element_id_up->center_derivative(p));
            index = mesh->dim * derivative_vertice_ids[p->natural_id-1];
            #pragma omp atomic update
            grad[index] -= prefactor * deriv.x;
            #pragma omp atomic update
            grad[index+1] -= prefactor * deriv.y;
            #pragma omp atomic update
            grad[index+2] -= prefactor * deriv.z;
        }
        // add contribution of vertices of element up
        for (Vertice* p : con->element_id_up->vertices) { 
            if (p->fixed) {continue;}
            if (con->vertices->contains(p)) {
                continue; //already treated
            }
            deriv = temp * con->element_id_up->center_derivative(p));
            index = mesh->dim * derivative_vertice_ids[p->natural_id-1];
            #pragma omp atomic update
            grad[index] -= prefactor * deriv.x;
            #pragma omp atomic update
            grad[index+1] -= prefactor * deriv.y;
            #pragma omp atomic update
            grad[index+2] -= prefactor * deriv.z;
        }
    }
}


void OrthOpt::update_vertices_position(const Eigen::VectorXd &x) {
    unsigned int index = 0;
    for (Vertice* v : mesh->vertices) {
        if (v->fixed == false) {
            //vertice is not fixed, update position
            v->coor->x = x[index];
            v->coor->y = x[index+1];
            if (mesh->dim == 3) v->coor->z = x[index+2];
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

