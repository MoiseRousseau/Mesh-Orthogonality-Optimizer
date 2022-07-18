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



void OrthOpt::computeCostFunction() {
    cost_function_value = 0;
    #pragma omp parallel for
    for (unsigned int count=0; count != mesh->n_connections_internal; count++) {
    	mesh->connections_internal[count]->compute_orthogonality();
        face_error[count] = face_weight[count] * \
                       Ef->get_value(mesh->connections_internal[count]->orthogonality);
    }
    #pragma omp parallel for reduction (+:cost_function_value)
    for (double err : face_error) {
        cost_function_value += err;
    }
}


void OrthOpt::computeCostDerivative()
{
    computeCostFunction(); //update normal, cell center vector and other
    unsigned int index;
    double prefactor;
    Point deriv;
    Connection* con;
    for (auto p=face_error_derivative.begin();
              p!=face_error_derivative.end(); p++) {
        p->x = 0.; p->y = 0.; p->z = 0.;
    }
    #pragma omp parallel for private(con, deriv, index, prefactor) shared(face_error_derivative)
    for (unsigned int count=0; count != mesh->n_connections_internal; count++) {
        con = mesh->connections_internal[count];
        prefactor =  face_weight[count] * Ef->get_derivative(con->orthogonality);
        if (con->element_id_dn->type + con->element_id_up->type == 8) { 
            //two tet case a
            for (Vertice* p : con->vertices) { //A position
                if (p->fixed) {continue;}
                index = derivative_vertice_ids[p->natural_id-1];
                deriv = derivative_A_position(con, p) * prefactor;
                #pragma omp atomic update
                face_error_derivative[index].x += deriv.x;
                #pragma omp atomic update
                face_error_derivative[index].y += deriv.y;
                #pragma omp atomic update
                face_error_derivative[index].z += deriv.z;
            }
            deriv = derivative_E_position(con) * prefactor;
            //E position
            if (con->vertice_dn->fixed == false) {
                index = derivative_vertice_ids[con->vertice_dn->natural_id-1];
                #pragma omp atomic update
                face_error_derivative[index].x -= deriv.x;
                #pragma omp atomic update
                face_error_derivative[index].y -= deriv.y;
                #pragma omp atomic update
                face_error_derivative[index].z -= deriv.z;
            }
            //F position
            if (con->vertice_up->fixed == false) {
                index = derivative_vertice_ids[con->vertice_up->natural_id-1];
                #pragma omp atomic update
                face_error_derivative[index].x += deriv.x;
                #pragma omp atomic update
                face_error_derivative[index].y += deriv.y;
                #pragma omp atomic update
                face_error_derivative[index].z += deriv.z;
            }
        }
        if (con->element_id_dn->type + con->element_id_up->type == 10 and
            con->type == 3) { 
            //triangle between two pyr
            //TODO
        }
        if (con->element_id_dn->type + con->element_id_up->type == 10 and
            con->type == 4) { 
            //quad between two pyr
            //TODO
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

void OrthOpt::save_face_error_derivative(std::string f_out) {
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
}

