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


void OrthOpt::computeCostFunction() {
    cost_function_value = 0;
    #pragma omp parallel for reduction (+:cost_function_value)
    for (unsigned int count=0; count != mesh->n_connections_internal; count++) {
    	mesh->connections_internal[count]->compute_orthogonality();
        cost_function_value += face_weight[count] * \
                       Ef->get_value(mesh->connections_internal[count]->orthogonality);
    }
    #pragma omp parallel for reduction (+:cost_function_value)
    for (auto elem : mesh->elements) {
        if (elem->is_inverted()) cost_function_value = std::nan("");
    }
}


void OrthOpt::computeCostDerivative(Eigen::VectorXd& grad)
{
    computeCostFunction(); //update normal, cell center vector and other
    for (size_t i=0; i!=grad.size(); i++) {
        grad[i] = 0;
    }
    #pragma omp parallel for shared(grad)
    for (unsigned int count=0; count != mesh->n_connections_internal; count++) {
        Connection* con = mesh->connections_internal[count];
        double prefactor =  face_weight[count] * Ef->get_derivative(con->orthogonality); //w_f E_f'
        // Configuration I, point not on the face
        // if this configuration, only the cell center derivative is non-null
        // temp = (n_f^T - (n_f^T . r_f) r_f^T) / | R_f | (comes from unit vector derivative)
        Eigen::VectorXd temp = con->normal - con->cell_center_vector * (con->orthogonality);
        temp /= con->cell_center_vector_norm;
        Eigen::RowVectorXd deriv = Eigen::RowVectorXd::Zero(mesh->dim);
        // add contribution of vertices of element dn
        for (Vertice* p : con->element_id_dn->vertices) {
            if (std::find(con->vertices.begin(), con->vertices.end(), p) != con->vertices.end()) continue;
            if (p->fixed) {continue;}
            deriv = temp.transpose() * con->element_id_up->center_derivative(p);
            size_t index = mesh->dim * derivative_vertice_ids[p->natural_id-1];
            for (size_t i=0; i<mesh->dim; i++) {
                #pragma omp atomic update
                grad[index+i] -= prefactor * deriv[i];
            }
        }
        // add contribution of vertices of element up
        for (Vertice* p : con->element_id_up->vertices) { 
            if (std::find(con->vertices.begin(), con->vertices.end(), p) != con->vertices.end()) continue;
            if (p->fixed) {continue;}
            deriv = temp.transpose() * con->element_id_up->center_derivative(p);
            size_t index = mesh->dim * derivative_vertice_ids[p->natural_id-1];
            for (size_t i=0; i<mesh->dim; i++) {
                #pragma omp atomic update
                grad[index+i] -= prefactor * deriv[i];
            }
        }
        
        // configuration II, point on the face
        // 1. contribution of the normal part
        deriv = Eigen::RowVectorXd::Zero(mesh->dim);
        temp = con->cell_center_vector - con->normal * (con->orthogonality);
        temp /= con->cell_center_vector_norm; 
        //we do not divide by the area because it's done in the derivative_normal routine
        for (Vertice* p : con->vertices) {
            if (p->fixed) {continue;}
            //normal derivative
            Eigen::MatrixXd deriv_dA = con->derivative_A_position_normal(p) * temp;
            for (int i=0; i<mesh->dim; i++) deriv[i] = deriv_dA(i,0);
            size_t index = mesh->dim * derivative_vertice_ids[p->natural_id-1];
            for (size_t i=0; i<mesh->dim; i++) {
                #pragma omp atomic update
                grad[index+i] -= prefactor * deriv[i];
            }
        }
        
        // 2. contribution of the cell center part
        // only if there is different element type
        if (con->element_id_dn->type == -3 and con->element_id_up->type != -3 or con->element_id_dn->type == 4 and con->element_id_up->type != 4) {
            deriv = Eigen::RowVectorXd::Zero(mesh->dim);
            temp = con->normal - con->cell_center_vector * (con->orthogonality);
            temp /= con->cell_center_vector_norm;
            for (Vertice* p : con->vertices) {
                if (p->fixed) {continue;}
                //std::cout << temp << std::endl;
                //std::cout << con->element_id_up->center_derivative(p) << std::endl;
                //std::cout << con->element_id_dn->center_derivative(p) << std::endl;
                deriv = temp.transpose() * (con->element_id_up->center_derivative(p) - con->element_id_dn->center_derivative(p)); 
                //std::cout << deriv << std::endl;
                size_t index = mesh->dim * derivative_vertice_ids[p->natural_id-1];
                for (size_t i=0; i<mesh->dim; i++) {
                    #pragma omp atomic update
                    grad[index+i] -= prefactor * deriv[i];
                }
            }
        }
    }
}


/*
Compute the derivative of the cost function using finite difference
For debug purpose
*/
void OrthOpt::computeCostDerivative_FD(Eigen::VectorXd& grad, double pertub)
{
    for (size_t i=0; i!=grad.size(); i++) {
        grad[i] = 0;
    }
    for (Vertice* v : mesh->vertices) {
         if (v->fixed) {continue;}
         computeCostFunction();
         double ref_cost = cost_function_value;
         size_t index = mesh->dim * derivative_vertice_ids[v->natural_id-1];
         for (size_t i=0; i<mesh->dim; i++) {
             double coor_ini = (*v->coor)[i];
             (*v->coor)[i] *= (1+pertub);
             computeCostFunction();
             double opt_pertub = cost_function_value;
             grad[index+i] = (opt_pertub - ref_cost) / pertub;
             (*v->coor)[i] = coor_ini;
         }
         computeCostFunction();
    }
}


void OrthOpt::update_vertices_position(const Eigen::VectorXd &x) {
    unsigned int index = 0;
    for (Vertice* v : mesh->vertices) {
        if (v->fixed == false) {
            //vertice is not fixed, update position
            for (size_t i=0; i<mesh->dim; i++) {
                (*v->coor)[i] = x[index+i];
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

