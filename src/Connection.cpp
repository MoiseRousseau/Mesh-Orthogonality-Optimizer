#include "Connection.h"
#include <Eigen/Core>
#include <iostream>

void Connection::check_orientation() {
    compute_orthogonality();
    if (orthogonality < 0.) {
        auto temp = element_id_up;
        element_id_up = element_id_dn;
        element_id_dn = temp;
        //auto temp2 = vertice_up;
        //vertice_up = vertice_dn;
        //vertice_dn = temp2;
        orthogonality = -orthogonality;
    }
}

double Connection::compute_orthogonality() {
    compute_cell_center_vector();
    compute_normal();
    orthogonality = cell_center_vector.dot(normal);
    if (std::abs(orthogonality) < 1e-6) {orthogonality = 1e-6;}
    //if (orthogonality > 1) {orthogonality = 1;}
    return orthogonality;
}

void Connection::compute_cell_center_vector() {
    cell_center_vector << element_id_up->center() - element_id_dn->center();
    cell_center_vector_norm = cell_center_vector.norm();
    cell_center_vector /= cell_center_vector_norm;
}

void Connection::compute_normal() {
    if (vertices.size() == 2) { //2D case
        normal[0] = (*vertices[1]->coor)[1] - (*vertices[0]->coor)[1];
        normal[1] = - (*vertices[1]->coor)[0] + (*vertices[0]->coor)[0];
        area = normal.norm(); //length in this case
        normal /= area;
    }
    else { //3D case
        Eigen::VectorXd u(3);
        u << (*vertices[1]->coor) - (*vertices[0]->coor);
        Eigen::VectorXd v(3);
        v << (*vertices[2]->coor) - (*vertices[1]->coor);
        //normal = u.cross(v);
        normal << u[1]*v[2]-u[2]*v[1], u[0]*v[2]-u[2]*v[0], u[0]*v[1]-u[1]*v[0];
        area = normal.norm();
        normal /= area;
        area /= 2;
        if (vertices.size() == 4) {
            Eigen::Vector3d u = *(vertices[3]->coor)-*(vertices[0]->coor);
            Eigen::Vector3d v = *(vertices[2]->coor)-*(vertices[3]->coor);
            area += u.cross(v).norm()/2;
        }
    }
}

Eigen::MatrixXd Connection::derivative_A_position(Vertice* A) {
    Eigen::MatrixXd deriv;
    if (vertices.size() == 2) { //2D case
        deriv = Eigen::MatrixXd::Zero(2,2);
        double x = 1./area; //area is length
        deriv << 0, -x, x, 0;
    }
    else {
        //determine B and C point in various cases
        Vertice* B = nullptr;
        Vertice* C = nullptr;
        unsigned int index = 0;
        for (Vertice* v : vertices) {
            if (v == A) break;
            index += 1;
        }
        if (index == 0) {B = vertices[1]; C = vertices[2];}
        else if (index == 1) {
            B = vertices[2];
            if (vertices.size() == 3) {C = vertices[0];}
            else {C = vertices[3];}
        }
        else if (index == 2) {
            if (vertices.size() == 3) {B = vertices[0]; C = vertices[1];}
            else {B = vertices[3]; C = vertices[0];}
        }
        else {B = vertices[0]; C = vertices[1];}
        //compute normal derivative
        Eigen::Vector3d temp, BC;
        BC = *(C->coor) - *(B->coor);
        if (vertices.size() == 3) BC *= 0.5 / area; //face is a triangle
        else BC *= 1. / area;
        //build matrix
        deriv = Eigen::MatrixXd::Zero(3,3);
        deriv << 0, -BC[2], BC[1],
                 BC[2], 0, -BC[0],
                 -BC[1], BC[0], 0;
    }
    return deriv;
}
