#include "Connection.h"
#include <Eigen/Core>


void Connection::check_orientation() {
    compute_orthogonality();
    if (orthogonality < 0.) {
        auto temp = element_id_up;
        element_id_up = element_id_dn;
        element_id_dn = temp;
        auto temp2 = vertice_up;
        vertice_up = vertice_dn;
        vertice_dn = temp2;
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
    cell_center_vector = element_id_up->center() - element_id_dn->center();
    cell_center_vector_norm = cell_center_vector.norm();
    cell_center_vector /= cell_center_vector_norm;
}

void Connection::compute_normal() {
    if (vertices.size() == 2) { //2D case
        Eigen::Vector3d u = *vertices[1]->coor-*vertices[0]->coor;
        normal[0] = u[1];
        normal[1] = -u[0];
        area = normal.norm(); //length in this case
        normal /= area;
    }
    else { //3D case
        Eigen::Vector3d u = *vertices[1]->coor-*vertices[0]->coor;
        Eigen::Vector3d v = *vertices[2]->coor-*vertices[1]->coor;
        normal = u.cross(v);
        area = normal.norm();
        normal /= area;
        area /= 2;
        if (vertices.size() == 4) {
            Eigen::Vector3d u = *vertices[3]->coor-*vertices[0]->coor;
            Eigen::Vector3d v = *vertices[2]->coor-*vertices[3]->coor;
            area += u.cross(v).norm()/2;
        }
    }
}
