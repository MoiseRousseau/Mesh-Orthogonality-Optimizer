#ifndef CONNECTION_H
#define CONNECTION_H

#include <vector>
#include "Point.h"
#include "Vertice.h"
#include "Element.h"

class Connection
{
    public:

        Element* element_id_up = nullptr; //element in the direction of face normal
        Element* element_id_dn = nullptr; //the other
        std::vector<Vertice*> vertices;
        Vertice* vertice_up = nullptr;
        Vertice* vertice_dn = nullptr;
        unsigned int type = 0; //face type: triangle = 3 or quad = 4
        double area = -1; //face area
        Point normal; //face normal, norm = area
        Point cell_center_vector;  //vector linking the two cell center
        double cell_center_vector_norm = 0.; // and its norm
        double orthogonality = -1; //(r_f*n_f), 1 mean no error, 0 full error
        double skewness = -1; //TODO
        //double weight = 1; //weighting factor on this connection

        Connection() {}
        virtual ~Connection(void) {};


        void set_vertices(std::vector<Vertice*> vec_p) {
            for (auto v : vec_p) {
                vertices.push_back(v);
            }
            type = vertices.size();
        }
        void set_up_info(Element* id_up, Vertice* v_up) {
            element_id_up = id_up;
            vertice_up = v_up;
        }
        void set_dn_info(Element* id_dn, Vertice* v_dn) {
            element_id_dn = id_dn;
            vertice_dn = v_dn;
        }
        void check_orientation() {
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

        double compute_orthogonality() {
            compute_cell_center_vector();
            compute_normal();
            orthogonality = cell_center_vector.dot(normal);
            if (std::abs(orthogonality) < 1e-6) {orthogonality = 1e-6;}
            //if (orthogonality > 1) {orthogonality = 1;}
            return orthogonality;
        }
        double compute_skewness() {
            skewness = -1; //TODO
            return skewness;
        }
        void compute_cell_center_vector() {
            cell_center_vector = element_id_up->center() - element_id_dn->center();
            cell_center_vector_norm = cell_center_vector.norm();
            cell_center_vector /= cell_center_vector_norm;
        }
        void compute_normal() {
            Point u = *vertices[1]->coor-*vertices[0]->coor;
            Point v = *vertices[2]->coor-*vertices[1]->coor;
            normal = u.cross(v);
            area = normal.norm();
            normal /= area;
            area /= 2;
            if (type == 4) {
                Point u = *vertices[3]->coor-*vertices[0]->coor;
                Point v = *vertices[2]->coor-*vertices[3]->coor;
                area += u.cross(v).norm()/2;
            }
        }

    protected:

    private:
};

#endif // CONNECTION_H
