#ifndef CONNECTION_H
#define CONNECTION_H

#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>
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
        double area = -1; //face area
        Eigen::Vector3d normal; //face normal, norm = area, or length in 2D
        Eigen::Vector3d cell_center_vector;  //vector linking the two cell center
        double cell_center_vector_norm = 0.; // and its norm
        double orthogonality = -1; //(r_f*n_f), 1 mean no error, 0 full error
        //double skewness = -1;
        //double weight = 1; //weighting factor on this connection

        Connection() {}
        virtual ~Connection(void) {};


        void set_vertices(std::vector<Vertice*> vec_p) {
            for (auto v : vec_p) {
                vertices.push_back(v);
            }
            //type = vertices.size();
        }
        void set_up_info(Element* id_up, Vertice* v_up) {
            element_id_up = id_up;
            vertice_up = v_up;
        }
        void set_dn_info(Element* id_dn, Vertice* v_dn) {
            element_id_dn = id_dn;
            vertice_dn = v_dn;
        }
        
        void check_orientation();
        double compute_orthogonality();
        //double compute_skewness();
        void compute_cell_center_vector();
        void compute_normal();
        
        friend std::ostream& operator<<(std::ostream& os, const Connection& con) {
            //to print connection information (debug purpose)
            os << "vertices: ";
            for (auto v : con.vertices) os << v->natural_id << ' ';
            os << std::endl;
            if (con.element_id_up != nullptr) os << "element_id_up: " << con.element_id_up->natural_id << std::endl;
            if (con.element_id_dn != nullptr) os << "element_id_dn: " << con.element_id_dn->natural_id << std::endl;
            if (con.vertice_up != nullptr) os << "vertice_up: " << con.vertice_up->natural_id << std::endl;
            if (con.vertice_dn != nullptr) os << "vertice_dn: " << con.vertice_dn->natural_id << std::endl;
            os << "area: " << con.area << std::endl;
            os << "normal: " << con.normal << std::endl;
            os << "cell_center_vector: " << con.cell_center_vector << std::endl;
            os << "orthogonality: " << con.orthogonality << std::endl;
            return os;
        }

    protected:

    private:
};

#endif // CONNECTION_H
