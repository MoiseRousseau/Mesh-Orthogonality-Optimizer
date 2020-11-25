#ifndef ORTHOPT_H
#define ORTHOPT_H

#include <array>
#include "Mesh.h"
#include "Point.h"
#include "Connection.h"
#include <Eigen/Core>
#include <map>


class OrthOpt
{
    public:
        //TODO constructor when mesh is passed by numpy
        OrthOpt() {};

        OrthOpt(Mesh* m) {
            mesh = m;
            n_vertices_to_opt = mesh->n_vertices; //initiate to n_vertices and substract the fixed point);
            populate_connections();
            face_error.resize(n_connections);
            face_error_derivative.resize(n_vertices_to_opt);
            unsigned int index = 0;
            for (Vertice* v : mesh->vertices) {
                if (v->fixed == false) {
                    index += 1;
                    derivative_vertice_ids.push_back(index);
                }
                else {derivative_vertice_ids.push_back(0);}
            }
        }

        virtual ~OrthOpt() {
            for (auto con = connections.begin(); con != connections.end(); con++) {
                delete *con;
            }
            connections.clear();
        };

        // INPUT //
        Mesh* mesh;
        std::vector<Connection*> connections;
        std::vector<Connection*> boundary_connections;
        double penalizing_power = 1.;
        unsigned int n_vertices_to_opt = 0;
        unsigned int n_connections = 0;

        // OUTPUT //
        double cost_function_value;
        std::vector<double> face_error;
        std::vector<Point> face_error_derivative;
        //for each point, store the position of its derivative in face_error_derivative
        // 0 if fixed
        std::vector<unsigned int> derivative_vertice_ids;


        void set_penalizing_power(double x) {penalizing_power = x;}
        double getCostFunctionValue() {
            computeCostFunction();
            return cost_function_value;
        }
        void weight_by_area() {
            if (connections[0]->area < 0) {computeCostFunction();}
            for (Connection* con : connections) {
                con->weight = con->area;
            }
        }
        void weight_by_volume(int method) {};

        std::vector<double> getFaceError() {
            computeCostFunction();
            return face_error;
        }
        std::vector<Point> getCostFunctionDerivative() {
            computeCostDerivative();
            return face_error_derivative;
        }

        void computeCostFunction();
        void computeCostDerivative();

        //update vertices position with eigen vector
        void update_vertices_position(const Eigen::VectorXd &x);

        //output function
        void save_face_non_orthogonality_angle(std::string f_out);
        void save_face_error_derivative(std::string f_out);


    protected:

    private:
        void build_connection(std::array<unsigned int, 4> key, \
                               Element* elem_key, unsigned int opposite, \
                               std::map<std::array<unsigned int, 4>, Connection*> &unique_id_map);
        void populate_connections();


        Point derivative_E_position(Connection* con) {
            return ((con->normal - con->cell_center_vector * (1-con->error)) \
                    / con->cell_center_vector_norm * 0.25);
        }

        Point derivative_A_position(Connection* con, Vertice* A) {
            //determine B and C point in various cases
            Vertice* B = nullptr;
            Vertice* C = nullptr;
            unsigned int index = 0;
            for (Vertice* v : con->vertices) {
                if (v == A) break;
                index += 1;
            }
            if (index == 0) {B = con->vertices[1]; C = con->vertices[2];}
            else if (index == 1) {
                B = con->vertices[2];
                if (con->type == 3) {C = con->vertices[0];}
                else {C = con->vertices[3];}
            }
            else if (index == 2) {
                if (con->type == 3) {B = con->vertices[0]; C = con->vertices[1];}
                else {B = con->vertices[3]; C = con->vertices[0];}
            }
            else {B = con->vertices[0]; C = con->vertices[1];}
            //compute normal derivative
            if (con->type == 3) {
                return ((con->cell_center_vector - con->normal*(1-con->error)).cross( \
                        *(C->coor)-*(B->coor)) / (2*con->area));
            }
            else {
                return ((con->cell_center_vector - con->normal*(1-con->error)).cross( \
                        *(C->coor)-*(B->coor)) / (con->area));
            }
        }

};

#endif // ORTHOPT_H
