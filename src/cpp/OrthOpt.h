#ifndef ORTHOPT_H
#define ORTHOPT_H

#include <array>
#include "Mesh.h"
#include "Point.h"
#include "Connection.h"
#include <Eigen/Core>



class OrthOpt
{
    public:

        // INPUT //
        Mesh* mesh;
        double penalizing_power = 1.;
        unsigned int n_vertices_to_opt = 0;

        // OUTPUT //
        double cost_function_value;
        std::vector<double> face_weight;
        std::vector<double> face_error;
        std::vector<Point> face_error_derivative;
        //for each point, store the position of its derivative in face_error_derivative
        // 0 if fixed
        std::vector<unsigned int> derivative_vertice_ids;
        
        //TODO constructor when mesh is passed by numpy
        OrthOpt() {};

        OrthOpt(Mesh* m) {
            mesh = m;
            if (mesh->n_connections_internal == 0) {mesh->decompose();};
            n_vertices_to_opt = mesh->n_vertices; //initiate to n_vertices and substract the fixed point);
            face_error.resize(mesh->n_connections_internal);
            face_weight.resize(mesh->n_connections_internal);
            derivative_vertice_ids.resize(mesh->n_vertices);
            unsigned int index = 0;
            size_t count = -1;
            for (Vertice* v : mesh->vertices) {
                count++;
                if (v->fixed == false) {
                    index += 1;
                    derivative_vertice_ids[count] = index;
                }
                else {
                    derivative_vertice_ids[count] = -1;
                    n_vertices_to_opt -= 1;
                }
            }
            face_error_derivative.resize(n_vertices_to_opt);
        }

        virtual ~OrthOpt() {};


        void set_penalizing_power(double x) {penalizing_power = x;}
        double getCostFunctionValue() {
            computeCostFunction();
            return cost_function_value;
        }
        
        void weight_unitary() {
            for (size_t count=0; count!=mesh->n_connections_internal; count++) {
                face_weight[count] = 1.;
            }
        }
        void weight_by_area() {
            size_t count = 0;
            for (Connection* con : mesh->connections_internal) {
                face_weight[count] = con->area;
                count++;
            }
        }
        void weight_by_area_inverse() {
            size_t count = 0;
            for (Connection* con : mesh->connections_internal) {
                face_weight[count] = 1./con->area;
                count++;
            }
        }
        void weight_by_volume_inverse() {};
        
        void set_fixed_vertices(std::vector<unsigned int> ids) {};

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
        void save_face_error(std::string f_out);
        void save_face_error_derivative(std::string f_out);


    protected:

    private:
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
