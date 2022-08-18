#ifndef ORTHOPT_H
#define ORTHOPT_H

#include <array>
#include "Mesh.h"
#include "Connection.h"
#include "Error_Functions.h"
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
        std::vector<Connection*> faces;
        std::vector<double> face_weight;
        //for each point, store the position of its derivative in face_error_derivative
        // 0 if fixed
        std::vector<unsigned int> derivative_vertice_ids;
        
        // WEIGHTING FUNCTION //
        Error_Function* Ef;
        
        //TODO constructor when mesh is passed by numpy
        OrthOpt() {};

        OrthOpt(Mesh* m, Error_Function* Ef_=nullptr) {
            if (Ef_ == nullptr) {Ef_ =  new Id_Function();}
            Ef = Ef_;
            mesh = m;
            if (mesh->n_connections_internal == 0) {mesh->decompose();};
            n_vertices_to_opt = mesh->n_vertices; //initiate to n_vertices and substract the fixed point);
            face_weight.resize(mesh->n_connections_internal);
            derivative_vertice_ids.resize(mesh->n_vertices);
            unsigned int index = -1;
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
        }
        

        virtual ~OrthOpt() {};


        void set_error_function(Error_Function* Ef_) {Ef = Ef_;}
        
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
        
        void set_fixed_vertices(std::vector<unsigned int> ids) {};

        void computeCostFunction();
        void computeCostDerivative(Eigen::VectorXd& grad); //in place assignement
        
        void decompose_mesh();

        //update vertices position with eigen vector
        void update_vertices_position(const Eigen::VectorXd &x);

        //output function
        void save_face_error(std::string f_out);
        void save_face_error_derivative(std::string f_out);


    protected:

    private:
        Eigen::Vector3d derivative_E_position(Connection* con) {
            return ((con->normal - con->cell_center_vector * (con->orthogonality)) \
                    / con->cell_center_vector_norm * 0.25);
        }

        Eigen::Vector3d derivative_A_position(Connection* con, Vertice* A) {
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
                if (con->vertices.size() == 3) {C = con->vertices[0];}
                else {C = con->vertices[3];}
            }
            else if (index == 2) {
                if (con->vertices.size() == 3) {B = con->vertices[0]; C = con->vertices[1];}
                else {B = con->vertices[3]; C = con->vertices[0];}
            }
            else {B = con->vertices[0]; C = con->vertices[1];}
            //compute normal derivative
            Eigen::Vector3d temp, BC;
            temp = con->cell_center_vector - con->orthogonality*con->normal;
            BC = *(C->coor)-*(B->coor);
            temp = temp.cross(BC) / con->area;
            if (con->vertices.size() == 3) { //face is a triangle
                return 0.5*temp;
            }
            else return temp;
        }
};

#endif // ORTHOPT_H
