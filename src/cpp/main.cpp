#include <iostream>
#include <string>

#include "Connection.h"
#include "Mesh.h"
#include "OrthOpt.h"
#include "Point.h"

#include <Eigen/Core>
#include <include/LBFGS.h>


using namespace std;


void print_program_information() {
    cout << endl;
    cout << " ###########################################" << endl;
    cout << " #                                         #" << endl;
    cout << " #  OrthOpt: Mesh Orthogonality Optimizer  #" << endl;
    cout << " #                                         #" << endl;
    cout << " #        By: Moise Rousseau (2020)        #" << endl;
    cout << " #                                         #" << endl;
    cout << " ###########################################" << endl;
    cout << endl;
}

void print_program_help(char* argv[]) {
    std::cout << "Usage: " << argv[0] << " [input_mesh] [output_optimized_vertices] [parameters]" << std::endl;
    std::cout << "Specify the input mesh using:" << std::endl;
    std::cout << "\t-v <path to mesh coordinates in xyz format>" << std::endl;
    std::cout << "\t-c <path to mesh elements defined by their vertices id>" << std::endl;
    std::cout << "Specify the optimized vertices position (xyz format) using" << std::endl;
    std::cout << "\t-o <output file>" << std::endl;
    std::cout << std::endl;
}



class Wrapper_for_LBFGS
{
    public:

        OrthOpt* opt;
        LBFGSpp::LBFGSParam<double> param;
        unsigned int iteration = 0;

        Wrapper_for_LBFGS(OrthOpt* opt_in) {
            opt = opt_in;
        }

        double operator() (const Eigen::VectorXd& x,
                            Eigen::VectorXd& grad) {
            //print info to user
            iteration++;
            cout << "Iteration " << iteration << '\n';
            //update vertices position from x first
            opt->update_vertices_position(x);
            //then compute cost and derivative
            opt->computeCostDerivative();
            unsigned int index = 0;
            for (Point p : opt->face_error_derivative) {
                //cout << p.x << " " << p.y << ' ' << p.z << '\n';
                grad[index] = -p.x;
                grad[index+1] = -p.y;
                grad[index+2] = -p.z; //TODO investigate the minus
                //TODO some discrepancy with the python code
                index += 3;
            }
            cout << "Cost function: " << opt->cost_function_value << '\n';
            return opt->cost_function_value;
        }

    private:

};



int main(int argc, char* argv[])
{
    print_program_information();
    //parse input
    if (argc < 2) {
    // Tell the user how to run the program
        std::cerr << "No input mesh, run " << argv[0] << " -h for help" << std::endl;
        std::cerr << std::endl;
        //return 1;
    }
    if (argc == 2) {
        print_program_help(argv);
        return 0;
    }
    std::string f_vertices = "./588_pts.node"; //TODO nan with 10k
    std::string f_elements = "./588_pts.ele";
    Mesh mesh;
    //mesh.load_vertices_xyz(f_vertices);
    //mesh.load_elements_with_vertice_ids(f_elements);
    mesh.load_vertices_tetgen(f_vertices);
    mesh.load_elements_tetgen(f_elements);
    cout << "Number of vertices in mesh: " << mesh.vertices.size() << '\n';
    cout << "Number of elements in mesh: " << mesh.elements.size() << '\n';
    OrthOpt opt(&mesh);
    opt.set_penalizing_power(1);
    cout << "Number of internal connections in mesh: " << opt.connections.size() << '\n';
    cout << "Number of vertices to optimize: " << opt.n_vertices_to_opt << '\n';
    opt.save_face_non_orthogonality_angle("./face_error_initial.txt");
    //optimizer
    Wrapper_for_LBFGS wrapper(&opt);

    //create initial guess to pass to the solver
    //the vertice position actually
    Eigen::VectorXd x = Eigen::VectorXd::Zero(3*opt.n_vertices_to_opt);
    unsigned int index=0;
    for (Vertice* v : opt.mesh->vertices) {
        if (v->fixed == false) {
            x[index] = v->coor->x;
            x[index+1] = v->coor->y;
            x[index+2] = v->coor->z;
            index += 3;
        }
    }
    //create solver
    LBFGSpp::LBFGSParam<double> param;
    param.max_iterations = 20;
    param.epsilon = 1e-6;
    LBFGSpp::LBFGSSolver<double> solver(param);
    //solve
    double fx;
    unsigned int niter;
    try {
        niter = solver.minimize(wrapper, x, fx);
    }
    catch (runtime_error) { }
    //what to do the x optimized
    cout << "Optimized cost function: " << fx << '\n';

    opt.mesh->save_vertices_tetgen("./10k_pts_opt.node");
    opt.save_face_non_orthogonality_angle("./face_error_final.txt");
    return 0;
}
