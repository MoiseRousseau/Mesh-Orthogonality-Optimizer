#include <iostream>
#include <string>
#include <chrono>
#include <omp.h>

#include "Connection.h"
#include "Mesh.h"
#include "OrthOpt.h"
#include "Point.h"

#include <Eigen/Core>
#include <include/LBFGS.h>

using namespace std;


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
                grad[index] = p.x;
                grad[index+1] = p.y;
                grad[index+2] = p.z;
                //TODO some discrepancy with the python code
                index += 3;
            }
            cout << "Cost function: " << opt->cost_function_value << '\n';
            //std::ostringstream f_out;
            //f_out << "derivative_it" << iteration << ".dat";
            //opt->save_face_error_derivative(f_out.str());
            return opt->cost_function_value;
        }

        Eigen::VectorXd initialize_solution_vec() {
            Eigen::VectorXd x = Eigen::VectorXd::Zero(3*opt->n_vertices_to_opt);
            unsigned int index=0;
            for (Vertice* v : opt->mesh->vertices) {
                if (v->fixed == false) {
                    x[index] = v->coor->x;
                    x[index+1] = v->coor->y;
                    x[index+2] = v->coor->z;
                    index += 3;
                }
            }
            return x;
        }

    private:

};


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
    std::cout << "Usage: " << argv[0] << " [input_mesh] [output_optimized_vertices] [opt_parameters]" << endl << endl;
    std::cout << "Specify the input mesh using:" << endl;
    std::cout << "\t-v <path to mesh coordinates in xyz format>" << endl;
    std::cout << "\t-e <path to mesh elements defined by their vertices id>" << endl;
    std::cout << "Specify the optimized vertices position output (xyz format) using:" << endl;
    std::cout << "\t-o <output file> (Defaut \"out.xyz\")" << endl;
    std::cout << "Optimization parameters:" << endl;
    std::cout << "\t-penalizing_power <value> (Penalize face error with the specified power, default = 1)" << endl;
    std::cout << "\t-maxit <value> (Maximum number of iteration, default = 20)" << endl;
    std::cout << "\t-eps <value> ()" << endl;
    std::cout << "Miscellaneous parameters:" << endl;
    std::cout << "\t-q (Quiet operation)" << endl;
    cout << "\t-n_threads (Number of thread to run in parallel, defaut = OpenMP decide)" << endl;
    std::cout << endl;
}

int get_argument(int argc, char* argv[], string arg) {
    int count;
    for (count=0; count!=argc; count++) {
        if (arg.compare(string(argv[count])) == 0) {return count;}
    }
    return 0;
}





int main(int argc, char* argv[])
{
    print_program_information();
    int index;
    std::string f_vertices;
    std::string f_elements;
    std::string f_output;
    double penalizing_power = 1.;
    int maxit = 20;

    //parse input
    if (argc < 2 or get_argument(argc, argv, "-h") != 0) {
        // Tell the user how to run the program
        print_program_help(argv);
        return 0;
    }
    //vertices
    index = get_argument(argc, argv, "-v");
    if (index == 0) {cerr << "Vertice coordinates not provided" << endl; return 1;}
    try {f_vertices = argv[index +1];}
    catch(...) {cerr << "Error while reading vertice coordinates" << endl; return 1;}
    //elements
    index = get_argument(argc, argv, "-e");
    if (index == 0) {cerr << "Elements topology not provided" << endl; return 1;}
    try {f_elements = argv[index +1];}
    catch(...) {cerr << "Error while reading element topology" << endl; return 1;}
    //output
    index = get_argument(argc, argv, "-o");
    if (index == 0) {f_output = "out.xyz";}
    else {
        try {f_output = argv[index +1];}
        catch(...) {cerr << "Error while reading output file" << endl; return 1;}
    }

    //penalizing power
    index = get_argument(argc, argv, "-penalizing_power");
    try {
        if (index != 0) {penalizing_power = atof(argv[index +1]);}
    }
    catch(...) {cerr << "Error while reading penalizing power" << endl; return 1;}
    //max iteratin
    index = get_argument(argc, argv, "-maxit");
    try {
        if (index != 0) {maxit = atof(argv[index +1]);}
    }
    catch(...) {cerr << "Error while reading maximum iteration" << endl; return 1;}
    //thread number
    index = get_argument(argc, argv, "-t");
    try {
        if (index != 0) {omp_set_num_threads(atoi(argv[index +1]));}
    }
    catch(...) {cerr << "Error while reading number of thread" << endl; return 1;}


    auto t1 = chrono::high_resolution_clock::now();

    Mesh mesh;
    //mesh.load_vertices_xyz(f_vertices);
    //mesh.load_elements_with_vertice_ids(f_elements);
    cout << "Read mesh" << endl;
    mesh.load_vertices_tetgen(f_vertices);
    mesh.load_elements_tetgen(f_elements);
    cout << "Number of vertices in mesh: " << mesh.vertices.size() << endl;
    cout << "Number of elements in mesh: " << mesh.elements.size() << endl;

    cout << endl << "Build internal connections" << endl;
    OrthOpt opt(&mesh);
    opt.set_penalizing_power(penalizing_power);
    cout << "Number of internal connections in mesh: " << opt.connections.size() << endl;
    cout << "Number of vertices to optimize: " << opt.n_vertices_to_opt << endl;
    opt.save_face_non_orthogonality_angle("./face_error_initial.txt");
    auto t2 = chrono::high_resolution_clock::now();
    cout << "Time elapsed: " << (t2-t1).count()/1e9 << " s"<< endl;

    cout << endl << "Optimize" << endl;
    //optimizer
    Wrapper_for_LBFGS wrapper(&opt);
    //create initial guess to pass to the solver
    //the vertice position actually
    Eigen::VectorXd x = wrapper.initialize_solution_vec();
    //create solver
    LBFGSpp::LBFGSParam<double> param;
    param.max_iterations = maxit;
    param.epsilon = 1e-6;
    LBFGSpp::LBFGSSolver<double> solver(param);
    //solve
    double fx;
    unsigned int niter;
    try {
        niter = solver.minimize(wrapper, x, fx);
    }
    catch (runtime_error) {
        cout << "Maximum iteration reached, stop..." << endl;
    }
    catch (logic_error) {
        cout << endl << "Logic error" << endl << endl;
        //opt.save_face_error_derivative("derivative_error.dat");
        return 1;
    }
    //what to do the x optimized
    cout << "Optimized cost function: " << fx << endl;

    opt.mesh->save_vertices_tetgen(f_output);
    opt.save_face_non_orthogonality_angle("./face_error_final.txt");

    t2 = chrono::high_resolution_clock::now();
    cout << "Time elapsed: " << (t2-t1).count()/1e9 << " s"<< endl;
    return 0;
}
