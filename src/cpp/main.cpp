#include <iostream>
#include <string>
#include <chrono>
#include <omp.h>

#include "Connection.h"
#include "Mesh.h"
#include "OrthOpt.h"
#include "Point.h"
#include "IO.h"


#include <Eigen/Core>
#include "include/LBFGS.h"

using namespace std;

//Python binding

#if 0
extern "C" {
    OrthOpt* OrthOpt() {
        OrthOpt* res = new OrthOpt()
        return res;
    }
    OrthOpt* OrthOpt(Mesh* m) {return new OrthOpt(m); }

    Mesh* Mesh() {return new Mesh(); }

    void load_vertices_array(Mesh mesh, double vertices[], size_t s) {
        Mesh::load_elements_array(vertices, s);
    }

    void load_elements_array(Mesh mesh, double vertices[], size_t s) {
        Mesh::load_elements_array(vertices, s);
    }

    unsigned int optimize(OrthOpt* obj, double penalization_power,
                          unsigned int max_it) {

        return code
    }
}
#endif


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


int optimize_mesh(Mesh* mesh, 
                  int penalizing_power, int weighting_method,
                  int maxit, 
                  std::string f_output, std::string f_fixed) { 
    
    
    //set up optimization problem class
    OrthOpt opt(mesh);
    opt.set_penalizing_power(penalizing_power);
    switch (weighting_method) {
        case 0:
            opt.weight_unitary();
            break;
        case 1: 
            opt.weight_by_area();
            break;
        case 2: 
            opt.weight_by_area_inverse();
            break;
        case 3:
            opt.weight_by_volume_inverse();
            break;
    }
    //fix point
    if (f_fixed.length() >= 1) {
        std::vector<unsigned int> ids_fixed;
        //TODO
        opt.set_fixed_vertices(ids_fixed);
    }
    
    //setup optimizer
    Wrapper_for_LBFGS wrapper(&opt);
    //create initial guess to pass to the solver
    //the vertice position actually
    Eigen::VectorXd x = wrapper.initialize_solution_vec();
    //create solver
    LBFGSpp::LBFGSParam<double> param;
    param.max_iterations = maxit;
    param.epsilon = 1e-6;
    LBFGSpp::LBFGSSolver<double> solver(param);
    
    //pre-solve
    cout << "Using penalization power " << opt.penalizing_power << endl;
    cout << "Number of vertices to optimize: " << opt.n_vertices_to_opt << endl;
    mesh->save_face_non_orthogonality_angle("./face_error_initial.txt");
    cout << endl << "Optimize" << endl;
    
    //solve
    double fx;
    unsigned int niter;
    try {
        niter = solver.minimize(wrapper, x, fx);
    }
    catch (runtime_error) {
        cout << "Maximum iteration reached, stop..." << endl;
        return 1;
    }
    catch (logic_error) {
        cout << endl << "Logic error" << endl;
        cout << "You are not supposed to be here and probably find a bug..." << endl;
        cout << "You can report it, or ignore it." << endl;
        cout << "Or you ca change the penalization power to make thing work" << endl;
        cout << endl;
        //opt.save_face_error_derivative("derivative_error.dat");
        return 1;
    }
    cout << "Optimized cost function: " << fx << endl;

    mesh->save_face_non_orthogonality_angle("./face_error_final.txt");
    return 0;
}

int scan_mesh(Mesh &m, std::string f_output) {
    cout << endl << "Build internal connections" << endl;
    m.decompose();
    cout << "Number of internal connections in mesh: " << m.connections_internal.size() << endl;
    cout << "Save face informations" << endl;
    m.save_face_informations(f_output);
    m.display_stats();
    return 0;
}


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
    cout << "Usage: " << argv[0] << "[mode] [input_mesh] [output_optimized_vertices] [opt_parameters]" << endl << endl;
    cout << "Mode:" << endl;
    cout << "-scan (read the mesh and return statistics on mesh area, center, non-orthogonality and skewness)"  << endl;
    cout << "-optimize (optimize the mesh vertice position - default)" << endl << endl;
    cout << "Input mesh:" << endl;
    cout << "\t-v <path to mesh coordinates in xyz format>" << endl;
    cout << "\t-e <path to mesh elements defined by their vertices id>" << endl;
    cout << "\t-fixed <path to fixed vertices id> (optional)" << endl << endl;
    cout << "Optimized vertices position output (xyz format):" << endl;
    cout << "\t-o <output file> (Defaut \"out.xyz\")" << endl << endl;
    cout << "Optimization parameters:" << endl;
    cout << "\t-penalizing_power <float> (Penalize face error with the specified power, default = 1)" << endl;
    cout << "\t-maxit <int> (Maximum number of iteration, default = 20)" << endl;
    cout << "\t-face_weighting <int> (Method to weight face error, 0 = no weighting, ";
    cout << "1 = weight by face area, 2 = weight by face area inverse, default = 0)"<< endl << endl;
    //cout << "\t-eps <value> ()" << endl;
    cout << "Miscellaneous parameters:" << endl;
    cout << "\t-q (Quiet operation)" << endl;
    cout << "\t-n_threads <int> (Number of thread to run in parallel, defaut = OpenMP decide)" << endl;
    cout << endl;
}

int get_argument(int argc, char* argv[], string arg) {
    int count;
    for (count=0; count!=argc; count++) {
        if (arg.compare(string(argv[count])) == 0) {return count;}
    }
    return 0;
}





int main(int argc, char* argv[]) {

    print_program_information();
    //defaut parameter
    std::string f_vertices = "";
    std::string f_elements = "";
    std::string f_fixed = "";
    std::string f_output = "out.xyz";
    double penalizing_power = 1.;
    int mode = 0; //0=optimize, 1=scan
    int maxit = 20;
    int weighting_method = 0; //0=no weighting, 1=face area, 2=face area inverse
    bool quiet = false;
    bool surface = false;

    //parse input
    if (argc < 2 or get_argument(argc, argv, "-h") != 0) {
        // Tell the user how to run the program
        print_program_help(argv);
        return 0;
    }

    int iarg = 1;
    auto arg = argv[0];
    while (iarg < argc) {
        arg = argv[iarg];
        if (!strcmp(arg, "-v")) {
            iarg++; f_vertices = argv[iarg];
        }
        else if (!strcmp(arg, "-e")) {
            iarg++; f_elements = argv[iarg];
        }
        else if (!strcmp(arg, "-fixed")) {
            iarg++; f_fixed = argv[iarg];
        }
        else if (!strcmp(arg, "-o")) {
            iarg++; f_output = argv[iarg];
        }
        else if (!strcmp(arg, "-maxit")) {
            iarg++; maxit = atof(argv[iarg]);
        }
        else if (!strcmp(arg, "-penalizing_power")) {
            iarg++; penalizing_power = atof(argv[iarg]);
        }
        else if (!strcmp(arg, "-face_weighting")) {
            iarg++; weighting_method = atoi(argv[iarg]);
        }
        else if (!strcmp(arg, "-n_threads")) {
            iarg++; omp_set_num_threads(atoi(argv[iarg]));
        }
        else if (!strcmp(arg, "-q")) {
            quiet = true;
        }
        else if (!strcmp(arg, "-optimize")) {
            mode=0;
        }
        else if (!strcmp(arg, "-scan")) {
            mode=1;
        }
        else {
            cerr << "Argument not recognized " << arg << endl;
            return 1;
        }
        iarg++;
    }
    //check mandatory input
    if (f_vertices.size() == 0) {cerr << "Vertice coordinates not provided" << endl; return 1;}
    if (f_elements.size() == 0) {cerr << "Elements topology not provided" << endl; return 1;}

    //load mesh
    auto t1 = chrono::high_resolution_clock::now();
    cout << "Read mesh" << endl;
    Mesh mesh;
    IO io_mesh(&mesh);
    //IO::load_mesh_auto(mesh, )
    io_mesh.load_vertices_tetgen(f_vertices);
    io_mesh.load_elements_tetgen(f_elements);
    cout << "Number of vertices in mesh: " << mesh.vertices.size() << endl;
    cout << "Number of elements in mesh: " << mesh.elements.size() << endl;
    
    int ret = 0;
    switch (mode) {
        case 0:
            ret = optimize_mesh(&mesh, penalizing_power, weighting_method, maxit, f_output, f_fixed);
            io_mesh.save_vertices_tetgen(f_output);
            break;
        case 1: 
            ret = scan_mesh(mesh, f_output);
            break;
        default:
            break;
    }
    
    auto t2 = chrono::high_resolution_clock::now();
    cout << "Time elapsed: " << (t2-t1).count()/1e9 << " s"<< endl;
    return ret;
}
