#include <iostream>
#include <string>
#include <chrono>
#include <omp.h>

#include "Connection.h"
#include "Mesh.h"
#include "OrthOpt.h"
#include "Point.h"
#include "IO.h"
#include "Error_Functions.h"


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

        Wrapper_for_LBFGS(OrthOpt* opt_in) {
            opt = opt_in;
        }

        double operator() (const Eigen::VectorXd& x,
                            Eigen::VectorXd& grad) {
            //print info to user
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
                index += 3;
            }
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
                  Error_Function* Ef, int weighting_method,
                  int maxit, double eps,
                  std::string f_output, std::string f_fixed) { 
    
    
    //set up optimization problem class
    OrthOpt opt(mesh, Ef);
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
    param.epsilon = eps;
    param.max_linesearch = 100;
    LBFGSpp::LBFGSSolver<double> solver(param);
    
    //pre-solve
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
        mesh->save_face_non_orthogonality_angle("./face_error_final.txt");
        throw;
        return 1;
    }
    catch (logic_error) {
        cout << endl << "Logic error" << endl;
        cout << "You are not supposed to be here and probably find a bug..." << endl;
        cout << "You can report it, or ignore it." << endl;
        cout << "Or you ca change the penalization power to make thing work" << endl;
        cout << endl;
        //opt.save_face_error_derivative("derivative_error.dat");
        throw;
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
    
    m.save_face_detailed_informations(f_output);
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
    cout << "Usage: " << argv[0] << " [mode] [input_mesh] \
[output_optimized_vertices] [opt_parameters]" << endl;
    cout << endl;
    cout << "Mode:" << endl;
    cout << "\t-scan (read the mesh and return statistics on mesh area, center, \
non-orthogonality and skewness)"  << endl;
    cout << "\t-optimize (optimize the mesh vertice position - default)" << endl;
    cout << endl;
    cout << "Input mesh:" << endl;
    cout << "\t-v <path to mesh coordinates in xyz format>" << endl;
    cout << "\t-e <path to mesh elements defined by their vertices id>" << endl;
    cout << "\t-fixed <path to fixed vertices id> (optional)" << endl;
    cout << endl;
    cout << "Optimized vertices position output (xyz format):" << endl;
    cout << "\t-o <output file> (Defaut \"out.xyz\")" << endl;
    cout << endl;
    cout << "Error function parameters:" << endl;
    cout << "\t-function_type <int> (Method to compute penalize error, 0 = power function, 1 = inverse function, 2 = log function, default = 0)" << endl;
    cout << "\t-penalizing_power <float> (Power or inverse function, penalize\
face error with the specified power, default = 1.)" << endl;
    cout << "\t-face_weighting <int> (Method to weight face error, \
0 = no weighting, 1 = weight by face area,\
2 = weight by face area inverse, default = 0)"<< endl;
    cout << "\t-face_weighting_file <path>" << endl;
    cout << endl;
    cout << "Optimization parameters:" << endl;
    cout << "\t-maxit <int> (Maximum number of iteration, default = 20)" << endl;
    //cout << "\t-eps <value> ()" << endl;
    cout << "Miscellaneous parameters:" << endl;
    cout << "\t-q (Quiet operation)" << endl;
    cout << "\t-n_threads <int> (Number of thread to run in parallel,\
defaut = OpenMP decide)" << endl;
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
    int function_type = 0;
    int maxit = 20;
    double eps = 1e-6;
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
        // MODE //
        if (!strcmp(arg, "-optimize")) {
            mode=0;
        }
        else if (!strcmp(arg, "-scan")) {
            mode=1;
        }
        // INPUT //
        else if (!strcmp(arg, "-v")) {
            iarg++; f_vertices = argv[iarg];
        }
        else if (!strcmp(arg, "-e")) {
            iarg++; f_elements = argv[iarg];
        }
        else if (!strcmp(arg, "-fixed")) {
            iarg++; f_fixed = argv[iarg];
        }
        // OUTPUT //
        else if (!strcmp(arg, "-o")) {
            iarg++; f_output = argv[iarg];
        }
        // ERROR FUNCTION //
        else if (!strcmp(arg, "-function_type")) {
            iarg++; function_type = atoi(argv[iarg]);
        }
        else if (!strcmp(arg, "-penalizing_power")) {
            iarg++; penalizing_power = atof(argv[iarg]);
        }
        else if (!strcmp(arg, "-face_weighting")) {
            iarg++; weighting_method = atoi(argv[iarg]);
        }
        // OPTIMIZATION PARAMETER //
        else if (!strcmp(arg, "-maxit")) {
            iarg++; maxit = atoi(argv[iarg]);
        }
        else if (!strcmp(arg, "-eps")) {
            iarg++; eps = atof(argv[iarg]);
        }
        else if (!strcmp(arg, "-n_threads")) {
            iarg++; omp_set_num_threads(atoi(argv[iarg]));
        }
        else if (!strcmp(arg, "-q")) {
            quiet = true;
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
    
    if ((function_type < 0) or (function_type) > 2) {
        cerr << "Error! function type not recognized: " << function_type << endl;
        return 1;
    }
    
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
    
    Error_Function *Ef;
    if (mode == 0) {
        if (mode == 0) {
            switch (function_type) {
                case 0:
                    Ef = new Power_Function(penalizing_power);
                    cout << "Penalize orthogonality error with power\
 function at power " << penalizing_power << endl;
                    break;
                case 1:
                    Ef = new Inverse_Function (penalizing_power);
                    cout << "Penalize orthogonality error with inverse\
 function at power " << penalizing_power << endl;
                    break;
                case 2:
                    Ef = new Log_Function();
                    cout << "Penalize orthogonality error with logarithm\
 function" << endl;
                    break;
                case 3: 
                    Ef = new Tan_Function();
                    cout << "Penalize orthogonality error with tangent\
 function" << endl;
                    break;
                default:
                    break;
            }
        }
    }
    
    int ret = 0;
    switch (mode) {
        case 0:
            ret = optimize_mesh(&mesh, Ef, weighting_method, 
                                maxit, eps, f_output, f_fixed);
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
