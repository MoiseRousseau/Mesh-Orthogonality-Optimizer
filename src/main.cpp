#include <iostream>
#include <string>
#include <chrono>
#include <omp.h>

#include "Connection.h"
#include "Mesh.h"
#include "OrthOpt.h"
#include "IO.h"
#include "Error_Functions.h"


#include <Eigen/Core>
#include "LBFGS.h"
#include <fstream>

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
        bool fd_grad = false;

        Wrapper_for_LBFGS(OrthOpt* opt_in) {
            opt = opt_in;
        }

        double operator() (const Eigen::VectorXd& x,
                            Eigen::VectorXd& grad) {
            //print info to user
            //update vertices position from x first
            opt->update_vertices_position(x);
            //then compute cost and derivative
            if (fd_grad) opt->computeCostDerivative_FD(grad);
            else opt->computeCostDerivative(grad);
            
            //std::ostringstream f_out;
            //f_out << "derivative_it" << iteration << ".dat";
#ifdef DEBUG_MODE
            std::cout << "Analytic gradient:" << std::endl;
            std::cout << grad << std::endl;
            std::cout << "FD gradient:" << std::endl;
            opt->computeCostDerivative_FD(grad);
            std::cout << grad << std::endl;
            std::cout << "position" << std::endl;
            std::cout << x << std::endl;
#endif
            return opt->cost_function_value;
        }

        Eigen::VectorXd initialize_solution_vec() {
            Eigen::VectorXd x = Eigen::VectorXd::Zero(opt->mesh->dim * opt->n_vertices_to_opt);
            unsigned int index=0;
            for (Vertice* v : opt->mesh->vertices) {
                if (v->fixed == false) {
                    for (size_t i=0; i<opt->mesh->dim; i++) {
                        x[index+i] = (*(v->coor))[i];
                    }
                    index += opt->mesh->dim;
                }
            }
            return x;
        }

    private:

};


int optimize_mesh(
    Mesh* mesh, 
    Error_Function* Ef, 
    int weighting_method,
    int maxit, 
    double eps, 
    double eps_rel,
    double x_rel,
    std::string f_output, 
    std::vector<unsigned int> ids_fixed,
    std::vector<double> aniso_coeff,
    bool fd_grad
) { 
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
    }
    
    //fix point
    opt.set_fixed_vertices(ids_fixed);
    
    //setup optimizer
    Wrapper_for_LBFGS wrapper(&opt);
    wrapper.fd_grad = fd_grad;
    //create initial guess to pass to the solver
    //the vertice position actually
    Eigen::VectorXd x = wrapper.initialize_solution_vec();
    //create solver
    LBFGSpp::LBFGSParam<double> param;
    param.max_iterations = maxit;
    param.epsilon = eps;
    param.epsilon_rel = eps_rel;
    param.max_linesearch = 100;
    LBFGSpp::LBFGSSolver<double> solver(param);
    
    //pre-solve
    cout << "Number of vertices to optimize: " << opt.n_vertices_to_opt << endl;
    //mesh->save_face_non_orthogonality_angle("./face_error_initial.txt");
    cout << endl << "Optimize" << endl;
    
    //solve
    double fx;
    unsigned int niter;
    try {
        niter = solver.minimize(wrapper, x, fx);
    }
    catch (runtime_error) {
        //mesh->save_face_non_orthogonality_angle("./face_error_final.txt");
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

    //mesh->save_face_non_orthogonality_angle("./face_error_final.txt");
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
    cout << " #        By: Moise Rousseau (2022)        #" << endl;
    cout << " #                                         #" << endl;
    cout << " ###########################################" << endl;
    cout << endl << "Publication DOI: TODO" << endl;
    cout << endl;
}

void print_program_help(char* argv[]) {
    cout << "Usage: " << argv[0] << " [mode] [input_mesh] \
[output_optimized_vertices] [opt_parameters]" << endl;
    cout << endl;
    cout << "Mode:" << endl;
    cout << "\t--scan (read the mesh and return statistics on mesh area, center, \
non-orthogonality and skewness)"  << endl;
    cout << "\t--optimize (optimize the mesh vertice position - default)" << endl;
    cout << endl;
    cout << "Input mesh:" << endl;
    cout << "\t--v-tetgen <path to tetgen mesh coordinates>" << endl;
    cout << "\t--e-tetgen <path to tetgen mesh elements>" << endl;
    cout << "\t--mesh / -m <path to mesh file (Medit / Salome DAT / PFLOTRAN)>" << endl;
    cout << "\t--fixed <path to fixed vertices id> (optional)" << endl;
    cout << "\t--dimension <2/3> (Mesh dimension, default 3 if not specified in mesh format, e.g. Medit)" << endl;
    cout << "\t--anisotropic_diff_coef <path to file> (Optimize a 2D considering a anisotropic and heterogeneous diffusion coefficient to respect the K-orthogonality condition, optional)" << endl;
    cout << endl;
    cout << "Optimized vertices position output (xyz format):" << endl;
    cout << "\t--output / -o <output file> (Default \"out.xyz\", output format deduced from extension)" << endl;
    cout << endl;
    cout << "Error function parameters:" << endl;
    cout << "\t--function_type <int> (Method to compute penalize error, 0 = identity, 1 = power function, 2 = inverse function, 3 = log function, 4 = exp function, default = 0)" << endl;
    cout << "\t--penalizing_power <float> (Power or inverse function, penalize\
face error with the specified power, default = 1.)" << endl;
    cout << "\t--face_weighting <int> (Method to weight face error, \
0 = no weighting, 1 = weight by face area,\
2 = weight by face area inverse, default = 0)"<< endl;
    //cout << "\t-face_weighting_file <path>" << endl;
    cout << endl;
    cout << "Optimization parameters:" << endl;
    cout << "\t--maxit <int> (Maximum number of iteration, default = 100)" << endl;
    cout << "\t--eps <value> (Absolute value of gradient norm for termination, default 1e-6)" << endl;
    cout << "\t--eps_rel <value> (Relative value of gradient norm for termination, default 1e-6)" << endl;
    //cout << "\t--x_rel <value> (Relative value of solution norm for termination, default 1e-6)" << endl;
    cout << "\t--fd_gradient (Use finite difference to compute the gradient. This is very inefficient and should be used for debug purpose)" << endl;
    cout << "Miscellaneous parameters:" << endl;
    cout << "\t--quiet / -q (Quiet operation)" << endl;
    cout << "\t--n_threads <int> (Number of thread to run in parallel,\
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
    std::string f_vertices = ""; //tetgen vertices
    std::string f_elements = ""; //tetgen element
    std::string f_mesh = ""; //general mesh (medit / ...)
    std::string f_fixed = ""; //path to fixed point in the mesh
    std::string f_output = "out.xyz"; //output name
    std::string f_aniso = ""; //anisotropic diffusion coefficient file
    double penalizing_power = 1.;
    int mode = 0; //0=optimize, 1=scan
    int function_type = 0; //identity
    int maxit = 100;
    double eps = 1e-6, eps_rel = 1e-6, x_rel = 1e-6;
    int weighting_method = 0; //0=no weighting, 1=face area, 2=face area inverse
    bool quiet = false;
    bool fd_grad = false;
    int dim = 3;

    //parse input
    if (argc < 2 or get_argument(argc, argv, "--help") != 0 or get_argument(argc, argv, "-h")) {
        // Tell the user how to run the program
        print_program_help(argv);
        return 0;
    }

    int iarg = 1;
    auto arg = argv[0];
    while (iarg < argc) {
        arg = argv[iarg];
        // MODE //
        if (!strcmp(arg, "--optimize")) {
            mode=0;
        }
        else if (!strcmp(arg, "--scan")) {
            mode=1;
        }
        // INPUT //
        else if (!strcmp(arg, "--v-tetgen")) {
            iarg++; f_vertices = argv[iarg];
        }
        else if (!strcmp(arg, "-m") or (!strcmp(arg, "--mesh"))) {
            iarg++; f_mesh = argv[iarg];
        }
        else if (!strcmp(arg, "--e-tetgen")) {
            iarg++; f_elements = argv[iarg];
        }
        else if (!strcmp(arg, "--fixed")) {
            iarg++; f_fixed = argv[iarg];
        }
        else if (!strcmp(arg, "--dimension")) {
            iarg++; dim = atoi(argv[iarg]);
        }
        else if (!strcmp(arg, "--anisotropic_diff_coef")) {
            iarg++; f_aniso = argv[iarg];
        }
        // OUTPUT //
        else if (!strcmp(arg, "-o") or !strcmp(arg, "--output")) {
            iarg++; f_output = argv[iarg];
        }
        // ERROR FUNCTION //
        else if (!strcmp(arg, "--function_type")) {
            iarg++; function_type = atoi(argv[iarg]);
        }
        else if (!strcmp(arg, "--penalizing_power")) {
            iarg++; penalizing_power = atof(argv[iarg]);
        }
        else if (!strcmp(arg, "--face_weighting")) {
            iarg++; weighting_method = atoi(argv[iarg]);
        }
        // OPTIMIZATION PARAMETER //
        else if (!strcmp(arg, "--maxit")) {
            iarg++; maxit = atoi(argv[iarg]);
        }
        else if (!strcmp(arg, "--eps")) {
            iarg++; eps = atof(argv[iarg]);
        }
        else if (!strcmp(arg, "--eps_rel")) {
            iarg++; eps_rel = atof(argv[iarg]);
        }
        else if (!strcmp(arg, "--x_rel")) {
            iarg++; x_rel = atof(argv[iarg]);
        }
        else if (!strcmp(arg, "--fd_gradient")) {
            fd_grad = true;
        }
        else if (!strcmp(arg, "--n_threads")) {
            iarg++; eps_rel = atof(argv[iarg]);
#if defined _OPENMP
#endif
        }
        else if (!strcmp(arg, "-q") or !strcmp(arg, "--quiet")) {
            quiet = true;
        }
        else {
            cerr << "Argument not recognized " << arg << endl;
            return 1;
        }
        iarg++;
    }
    
    //check mandatory input
    if ( (f_mesh.size() == 0) and
         (f_vertices.size() == 0) and
         (f_elements.size() == 0) ) {
         cerr << "Mesh not provided" << endl; 
         return 1;
    }
    
    if ((function_type < 0) or (function_type) > 4) {
        cerr << "Error! function type not recognized: " << function_type << endl;
        return 1;
    }
    
    //load mesh
    auto t1 = chrono::high_resolution_clock::now();
    Mesh mesh;
    IO io_mesh(&mesh);
    mesh.dim = dim;
    io_mesh.load_mesh_auto(f_mesh, f_vertices, f_elements);
    cout << "Mesh dimension: " << mesh.dim << endl;
    cout << "Number of vertices in mesh: " << mesh.vertices.size() << endl;
    cout << "Number of elements in mesh: " << mesh.elements.size() << endl;
    
    Error_Function *Ef;
    if (mode == 0) {
        if (mode == 0) {
            switch (function_type) {
                case 0:
                    Ef = new Id_Function();
                    cout << "Penalize orthogonality error with identity function" << endl;
                    break;
                case 1:
                    Ef = new Power_Function(penalizing_power);
                    cout << "Penalize orthogonality error with power\
 function at power " << penalizing_power << endl;
                    break;
                case 2:
                    Ef = new Inverse_Function (penalizing_power);
                    cout << "Penalize orthogonality error with inverse\
 function at power " << penalizing_power << endl;
                    break;
                case 3:
                    Ef = new Log_Function();
                    cout << "Penalize orthogonality error with logarithm\
 function" << endl;
                    break;
                case 4:
                    Ef = new Exp_Function(penalizing_power);
                    cout << "Penalize orthogonality error with exponential\
 function at power " << penalizing_power << endl;
                    break;
                case 5: 
                    Ef = new Tan_Function();
                    cout << "Penalize orthogonality error with tangent\
 function" << endl;
                    break;
                default:
                    break;
            }
        }
    }
    
    // retrieve fixed point
    std::vector<unsigned int> ids_fixed;
    unsigned int id_fixed;
    if (f_fixed.size() != 0) {
        std::ifstream src_fixed(f_fixed);
        while (src_fixed >> id_fixed) {
            ids_fixed.push_back(id_fixed);
        }
        src_fixed.close();
    }
    // get anisotropic coefficient
    // Kx Ky Kxy per line, 1 line per element
    std::vector<double> aniso_coeff;
    double Kx, Ky, Kxy;
    size_t count = 0;
    if (f_aniso.size() != 0) {
        std::ifstream src_fixed(f_fixed);
        while (src_fixed >> Kx >> Ky >> Kxy) {
            aniso_coeff.push_back(Kx);
            aniso_coeff.push_back(Ky);
            aniso_coeff.push_back(Kxy);
            count++;
        }
        if (count != mesh.elements.size()) {
            cout << "Number of diffusion coefficient does not match number of element in the mesh. You should define three diffusion coefficient components (Kx, Ky, Kxy) for each element." << endl; 
            exit(1);
        }
        src_fixed.close();
    }
    
    int ret = 0;
    switch (mode) {
        case 0:
            ret = optimize_mesh(&mesh, Ef, weighting_method, 
                                maxit, eps, eps_rel, x_rel,
                                f_output, ids_fixed,
                                aniso_coeff, fd_grad);
            io_mesh.save_mesh_auto(f_output);
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
