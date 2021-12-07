# OrthOpt: a tetrahedral mesh orthogonality optimizer

Optimize a tetrahedral mesh to minimize its mean non-orthogonality (possibly weighted).


## Features

* Decompose the mesh (i.e. finding its internal faces) and compute statistic on non-orthogonality and skewness.
* Mesh non orthogonality is defined by the dot product of face normal and vector defined by the two cell centers sharing the face.
* Optimize the non-orthogonality distribution in the mesh according to a gradient based method using the LBGFS algorithm (from [LBFGS++](https://lbfgspp.statr.me/) library) by moving its vertices.
* Optimize vertices shared by only tetrahedra (e.g. if the vertices belong to one hex element and 100 tetrahedra, it is NOT optimized).
* Control the non-orthogonality distribution in the mesh by applying penalization function to the non-orthogonality (power function, inverse function, log function, exponential function).
* Read and write mesh in various formats (Medit, TetGen, PFLOTRAN for instance)


## Installation

1. Install dependencies: `sudo apt install libeigen3-dev`
2. Clone this repository `git clone https://github.com/MoiseRousseau/Mesh-Orthogonality-Optimizer.git OrthOpt && cd OrthOpt`
3. Launch the makefile `make`


## Use

Please see OrthOpt launching OrthOpt without arguments (e.g. `./OrthOpt`)

For example, for a mesh in Medit format (extension `.mesh`), the command
```
./OrthOpt -optimize -m in.mesh -o out.mesh -function_type 0 -penalizing_power 5 -maxit 100 
```
optimize the vertices position of the mesh stored in the file `in.mesh` (Medit format) considering a power weigthing function with a penalization power of 5 with a maximum of 100 iteration.
Optimized mesh is saved in Medit format in the file `out.mesh`.


## Salome plugin

A plugin for the CAD software [Salome](https://salome-platform.org/) plugin is also available.
The plugin export the mesh in a PFLOTRAN-like format, automatically call the optimizer, and reimport the optimize mesh in Salome.

Installation is done in three step:
1. Copy the content or the folder `salome_plugin` in this repository to a newly created folder nammed `OrthOpt` in the directory `$HOME/.config/salome/Plugins/`
2. Add the following two lines at the end of `smesh_plugins.py` file located in the previous directory:
```
import OrthOpt
salome_pluginsmanager.AddFunction('OrthOpt', ' ',
                                  OrthOpt.OrthOpt_opt)
```
3. Copy the `OrthOpt` executable in the `$HOME/.config/salome/Plugins/OrthOpt/` folder.

You are ready to go. The plugin is located in the MESH module under the `Mesh/SMESH Plugins/OrthOpt`.


## Examples

Example of the optimization of a mesh composed of nearly 25000 nodes and 125000 tetrahedra for simulation the flow around a cylinder as this OpenFOAM [tutorial](https://www.openfoam.com/documentation/tutorial-guide/2-incompressible-flow/2.2-flow-around-a-cylinder#x7-390002.2).
Mesh was generated using [Salome](https://salome-platform.org/) and the NETGEN algorithm with default optimization parameters and refined near the cylinder.
The table below summarizes the mesh statistics and the total number of iteration for the linear solver with the orthogonal correction for the Laplacian to converge (with tolerance 1e-6 and 1e-12).

| Mesh | Mean non-orthogonality (°) | Max non-orthogonality (°) | Solver iteration (1e-6) | Solver iteration (1e-12) |
|---|---|---|---|---|
| OrthOpt n=0.5 | **11.2** | 70.8 | 34 | > 100 |
| OrthOpt n=1 | 11.7 | 56.9 | 31 | 78 |
| Netgen default | 12.4 | 50.7 | 29 | 68 |
| OrthOpt n=3 | 13.6 | 44.1 | 29 | 59 |
| OrthOpt n=5 | 14.5 | **42.6** | **28** | **52** |

The figure below also shows the distribution of non-orthogonality angle within the internal faces of the meshes.

![Mesh non-orthogonality distribution example](https://github.com/MoiseRousseau/Mesh-Orthogonality-Optimizer/blob/master/examples.png)

In details, increasing the penalization power increase the mean non-orthogonality, but decrease its maximum.
Decrease of the maximum resulted in fewer iteration for the solver to converge, especially for a tolerance of 1e-12.


## TODO

* Parallel decomposition of mesh (maybe use some dedicaced library)
* Hybrid parallelization (OpenMP and MPI) ? GPU computing ?
* Some code optimization (e.g. store vertices coordinate as one table, and not as a structure)
* More IO format
* Add more examples
* Tests
* Write a paper ?


### Python Wrapper

Rather some note for me than to be used.

1. Open makefile and add option `-fPIC to TARGET`
2. Compile the binder: `g++ -c -fPIC binding.cpp -o binding.o -I/usr/include/eigen3`
3. Create the shared library: 
`g++ -shared -Wl,--no-undefined,-soname,lib.so -o lib.so  ../../build/obj/Mesh.o ../../build/obj/OrthOpt.o binding.o -fopenmp`
