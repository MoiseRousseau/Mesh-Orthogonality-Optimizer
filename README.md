# OrthOpt: a tetrahedral mesh orthogonality optimizer

Optimize your tetrahedral mesh to minimize its non-orthogonality for solving elliptic PDE using the finite volume method with the TPFA approximation.

Uses LBFGS++ library

## Getting started

### Installation

0. Install dependencies: `sudo apt install libeigen3-dev`
1. Clone this repository `git clone XXX`
2. Launch the makefile `make`

### Python Wrapper

Highly experimental (and it doesn't work)

1. Open makefile and add option `-fPIC to TARGET`
2. Compile the binder: `g++ -c -fPIC binding.cpp -o binding.o -I/usr/include/eigen3`
3. Create the shared library: 
`g++ -shared -Wl,--no-undefined,-soname,lib.so -o lib.so  ../../build/obj/Mesh.o ../../build/obj/OrthOpt.o binding.o -fopenmp`

## Examples

## Road map

* Scan mode (compute the mesh non-orthogonality and display statistics)
* Parallel decomposition of mesh
* Code vectorization and optimization
* GPU computing (for fun and learning)
* More IO format
* Python binding
* Interface with [Salome](https://www.salome-platform.org) CAD


## Contributing
