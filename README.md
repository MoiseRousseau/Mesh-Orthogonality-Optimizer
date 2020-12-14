# OrthOpt: a tetrahedral mesh orthogonality optimizer

OrthOpt is a 

Uses LBFGS++ library

## Getting started

### Installation

0. Install dependencies: `sudo apt install libeigen3-dev`
1. Clone this repository `git clone XXX`
2. Launch the makefile `make`

### Python Wrapper

1. Open makefile and add option `-fPIC to TARGET`
2. Compile the binder: `g++ -c -fPIC binding.cpp -o binding.o -I/usr/include/eigen3`
3. Create the shared library: 
`g++ -shared -Wl,--no-undefined,-soname,lib.so -o lib.so  ../../build/obj/Mesh.o ../../build/obj/OrthOpt.o binding.o -fopenmp`

## Examples
