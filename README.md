# InductivelyCoupledPlasma 
A Fortran program that solves the two-dimensional induction equation for the electric field inside an inductively coupled plasma (ICP) torch.

This program solves the induction equation with a first order Galerkin finite element method. The program uses the parallel direct solver MUMPS to solve the system of equations resulting from the discretization. 
Further information can be found in the documentation [InductivelyCoupledPlasma.pdf](doc/InductivelyCoupledPlasma.pdf).

## Organization

* [**Doc**](https://github.com/xavierdechamps/InductivelyCoupledPlasma/tree/master/doc): this folder contains the PDF file that explains the theory and shows the results from the program on several test cases.
* [**SRC**](https://github.com/xavierdechamps/InductivelyCoupledPlasma/tree/master/SRC): this folder contains the source code

## Get the code

You can get the latest code by cloning the master branch:

```
git clone https://github.com/xavierdechamps/InductivelyCoupledPlasma.git
```

## Build the code

The code can be compiled by CMake as the necessary files ([CMakeLists.txt](CMakeLists.txt) and [CMake.config](CMake.config)) are provided.
The user can modify the content of [CMake.config](CMake.config) to his/her own configuration.

```
mkdir build
cd build
cmake ..
make
```
