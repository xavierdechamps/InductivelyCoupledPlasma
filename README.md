# InductivelyCoupledPlasma 
A Fortran program that solves the two-dimensional induction equation for the electric field inside an inductively coupled plasma (ICP) torch.

This program solves the induction equation with a first order Galerkin finite element method. The program uses the parallel direct solver MUMPS to solve the system of equations resulting from the discretization. 
Further information can be found in the documentation [InductivelyCoupledPlasma.pdf](doc/InductivelyCoupledPlasma.pdf). A [walkthrough](doc/UserGuide.pdf) PDF is also included in order to help
the user using the programs.

## Organization

* [**Doc**](https://github.com/xavierdechamps/InductivelyCoupledPlasma/tree/master/doc): this folder contains the PDF file that explains the theory and shows the results from the program on several test cases.
This folder also contains a walkthrough PDF that explains the differents steps to follow in order to use the programs.
* [**SRC**](https://github.com/xavierdechamps/InductivelyCoupledPlasma/tree/master/SRC): this folder contains the source code

## Get the code

You can get the latest code by cloning the master branch:

```
git clone https://github.com/xavierdechamps/InductivelyCoupledPlasma.git
```

## Build the code

The code can be compiled by CMake as the necessary files ([CMakeLists.txt](CMakeLists.txt) and [config.cmake](config.cmake)) are provided.
The user can modify the content of [config.cmake](config.cmake) to his/her own configuration.
The user must have the library MUMPS compiled and linked at hand in order to link with the present code.

```
mkdir build
cd build
cmake ..
make
```

## Run the code

The parameters of the computation are specified in an external file, see the file [parameters](parameters) for an example.
A mesh must be built from the .geo geometry file. The mesh generator must be Gmsh. The name of the mesh is then specified in the parameter file.
Run the computation through the command
```
EXE\icp.exe parameters
```