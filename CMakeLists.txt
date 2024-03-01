cmake_minimum_required(VERSION 3.15)

set ( CMAKE_Fortran_COMPILER      "ifort" )

# Declaration du projet ####################################################
project(MyProject
LANGUAGES Fortran
VERSION 1.0
DESCRIPTION "Finite element solver for the 2D electric field induction equation")

include( "config.cmake" )

# Include the flags in the compiler options #####################################
include_directories(${MUMPS_INCLUDE_DIRS})

if ( Windows )
   message ( "   Windows environment..." )
   set ( ENVI_FLAG  "-DWINDOWS" )
else ()
   message ( "   Linux/MAX environment..." )
   set ( ENVI_FLAG  "" )
endif()

# The executables are put in the source_dir/EXE directory ##########################
message ( "The executables are to be put in ${CMAKE_CURRENT_SOURCE_DIR}/EXE" )
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/EXE)

# Generation de la liste des fichiers sources ####################################
file(GLOB_RECURSE SRC_file SRC/ module_icp.f90 )

# Extract the base of the paths of the previous files #############################
get_filename_component (SRC_dir ${SRC_file} PATH )

# Construction of the executables ###########################################
add_library ( libmodules
				      ${SRC_dir}/module_icp.f90 )
                      
# Construction of the executables ###########################################
### The solver itself
add_executable(icp 
				      ${SRC_dir}/csr.f90
				      ${SRC_dir}/find_elem.f90
				      ${SRC_dir}/getEcoil.f90
				      ${SRC_dir}/getStiffness.f90
				      ${SRC_dir}/gmsh_operations.f90
				      ${SRC_dir}/handling_parameters.f90
				      ${SRC_dir}/init.f90
				      ${SRC_dir}/main.f90
				      ${SRC_dir}/mem_allocate.f90
				      ${SRC_dir}/setBC.f90
				      ${SRC_dir}/solving.f90 )
add_dependencies( icp libmodules)
set_target_properties ( icp PROPERTIES LINKER_LANGUAGE Fortran)
target_link_libraries( icp PUBLIC libmodules
   			         ${BLAS_LIB}
		 	 	     ${METIS_LIBRARIES}
				     ${MUMPS_LIBRARIES}
				     ${SCALAPACK_LIBRARIES} )
                     
### Build Initial Condition = Build sigma
add_executable(build_initial_solution 
				      ${SRC_dir}/Build_initial_condition.f90
				      ${SRC_dir}/gmsh_operations.f90
				      ${SRC_dir}/mem_allocate.f90
				      ${SRC_dir}/handling_parameters.f90 
				      ${SRC_dir}/find_elem.f90 )
add_dependencies( build_initial_solution libmodules)
set_target_properties ( build_initial_solution PROPERTIES LINKER_LANGUAGE Fortran)
target_link_libraries( build_initial_solution libmodules )
