message ( "Entering in the configuration file" )

unset ( Windows       CACHE )
unset ( MANGLING   CACHE )

# General options / Flags #########################################################
set ( Windows           ON CACHE BOOL "Is it a Windows-based environment?" )

# Fortran compiler
# set ( CMAKE_Fortran_COMPILER      "mpif90" )
# set ( CMAKE_Fortran_FLAGS         "-w -fno-range-check -m64 -m64" )
# set ( CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -fbounds-check" )
# set ( CMAKE_Fortran_FLAGS_RELEASE "-O3" )
set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} /heap-arrays:0 /QxSSE4.2")
set ( CMAKE_Fortran_FLAGS_RELEASE "/O3" CACHE STRING "Fortran options release" )
set ( CMAKE_EXE_LINKER_FLAGS "/NODEFAULTLIB:library" )

# C compiler / choice between    UPPERCASE_    lowercase_    lowercase    depending on the compilator
# set ( MANGLING              "lowercase_" )
set ( MANGLING "UPPERCASE" CACHE STRING "The standard rename procedure for C function during the linking phase. You may choose between UPPERCASE_    lowercase_    lowercase")

# set ( CMAKE_C_FLAGS         "-m64 -lguide -openmp -lpthread -D${MANGLING}" )
# set ( CMAKE_C_FLAGS_DEBUG   "-O0 -g -fbounds-check" )
# set ( CMAKE_C_FLAGS_RELEASE "-O3" )

# The kind of compilation : release or debug
set ( CMAKE_BUILD_TYPE      RELEASE )

# METIS library
# set ( METIS_NEEDED ON CACHE BOOL "Do you need Metis?" )
set ( METIS_PATH            "D:/bin/SFELES_CODE_MHD_parallel/metis-4.0/install" )
set ( METIS_INCLUDE_DIR     "${METIS_PATH}/include"  CACHE PATH "Path to the Metis file metis.h" )
set ( METIS_LIBRARIES       "${METIS_PATH}/lib/metis.lib" CACHE FILEPATH "Path to the Metis library file" )

# MUMPS library
set ( MUMPS_NEEDED                        TRUE )
if  ( MUMPS_NEEDED )
   set ( MUMPS_PATH                         "D:/bin/SFELES_CODE_MHD_parallel/Mumps/MUMPS_5.4.0/install" CACHE PATH "Path to MUMPS general folder" )
   set ( MUMPS_INCLUDE_DIRS                 "${MUMPS_PATH}/include" )
   set ( MUMPS_LIBRARIES                    "${MUMPS_PATH}/lib/dmumps.lib" )
   set ( MUMPS_LIBRARIES ${MUMPS_LIBRARIES} "${MUMPS_PATH}/lib/mumps_common.lib" )
   set ( MUMPS_LIBRARIES ${MUMPS_LIBRARIES} "${MUMPS_PATH}/lib/pord.lib" )
   set ( MUMPS_LIBRARIES ${MUMPS_LIBRARIES} "${MUMPS_PATH}/lib/smumps.lib" )
endif ()

# Scalapack ( needed by MUMPS ) / use either MKL scalapack or the one compiled by yourself
# if ( MUMPS_NEEDED )
  # set ( SCALAPACK_LIBRARIES   "/Users/Xavier/Documents/SFELES_CODE_MHD_parallel/Mumps/scalapack-2.0.2/libscalapack.a" CACHE FILE "Path to Scalapack library" )
# endif()

# MKL BLAS if Intel compilator ####################################################
set ( BLAS_PATH "C:/Program Files (x86)/Intel/oneAPI/mkl/latest"    CACHE PATH "Path to MKL BLAS libraries" )

# MKL BLAS if Intel compilator
set (BLA_VENDOR Intel10_64lp) 
find_package(BLAS)
#
IF (BLAS_FOUND)
  set ( BLAS_LIB "${BLAS_LIBRARIES}" CACHE FILEPATH "Set from FindBLAS.cmake BLAS_LIBRARIES." FORCE)
  message ("BLAS_LIB from BLAS_FOUND= ${BLAS_LIB}")
ELSE()

  IF (DEFINED BLAS_mkl_core_dll_LIBRARY)
    set ( BLAS_LIB "${BLAS_mkl_intel_lp64_dll_LIBRARY};${BLAS_mkl_intel_thread_dll_LIBRARY};${BLAS_mkl_core_dll_LIBRARY}" CACHE FILEPATH "Set from FindBLAS.cmake BLAS_LIBRARIES." FORCE)
    message ("BLAS_LIB = ${BLAS_LIB}")
  ENDIF()

ENDIF()

if (${CMAKE_C_COMPILER_ID} STREQUAL "Intel")
   # set ( BLAS_PATH            "/opt/intel/Compiler/11.1/080/Frameworks/mkl/lib/em64t" CACHE PATH "Path to MKL BLAS libraries" )
   # set ( BLAS_LIB             "${BLAS_PATH}/libmkl_core.a" )
   # set ( BLAS_LIB ${BLAS_LIB} "${BLAS_PATH}/libmkl_intel_lp64.a" )
   # set ( BLAS_LIB ${BLAS_LIB} "${BLAS_PATH}/libmkl_sequential.a" )
   # message ( "Intel compiler" )
   set ( HAVE_MKL TRUE )
else()
   set ( HAVE_MKL FALSE )
endif()

# MPI library #####################################################################
set ( MPI_NEEDED ON CACHE BOOL "Do you need MPI?" )
if ( MPI_NEEDED )
  # include(${CMAKE_CURRENT_SOURCE_DIR}/mpi.cmake)
  
#  set ( MPI_INCLUDE ${MPI_C_INCLUDE_DIRS} CACHE PATH "Path to MPI include folder" )

   set ( MPI_PATH     "C:/Program Files (x86)/Intel/oneAPI/mpi/latest" )
   set ( MPI_INCLUDE  "${MPI_PATH}/include" CACHE PATH "Path to MPI include folder" )
   
   list ( APPEND BLAS_LIB "${MPI_PATH}/lib/release/impi.lib" )
   
endif ()

# Scalapack ( needed by MUMPS ) ###################################################
# use either MKL scalapack or the one compiled by yourself
if ( MUMPS_NEEDED )
   
# Parallel version
  IF ( MPI_NEEDED )
    set ( SCALAPACK_LIBRARIES "${BLAS_PATH}/lib/intel64/mkl_scalapack_lp64.lib;${BLAS_PATH}/lib/intel64/mkl_blacs_intelmpi_lp64.lib" CACHE STRING "Path to the Scalapack library file")
  ELSE ()
    set ( SCALAPACK_LIBRARIES "${BLAS_PATH}/lib/intel64/mkl_scalapack_lp64.lib" CACHE STRING "Path to the Scalapack library file")
  ENDIF ()
  
#  set ( SCALAPACK_LIBRARIES   "/Users/Xavier/Documents/SFELES_CODE_MHD_parallel/Mumps/scalapack-2.0.2/libscalapack.a" CACHE FILE "Path to the Scalapack library file" )
else ()
  unset ( SCALAPACK_LIBRARIES CACHE )
endif()