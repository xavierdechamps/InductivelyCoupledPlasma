message ( "Entering in the configuration file" )

# General options / Flags #########################################################
unset ( Windows       CACHE )
set ( Windows           ON CACHE BOOL "Is it a Windows-based environment?" )

# Fortran compiler
# set ( CMAKE_Fortran_COMPILER      "mpif90" )
# set ( CMAKE_Fortran_FLAGS         "-w -fno-range-check -m64 -m64" )
# set ( CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -fbounds-check" )
# set ( CMAKE_Fortran_FLAGS_RELEASE "-O3" )
set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} /heap-arrays:0 /QxSSE4.2")
set ( CMAKE_Fortran_FLAGS_RELEASE "/O3" CACHE STRING "Fortran options release" )
set ( CMAKE_EXE_LINKER_FLAGS "/NODEFAULTLIB:library" )

# The kind of compilation : release or debug
set ( CMAKE_BUILD_TYPE      RELEASE )

# METIS library ##################################################################
# Required by MUMPS
set ( METIS_PATH            "D:/bin/SFELES_CODE_MHD_parallel/metis-4.0/install" )
set ( METIS_INCLUDE_DIR     "${METIS_PATH}/include"  CACHE PATH "Path to the Metis file metis.h" )
set ( METIS_LIBRARIES       "${METIS_PATH}/lib/metis.lib" CACHE FILEPATH "Path to the Metis library file" )

# MUMPS library ##################################################################
set ( MUMPS_PATH                         "D:/bin/SFELES_CODE_MHD_parallel/Mumps/MUMPS_5.4.0/install" CACHE PATH "Path to MUMPS general folder" )
set ( MUMPS_INCLUDE_DIRS                 "${MUMPS_PATH}/include" )
set ( MUMPS_LIBRARIES                    "${MUMPS_PATH}/lib/dmumps.lib" )
set ( MUMPS_LIBRARIES ${MUMPS_LIBRARIES} "${MUMPS_PATH}/lib/mumps_common.lib" )
set ( MUMPS_LIBRARIES ${MUMPS_LIBRARIES} "${MUMPS_PATH}/lib/pord.lib" )
set ( MUMPS_LIBRARIES ${MUMPS_LIBRARIES} "${MUMPS_PATH}/lib/smumps.lib" )

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

set ( MPI_PATH     "C:/Program Files (x86)/Intel/oneAPI/mpi/latest" )
set ( MPI_INCLUDE  "${MPI_PATH}/include" CACHE PATH "Path to MPI include folder" )
list ( APPEND BLAS_LIB "${MPI_PATH}/lib/release/impi.lib" )

# Scalapack ( needed by MUMPS ) ###################################################
set ( SCALAPACK_LIBRARIES "${BLAS_PATH}/lib/intel64/mkl_scalapack_lp64.lib;${BLAS_PATH}/lib/intel64/mkl_blacs_intelmpi_lp64.lib" CACHE STRING "Path to the Scalapack library file")
