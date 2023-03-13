# - Find the FFTW_MPI library
#
# Usage:
#   find_package(FFTW_MPI [REQUIRED] [QUIET] )
#
# It sets the following variables:
#   FFTW_MPI_FOUND ----------- true if fftw3_mpi is found on the system
#   FFTW_MPI_LIBRARIES ------- full path to fftw3_mpi library
#   FFTW_MPI_INCULDE_DIRS ---- fftw3_mpi include directory
#
# The following variables can be set to specify a search location
#   FFTW_ROOT ------------ if set, the libraries are exclusively searched under this path


#If environment variable FFTW_ROOT is specified
if(NOT FFTW_ROOT AND ENV{FFTWDIR})
    file(TO_CMAKE_PATH "$ENV{FFTW_ROOT}" FFTW_ROOT)
    set(FFTW_ROOT "${FFTW_ROOT}" CACHE PATH "Prefix for fftw3_mpi installation.")
endif()

# Check if we can use PkgConfig
find_package(PkgConfig QUIET)
#Determine from PKG
if(PKG_CONFIG_FOUND AND NOT FFTW_ROOT)
    pkg_check_modules(PKG_FFTW QUIET "fftw3_mpi")
endif()


if(FFTW_ROOT)
    #find libs
    find_library(
            FFTW_MPI_LIBRARIES
            NAMES "fftw3_mpi"
            PATHS ${FFTW_ROOT}
            PATH_SUFFIXES "lib" "lib64"
            NO_DEFAULT_PATH
    )
    #find includes
    find_path(
            FFTW_MPI_INCULDE_DIRS
            NAMES "fftw3-mpi.h"
            PATHS ${FFTW_ROOT}
            PATH_SUFFIXES "include"
            NO_DEFAULT_PATH
    )
else()
    find_library(
            FFTW_MPI_LIBRARIES
            NAMES "fftw3_mpi"
            PATHS ${PKG_FFTW_LIBRARY_DIRS} ${LIB_INSTALL_DIR}
            PATH_SUFFIXES "lib" "lib64"
    )
    find_path(
            FFTW_MPI_INCULDE_DIRS
            NAMES "fftw3-mpi.h"
            PATHS ${PKG_FFTW_INCLUDE_DIRS} ${INCLUDE_INSTALL_DIR}
    )
endif(FFTW_ROOT)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FFTW_MPI "fftw3_mpi library could not be found. Please specify FFTW_ROOT."
        FFTW_MPI_INCULDE_DIRS FFTW_MPI_LIBRARIES)

mark_as_advanced(FFTW_MPI_INCULDE_DIRS FFTW_MPI_LIBRARIES)
