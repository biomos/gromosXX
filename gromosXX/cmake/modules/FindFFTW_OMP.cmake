# - Find the FFTW_OMP library
#
# Usage:
#   find_package(FFTW_OMP [REQUIRED] [QUIET] )
#
# It sets the following variables:
#   FFTW_OMP_FOUND ----------- true if fftw3_omp is found on the system
#   FFTW_OMP_LIBRARIES ------- full path to fftw3_omp library
#   FFTW_OMP_INCLUDE_DIRS ---- fftw3_omp include directory
#
# The following variables can be set to specify a search location
#   FFTW_ROOT ------------ if set, the libraries are exclusively searched under this path


#If environment variable FFTW_ROOT is specified
if(NOT FFTW_ROOT AND ENV{FFTW_OMPDIR})
    file(TO_CMAKE_PATH "$ENV{FFTW_ROOT}" FFTW_ROOT)
    set(FFTW_ROOT "${FFTW_ROOT}" CACHE PATH "Prefix for fftw3 installation.")
endif()

# Check if we can use PkgConfig
find_package(PkgConfig QUIET)
#Determine from PKG
if(PKG_CONFIG_FOUND AND NOT FFTW_ROOT)
    pkg_check_modules(PKG_FFTW_OMP QUIET "fftw3")
endif()


if(FFTW_ROOT)
    #find libs
    find_library(
            FFTW_OMP_LIBRARIES
            NAMES "fftw3_omp"
            PATHS ${FFTW_ROOT}
            PATH_SUFFIXES "lib" "lib64"
            NO_DEFAULT_PATH
    )
    #find includes
    find_path(
            FFTW_OMP_INCLUDE_DIRS
            NAMES "fftw3.h"
            PATHS ${FFTW_ROOT}
            PATH_SUFFIXES "include"
            NO_DEFAULT_PATH
    )
else()
    find_library(
            FFTW_OMP_LIBRARIES
            NAMES "fftw3_omp"
            PATHS ${PKG_FFTW_OMP_LIBRARY_DIRS} ${LIB_INSTALL_DIR}
            PATH_SUFFIXES "lib" "lib64"
    )
    find_path(
            FFTW_OMP_INCLUDE_DIRS
            NAMES "fftw3.h"
            PATHS ${PKG_FFTW_INCLUDE_DIRS} ${INCLUDE_INSTALL_DIR}
    )
endif(FFTW_ROOT)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FFTW_OMP "fftw3_omp library could not be found. Please specify FFTW_ROOT."
        FFTW_OMP_INCLUDE_DIRS FFTW_OMP_LIBRARIES)

mark_as_advanced(FFTW_OMP_INCLUDE_DIRS FFTW_OMP_LIBRARIES)
