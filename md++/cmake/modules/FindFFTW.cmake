# - Find the FFTW library
#
# Usage:
#   find_package(FFTW [REQUIRED] [QUIET] )
#
# It sets the following variables:
#   FFTW_FOUND ----------- true if fftw3 is found on the system
#   FFTW_LIBRARIES ------- full path to fftw3 library
#   FFTW_INCLUDE_DIRS ---- fftw3 include directory
#
# The following variables can be set to specify a search location
#   FFTW_ROOT ------------ if set, the libraries are exclusively searched under this path


#If environment variable FFTW_ROOT is specified
if(NOT FFTW_ROOT AND ENV{FFTWDIR})
    file(TO_CMAKE_PATH "$ENV{FFTW_ROOT}" FFTW_ROOT)
    set(FFTW_ROOT "${FFTW_ROOT}" CACHE PATH "Prefix for fftw3 installation.")
endif()

# Check if we can use PkgConfig
find_package(PkgConfig QUIET)
#Determine from PKG
if(PKG_CONFIG_FOUND AND NOT FFTW_ROOT)
    pkg_check_modules(PKG_FFTW QUIET "fftw3")
endif()


if(FFTW_ROOT)
    #find libs
    message( ${FFTW_ROOT} )
    find_library(
            FFTW_LIBRARIES
            NAMES "fftw3"
            PATHS ${FFTW_ROOT}
            PATH_SUFFIXES "lib" "lib64"
            NO_DEFAULT_PATH
    )
    #find includes
    find_path(
            FFTW_INCLUDE_DIRS
            NAMES "fftw3.h"
            PATHS ${FFTW_ROOT}
            PATH_SUFFIXES "include"
            NO_DEFAULT_PATH
    )
else()
    find_library(
            FFTW_LIBRARIES
            NAMES "fftw3"
            PATHS ${PKG_FFTW_LIBRARY_DIRS} ${LIB_INSTALL_DIR}
            PATH_SUFFIXES "lib" "lib64"
    )
    find_path(
            FFTW_INCLUDE_DIRS
            NAMES "fftw3.h"
            PATHS ${PKG_FFTW_INCLUDE_DIRS} ${INCLUDE_INSTALL_DIR}
    )
endif(FFTW_ROOT)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FFTW "fftw3 library could not be found. Please specify FFTW_ROOT."
        FFTW_INCLUDE_DIRS FFTW_LIBRARIES)

mark_as_advanced(FFTW_INCLUDE_DIRS FFTW_LIBRARIES)
