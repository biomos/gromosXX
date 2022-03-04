
# - Find the GSL library
#
# Usage:
#   find_package(GSL [REQUIRED] [QUIET] )
#
# It sets the following variables:
#   GSL_FOUND ----------- true if gsl is found on the system
#   GSL_LIBRARIES ------- full path to gsl and cblas library
#   GSL_LIBRARY --------- full path to gsl library only
#   GSL_CBLAS_LIBRARY --- full path to cblas library only
#   GSL_INCLUDE_DIRS ---- gsl include directory
#
# The following variables can be set to specify a search location
#   GSL_ROOT ------------ if set, the libraries are exclusively searched under this path

#If environment variable GSL_ROOT is specified
if(NOT GSL_ROOT AND ENV{GSL_ROOT})
    file(TO_CMAKE_PATH "$ENV{GSL_ROOT}" GSL_ROOT)
    set(GSL_ROOT "${GSL_ROOT}" CACHE PATH "Prefix for GSL installation.")
endif()

# Check if we can use PkgConfig
find_package(PkgConfig QUIET)

#Determine from PKG
if(PKG_CONFIG_FOUND AND NOT GSL_ROOT)
    pkg_check_modules(PKG_GSL QUIET "gsl")
endif()

if(GSL_ROOT)
    #find libs
    find_library(
            GSL_LIBRARY
            NAMES "gsl"
            PATHS ${GSL_ROOT} ${GSL_ROOT}/lib
            PATH_SUFFIXES "Release" "Debug"
            NO_DEFAULT_PATH
    )
    find_library(
            GSL_CBLAS_LIBRARY
            NAMES "gslcblas" "cblas"
            PATHS ${GSL_ROOT} ${GSL_ROOT}/lib
            PATH_SUFFIXES "Release" "Debug"
            NO_DEFAULT_PATH
    )
    #find includes
    find_path(
            GSL_INCLUDE_DIRS
            NAMES "gsl_sf.h"
            PATHS ${GSL_ROOT} ${GSL_ROOT}/include
            PATH_SUFFIXES "gsl"
            NO_DEFAULT_PATH
    )
else()
    #find libs
    find_library(
            GSL_LIBRARY
            NAMES "gsl"
            PATHS ${PKG_GSL_LIBRARY_DIRS} ${LIB_INSTALL_DIR}
            PATH_SUFFIXES "Release" "Debug"
    )
    find_library(
            GSL_CBLAS_LIBRARY
            NAMES "gslcblas" "cblas"
            PATHS ${PKG_GSL_LIBRARY_DIRS} ${LIB_INSTALL_DIR}
            PATH_SUFFIXES "Release" "Debug"
    )
    #find includes
    find_path(
            GSL_INCLUDE_DIRS
            NAMES "gsl_sf.h"
            PATHS ${PKG_GSL_INCLUDE_DIRS} ${INCLUDE_INSTALL_DIR}
            PATH_SUFFIXES "gsl"
    )
endif(GSL_ROOT)

# check if found, handes REQUIRED, set GSL_FOUND
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GSL "GSL library could not be found. Please specify GSL_ROOT."
        GSL_INCLUDE_DIRS GSL_LIBRARY GSL_CBLAS_LIBRARY)

set(GSL_LIBRARIES ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY})
mark_as_advanced(GSL_INCLUDE_DIRS GSL_LIBRARIES GSL_LIBRARY GSL_CBLAS_LIBRARY)
