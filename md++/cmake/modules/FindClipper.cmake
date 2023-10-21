
# - Find the CLIPPER library
#
# Usage:
#   find_package(CLIPPER [REQUIRED] [QUIET] )
#
# It sets the following variables:
#   CLIPPER_FOUND ----------- true if clipper is found on the system
#   CLIPPER_LIBRARIES ------- full path to clipper
#   CLIPPER_INCLUDE_DIRS ---- clipper include directory
#
# The following variables can be set to specify a search location
#   CLIPPER_ROOT ------------ if set, the libraries are exclusively searched under this path


#If environment variable CLIPPER_ROOT is specified, it has same effect as CLIPPER_ROOT
if(NOT CLIPPER_ROOT AND ENV{CLIPPER_ROOT})
    file(TO_CMAKE_PATH "$ENV{GSL_ROOT}" CLIPPER_ROOT)
    set(GSL_ROOT "${CLIPPER_ROOT}" CACHE PATH "Prefix for Clipper installation.")
endif()

# Check if we can use PkgConfig
find_package(PkgConfig QUIET)
#Determine from PKG
if(PKG_CONFIG_FOUND AND NOT CLIPPER_ROOT)
    pkg_check_modules(PKG_FFTW QUIET "clipper" QUIET)
endif()


if(CLIPPER_ROOT)
    #find libs
    find_library(
            CLIPPER_CCP4_LIB
            NAMES "clipper-ccp4"
            PATHS ${CLIPPER_ROOT}
            PATH_SUFFIXES "lib" "lib64"
            NO_DEFAULT_PATH
    )
    find_library(
            CLIPPER_CCP4C_LIB
            NAMES "ccp4c"
            PATHS ${CLIPPER_ROOT}
            PATH_SUFFIXES "lib" "lib64"
            NO_DEFAULT_PATH
    )
    find_library(
            CLIPPER_CONTRIB_LIB
            NAMES "clipper-contrib"
            PATHS ${CLIPPER_ROOT}
            PATH_SUFFIXES "lib" "lib64"
            NO_DEFAULT_PATH
    )
    find_library(
            CLIPPER_CORE_LIB
            NAMES "clipper-core"
            PATHS ${CLIPPER_ROOT}
            PATH_SUFFIXES "lib" "lib64"
            NO_DEFAULT_PATH
    )
    #find includes
    find_path(
            CLIPPER_INCLUDE_DIRS
            NAMES "clipper.h"
            PATHS ${CLIPPER_ROOT}
            PATH_SUFFIXES "include clipper include/clipper"
            NO_DEFAULT_PATH
    )
else()
    #find libs
    find_library(
            CLIPPER_CCP4_LIB
            NAMES "clipper-ccp4"
            PATH_SUFFIXES "lib" "lib64"
    )
    find_library(
            CLIPPER_CCP4C_LIB
            NAMES "ccp4c"
            PATH_SUFFIXES "lib" "lib64"
    )
    find_library(
            CLIPPER_CONTRIB_LIB
            NAMES "clipper-contrib"
            PATH_SUFFIXES "lib" "lib64"
    )
    find_library(
            CLIPPER_CORE_LIB
            NAMES "clipper-core"
            PATH_SUFFIXES "lib" "lib64"
    )
    #find includes
    find_path(
            CLIPPER_INCLUDE_DIRS
            NAMES "clipper.h"
            PATH_SUFFIXES "include" "clipper" "include/clipper"
    )
endif()


include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Clipper "Clipper library could not be found. Please specify CLIPPER_ROOT."
        CLIPPER_INCLUDE_DIRS CLIPPER_CCP4_LIB CLIPPER_CCP4C_LIB CLIPPER_CONTRIB_LIB CLIPPER_CORE_LIB)

set(CLIPPER_LIBRARIES
        ${CLIPPER_CCP4_LIB}
        ${CLIPPER_CCP4C_LIB}
        ${CLIPPER_CONTRIB_LIB}
        ${CLIPPER_CORE_LIB}
        )

mark_as_advanced(CLIPPER_INCLUDE_DIRS CLIPPER_LIBRARIES)

