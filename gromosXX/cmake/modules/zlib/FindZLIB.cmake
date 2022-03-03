
# - Find the ZLIB library
#
# Usage:
#   find_package(ZLIB [REQUIRED] [QUIET] )
#
# It sets the following variables:
#   ZLIB_FOUND ----------- true if zlib is found on the system
#   ZLIB_LIBRARIES ------- full path to zlib library
#   ZLIB_INCLUDE_DIRS ---- zlib include directory
#
# The following variables can be set to specify a search location
#   ZLIB_ROOT -------- if set, the libraries are exclusively searched under this path

#If environment variable ZLIB_ROOT is specified, it has same effect as ZLIB_ROOT
if(NOT ZLIB_ROOT AND ENV{ZLIB_ROOT})
    file(TO_CMAKE_PATH "$ENV{ZLIB_ROOT}" ZLIB_ROOT)
    set(ZLIB_ROOT "${ZLIB_ROOT}" CACHE PATH "Prefix for ZLIB installation.")
endif()

# Check if we can use PkgConfig
find_package(PkgConfig QUIET)

#Determine from PKG
if(PKG_CONFIG_FOUND AND NOT ZLIB_ROOT)
    pkg_check_modules(PKG_ZLIB QUIET "zlib")
endif()

if(ZLIB_ROOT)
    #find libs
    find_library(
            ZLIB_LIBRARIES
            NAMES z zlib zdll zlib1 zlibd zlibd1
            PATHS ${ZLIB_ROOT}
            PATH_SUFFIXES "lib"
            NO_DEFAULT_PATH
    )
    #find includes
    find_path(
            ZLIB_INCLUDE_DIRS
            NAMES "zlib.h"
            PATHS ${ZLIB_ROOT}
            PATH_SUFFIXES "include"
            NO_DEFAULT_PATH
    )
else()
    #find libs
    find_library(
            ZLIB_LIBRARIES
            NAMES z zlib zdll zlib1 zlibd zlibd1
            PATHS ${PKG_ZLIB_LIBRARY_DIRS} ${LIB_INSTALL_DIR}
            PATH_SUFFIXES "lib"
    )
    #find includes
    find_path(
            ZLIB_INCLUDE_DIRS
            NAMES "zlib.h"
            PATHS ${PKG_ZLIB_INCLUDE_DIRS} ${INCLUDE_INSTALL_DIR}
            PATH_SUFFIXES "include"
    )
endif(ZLIB_ROOT)

# check if found, handes REQUIRED, set ZLIB_FOUND
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ZLIB "ZLIB library could not be found. Please specify ZLIB_ROOT."
        ZLIB_INCLUDE_DIRS ZLIB_LIBRARIES)

mark_as_advanced(ZLIB_INCLUDE_DIRS ZLIB_LIBRARIES)
