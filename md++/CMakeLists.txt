cmake_minimum_required(VERSION 3.13)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(USE_CUDA)
    project(gromosXX LANGUAGES CXX CUDA)
    include_directories(${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
else()
    project(gromosXX LANGUAGES CXX)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# get gromos version number
file(READ "VERSION" VERSION_FILE)
string(REGEX MATCH MAJOR.* MD_MAJOR_VERSION "${VERSION_FILE}")
string(REGEX MATCH [0-9]+ MD_MAJOR_VERSION "${MD_MAJOR_VERSION}")
string(REGEX MATCH MINOR.* MD_MINOR_VERSION "${VERSION_FILE}")
string(REGEX MATCH [0-9]+ MD_MINOR_VERSION "${MD_MINOR_VERSION}")
string(REGEX MATCH MICRO.* MD_MICRO_VERSION "${VERSION_FILE}")
string(REGEX MATCH [0-9]+ MD_MICRO_VERSION "${MD_MICRO_VERSION}")

# set module search path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake/modules/")

# check folder
enable_testing()

# options
include(${PROJECT_SOURCE_DIR}/cmake/options.cmake)
# find external libraries based on options
include(${PROJECT_SOURCE_DIR}/cmake/external_libraries.cmake)
# perform tests and configure according to results / options
include(${PROJECT_SOURCE_DIR}/cmake/configure.cmake)

# code
add_subdirectory(src)

# remove duplicate includes and libraries
list(REMOVE_DUPLICATES EXTERNAL_INCLUDES)
list(REMOVE_DUPLICATES EXTERNAL_LIBRARIES)

# md, md_mpi, eds, repex_mpi, ...
add_subdirectory(program)
