# set default build type to release for single generators
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING
            "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif()

# set options
option(OMP "enable OMP" OFF)
option(MPI "enable MPI" OFF)
option(USE_CUDA "enable CUDA" OFF)
option(XTB "enable XTB" OFF)

option(BUILD_SHARED_LIBS "Build shared libraries instead of static" ON)

# TODO do we still want to support these options?
option(FORCEGROUPS "enable forcegroups" OFF)
option(HEAVISIDE "enable heaviside" OFF)

# more pedantic compiler checks
option(PEDANTIC "enable pedantic compiler checks")

set(CMAKE_CXX_FLAGS "-Wall")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

if(PEDANTIC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra -Wshadow -Wnon-virtual-dtor -pedantic -Wunused -Wconversion -Wsign-conversion")
endif()

if(OMP AND MPI)
    message(FATAL_ERROR "OMP and MPI must NOT be enabled at the same time")
endif()

if(OMP)
    # if(CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    #     message(FATAL_ERROR "Apple Clang does not support OpenMP. Use Clang from Homebrew or GCC.")
    # endif()
    find_package(OpenMP REQUIRED)
    add_definitions(-DOMP)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

if(MPI)
    find_package(MPI REQUIRED)
    add_definitions(-DXXMPI)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS}")
    list(APPEND EXTERNAL_LIBRARIES ${MPI_CXX_LIBRARIES})
    list(APPEND EXTERNAL_INCLUDES ${MPI_CXX_INCLUDE_PATH})
endif()
