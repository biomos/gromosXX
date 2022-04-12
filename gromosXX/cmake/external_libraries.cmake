# use backported find module for Zlib library if version < 3.0
if(${CMAKE_VERSION} VERSION_LESS "3.0.0")
    set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake/modules/zlib")
endif()

# use backported find module for GSL library if version < 3.2
if(${CMAKE_VERSION} VERSION_LESS "3.2.0")
    set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake/modules/gsl")
endif()

# find always required libraries
find_package(FFTW REQUIRED)
find_package(Threads REQUIRED)
find_package(GSL REQUIRED)
find_package(ZLIB REQUIRED)
set(EXTERNAL_LIBRARIES
    ${EXTERNAL_LIBRARIES}
    ${CMAKE_THREAD_LIBS_INIT}
    ${GSL_LIBRARIES}
    ${FFTW_LIBRARIES}
    ${ZLIB_LIBRARIES}
    xtb
)

# xtb
include_directories(/home/fpultar/xtb/include)
link_directories(/home/fpultar/src/xtb/build)

# find optional libraries
find_package(Clipper QUIET)
if(CLIPPER_FOUND)
    message(STATUS "Clipper library found, enabling usage.")
    set(EXTERNAL_LIBRARIES ${EXTERNAL_LIBRARIES} ${CLIPPER_LIBRARIES})
    include_directories(${CLIPPER_INCLUDE_DIRS})
else()
    message(STATUS "Clipper usage disabled.")
endif()

# find options based libraries
if(OMP)
    find_package(FFTWomp REQUIRED)
    set(EXTERNAL_LIBRARIES ${EXTERNAL_LIBRARIES} ${FFTW_OMP_LIBRARIES})
    include_directories(${FFTW_OMP_INCLUDE_DIRS})
endif()

if(MPI)
    find_package(FFTWmpi REQUIRED)
    set(EXTERNAL_LIBRARIES ${EXTERNAL_LIBRARIES} ${FFTW_MPI_LIBRARIES})
    include_directories(${FFTW_MPI_INCLUDE_DIRS})
endif()
