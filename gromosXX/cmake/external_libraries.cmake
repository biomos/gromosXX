# find always required libraries
find_package(FFTW REQUIRED)
find_package(Threads REQUIRED)
find_package(GSL REQUIRED)
find_package(ZLIB REQUIRED)

# external includes (excluding special options like xtb -> they are added below)
# included by all targets
set(EXTERNAL_INCLUDES
    ${EXTERNAL_INCLUDES}
    ${GSL_INCLUDE_DIRS}
    ${ZLIB_INCLUDE_DIRS}
    ${FFTW_INCLUDE_DIRS}
)

# external libraries (excluding special options like xtb -> they are added below)
# included by all targets
set(EXTERNAL_LIBRARIES
    ${EXTERNAL_LIBRARIES}
    ${CMAKE_THREAD_LIBS_INIT}
    ${GSL_LIBRARIES}
    ${FFTW_LIBRARIES}
    ${ZLIB_LIBRARIES}
)

# find optional libraries
find_package(Clipper QUIET)
if(CLIPPER_FOUND)
    message(STATUS "Clipper library found, enabling usage.")
    set(EXTERNAL_LIBRARIES ${EXTERNAL_LIBRARIES} ${CLIPPER_LIBRARIES})
    set(EXTERNAL_INCLUDES ${EXTERNAL_INCLUDES} ${CLIPPER_INCLUDE_DIRS})
else()
    message(STATUS "Clipper usage disabled.")
endif()

if(XTB)
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(XTB REQUIRED IMPORTED_TARGET xtb)
    add_definitions(-DWITH_XTB)
    set(EXTERNAL_LIBRARIES ${EXTERNAL_LIBRARIES} ${XTB_LIBRARIES})
    set(EXTERNAL_INCLUDES ${EXTERNAL_INCLUDES} ${XTB_INCLUDE_DIRS})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${XTB_CFLAGS_OTHER}")
endif()

# find options based libraries
if(OMP)
    find_package(FFTW_OMP REQUIRED)
    set(EXTERNAL_LIBRARIES ${EXTERNAL_LIBRARIES} ${FFTW_OMP_LIBRARIES})
    set(EXTERNAL_INCLUDES ${EXTERNAL_INCLUDES} ${FFTW_OMP_INCLUDE_DIRS})
endif()

if(MPI)
    find_package(FFTW_MPI REQUIRED)
    set(EXTERNAL_LIBRARIES ${EXTERNAL_LIBRARIES} ${FFTW_MPI_LIBRARIES})
    set(EXTERNAL_INCLUDES ${EXTERNAL_INCLUDES}${FFTW_MPI_INCLUDE_DIRS})
endif()
