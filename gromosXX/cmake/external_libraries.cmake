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

# TODO: should be implemented with find_package() macro...
if(XTB)
    message(STATUS "Taking xtb library from: ${XTB}")
    set(XTB_LIBRARIES "${XTB}/lib/x86_64-linux-gnu/libxtb.so")
    set(XTB_INCLUDES "${XTB}/include")
    add_definitions(-DWITH_XTB)
    set(EXTERNAL_LIBRARIES ${EXTERNAL_LIBRARIES} ${XTB_LIBRARIES})
    set(EXTERNAL_INCLUDES ${EXTERNAL_INCLUDES} ${XTB_INCLUDES})
endif()

# find options based libraries
if(OMP)
    find_package(FFTWomp REQUIRED)
    set(EXTERNAL_LIBRARIES ${EXTERNAL_LIBRARIES} ${FFTW_OMP_LIBRARIES})
    set(EXTERNAL_INCLUDES ${EXTERNAL_INCLUDES} ${FFTW_OMP_INCLUDE_DIRS})
endif()

if(MPI)
    find_package(FFTWmpi REQUIRED)
    set(EXTERNAL_LIBRARIES ${EXTERNAL_LIBRARIES} ${FFTW_MPI_LIBRARIES})
    set(EXTERNAL_INCLUDES ${EXTERNAL_INCLUDES}${FFTW_MPI_INCLUDE_DIRS})
endif()
