set(SIM_CC
    simulation.cc
    mpiControl.cpp
    multibath.cc
)

add_library(grosimulation OBJECT ${SIM_CC})

target_include_directories(grosimulation PUBLIC
    ..
    ${PROJECT_BINARY_DIR}
    ${EXTERNAL_INCLUDES}
)