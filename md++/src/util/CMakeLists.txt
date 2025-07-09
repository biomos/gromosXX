set(UTIL_CC
    cycle_thread.cc
    bs_coordinate.cc
    bs_potentials.cc
    bs_subspace.cc
    bs_umbrella.cc
    bs_vector.cc
    debug.cc
    timing.cc
    usage.cc
    generate_velocities.cc
    parse_tcouple.cc
    parse_verbosity.cc
    prepare_virial.cc
    create_simulation.cc
    virtual_atom.cc
    coding.cc
    feature.cc
    le_coordinate.cc
    umbrella.cc
    umbrella_weight.cc
    thread.cc
    system_call.cc
    undirected_graph.cc
)

add_library(groutil OBJECT ${UTIL_CC})

target_include_directories(groutil PUBLIC
    ..
    ${PROJECT_BINARY_DIR}
    ${EXTERNAL_INCLUDES}
    ${PROJECT_BINARY_DIR}/src/util # config.h included as "../../config.h"
    ${PROJECT_SOURCE_DIR}/src
)
