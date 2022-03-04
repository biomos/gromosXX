set(ALADIP_CC
    aladip.t.cc
    check_forcefield.cc
    check_state.cc
)

add_executable(aladip ${ALADIP_CC})

target_include_directories(aladip PUBLIC
    ${PROJECT_BINARY_DIR}
    .
    ..
    interaction
    interaction/nonbonded
    interaction/nonbonded/pairlist
    interaction/nonbonded/interaction
    interaction/bonded
    interaction/forcefield
    interaction/special 
    algorithm
)