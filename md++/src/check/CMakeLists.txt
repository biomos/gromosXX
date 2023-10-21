set(TEST_PROGRAM_INCLUDES
    ${PROJECT_BINARY_DIR}
    .
    ..
    ../interaction
    ../interaction/nonbonded
    ../interaction/nonbonded/pairlist
    ../interaction/nonbonded/interaction
    ../interaction/bonded
    ../interaction/forcefield
    ../interaction/special
    ../interaction/qmmm
    ../topology
    ../simulation
    ../configuration
    ../algorithm
)

set(ALADIP_CC
    aladip.t.cc
    check_forcefield.cc
    check_state.cc
)

set(ALADIP_UNPERTURBED_CC
    aladip_unperturbed.t.cc
    check_forcefield.cc
    check_state.cc 
)

set(ALADIP_SPECIAL_CC
    aladip_special.t.cc
    check_forcefield.cc
)

set(ALADIP_ATOMIC_CC
    aladip_atomic.t.cc
    check_forcefield.cc
    check_state.cc
)

set(ALADIP_LS_CC
    aladip_ls.t.cc
    check_forcefield.cc
    check_state.cc
)

set(C16_CG_CC
    c16_cg.t.cc
    check_forcefield.cc
)

set(LAMBDAS_CC
    lambdas.t.cc
	check_forcefield.cc
	check_state.cc
)

# gromos object files
set(GROMOS_OBJECTS
    $<TARGET_OBJECTS:groalgorithm>
    $<TARGET_OBJECTS:groconfiguration>
    $<TARGET_OBJECTS:groio>
    $<TARGET_OBJECTS:grointeraction>
    $<TARGET_OBJECTS:groutil>
    $<TARGET_OBJECTS:grotopology>
    $<TARGET_OBJECTS:grosimulation>
    $<TARGET_OBJECTS:gromath>
    $<TARGET_OBJECTS:grorepex>
)

if(CUKERNEL)
    list(APPEND GROMOS_OBJECTS $<TARGET_OBJECTS:grocukernel>)
endif()

# executables
set(EXEC_NAMES aladip aladip_unperturbed aladip_special aladip_atomic aladip_ls c16_cg lambdas)

add_executable(aladip ${ALADIP_CC} ${GROMOS_OBJECTS})
add_executable(aladip_unperturbed ${ALADIP_UNPERTURBED_CC} ${GROMOS_OBJECTS})
add_executable(aladip_special ${ALADIP_SPECIAL_CC} ${GROMOS_OBJECTS})
add_executable(aladip_atomic ${ALADIP_ATOMIC_CC} ${GROMOS_OBJECTS})
add_executable(aladip_ls ${ALADIP_LS_CC} ${GROMOS_OBJECTS})

add_executable(c16_cg ${C16_CG_CC} ${GROMOS_OBJECTS})
add_executable(lambdas ${LAMBDAS_CC} ${GROMOS_OBJECTS})

# links and includes for tests
foreach(LOOP_VAR ${EXEC_NAMES})
    target_include_directories(${LOOP_VAR} PRIVATE ${TEST_PROGRAM_INCLUDES} ${EXTERNAL_INCLUDES})
    target_link_libraries(${LOOP_VAR} ${EXTERNAL_LIBRARIES})
    add_test(NAME ${LOOP_VAR} COMMAND $<TARGET_FILE:${LOOP_VAR}>)
endforeach()
