set(TEST_PROGRAM_INCLUDES
    ${PROJECT_BINARY_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/..
    ${CMAKE_CURRENT_SOURCE_DIR}/../interaction
    ${CMAKE_CURRENT_SOURCE_DIR}/../interaction/nonbonded
    ${CMAKE_CURRENT_SOURCE_DIR}/../interaction/nonbonded/pairlist
    ${CMAKE_CURRENT_SOURCE_DIR}/../interaction/nonbonded/interaction
    ${CMAKE_CURRENT_SOURCE_DIR}/../interaction/bonded
    ${CMAKE_CURRENT_SOURCE_DIR}/../interaction/forcefield
    ${CMAKE_CURRENT_SOURCE_DIR}/../interaction/special
    ${CMAKE_CURRENT_SOURCE_DIR}/../interaction/qmmm
    ${CMAKE_CURRENT_SOURCE_DIR}/../topology
    ${CMAKE_CURRENT_SOURCE_DIR}/../simulation
    ${CMAKE_CURRENT_SOURCE_DIR}/../configuration
    ${CMAKE_CURRENT_SOURCE_DIR}/../algorithm
    ${CMAKE_CURRENT_SOURCE_DIR}/../gpu
    ${CMAKE_CURRENT_SOURCE_DIR}/../gpu/cuda
    
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
set(GROMOS_LIBS
    grocuda
)

# executables
#set(EXEC_NAMES aladip aladip_unperturbed aladip_special aladip_atomic aladip_ls c16_cg lambdas)
set(EXEC_NAMES aladip aladip_unperturbed aladip_special aladip_atomic c16_cg lambdas)

add_executable(aladip ${ALADIP_CC})
add_executable(aladip_unperturbed ${ALADIP_UNPERTURBED_CC})
add_executable(aladip_special ${ALADIP_SPECIAL_CC})
add_executable(aladip_atomic ${ALADIP_ATOMIC_CC})
#add_executable(aladip_ls ${ALADIP_LS_CC})
add_executable(c16_cg ${C16_CG_CC})
add_executable(lambdas ${LAMBDAS_CC})

foreach(LOOP_VAR ${EXEC_NAMES})
    target_include_directories(${LOOP_VAR} PRIVATE ${TEST_PROGRAM_INCLUDES} ${EXTERNAL_INCLUDES}
    ${PROJECT_SOURCE_DIR}/src)
    target_link_libraries(${LOOP_VAR} PRIVATE ${GROMOS_OBJECTS} ${GROMOS_LIBS} ${EXTERNAL_LIBRARIES})
    add_test(NAME ${LOOP_VAR} COMMAND $<TARGET_FILE:${LOOP_VAR}>)
endforeach()