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

set(ALADIP_CUDA_CC
    aladip_cuda.t.cc
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

if(USE_CUDA)
    list(APPEND EXEC_NAMES aladip_cuda)
    add_executable(aladip_cuda ${ALADIP_CUDA_CC})

    list(APPEND EXEC_NAMES pairlist_cuda_equivalence)
    add_executable(pairlist_cuda_equivalence pairlist_cuda_equivalence.t.cc)

    list(APPEND EXEC_NAMES lj_crf_tile_kernel)
    add_executable(lj_crf_tile_kernel lj_crf_tile_kernel.t.cc)

    list(APPEND EXEC_NAMES cuda_nonbonded_interaction)
    add_executable(cuda_nonbonded_interaction cuda_nonbonded_interaction.t.cc)

    list(APPEND EXEC_NAMES cuda_skin_drift)
    add_executable(cuda_skin_drift cuda_skin_drift.t.cc)

    list(APPEND EXEC_NAMES leap_frog_gpu)
    add_executable(leap_frog_gpu leap_frog_gpu.t.cc)

    list(APPEND EXEC_NAMES cuda_atomic_cutoff_toggle)
    add_executable(cuda_atomic_cutoff_toggle cuda_atomic_cutoff_toggle.t.cc)

    list(APPEND EXEC_NAMES remove_com_motion_gpu)
    add_executable(remove_com_motion_gpu remove_com_motion_gpu.t.cc)

    list(APPEND EXEC_NAMES temperature_gpu)
    add_executable(temperature_gpu temperature_gpu.t.cc)

    list(APPEND EXEC_NAMES quartic_bond_gpu)
    add_executable(quartic_bond_gpu quartic_bond_gpu.t.cc)

    list(APPEND EXEC_NAMES angle_gpu)
    add_executable(angle_gpu angle_gpu.t.cc)

    list(APPEND EXEC_NAMES improper_dihedral_gpu)
    add_executable(improper_dihedral_gpu improper_dihedral_gpu.t.cc)

    list(APPEND EXEC_NAMES dihedral_gpu)
    add_executable(dihedral_gpu dihedral_gpu.t.cc)

    list(APPEND EXEC_NAMES shake_gpu)
    add_executable(shake_gpu shake_gpu.t.cc)

    list(APPEND EXEC_NAMES pressure_barostat_gpu)
    add_executable(pressure_barostat_gpu pressure_barostat_gpu.t.cc)

    list(APPEND EXEC_NAMES solute_shake_gpu)
    add_executable(solute_shake_gpu solute_shake_gpu.t.cc)

    list(APPEND EXEC_NAMES settle_gpu)
    add_executable(settle_gpu settle_gpu.t.cc)

    list(APPEND EXEC_NAMES lincs_gpu)
    add_executable(lincs_gpu lincs_gpu.t.cc)
endif()

foreach(LOOP_VAR ${EXEC_NAMES})
    target_include_directories(${LOOP_VAR} PRIVATE ${TEST_PROGRAM_INCLUDES} ${EXTERNAL_INCLUDES}
    ${PROJECT_SOURCE_DIR}/src)
    target_link_libraries(${LOOP_VAR} PRIVATE ${GROMOS_OBJECTS} ${GROMOS_LIBS} ${EXTERNAL_LIBRARIES})
    add_test(NAME ${LOOP_VAR} COMMAND $<TARGET_FILE:${LOOP_VAR}>)
endforeach()