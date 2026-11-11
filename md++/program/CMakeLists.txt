# include directories for executables

set(PROGRAM_INCLUDES
    ${PROJECT_SOURCE_DIR}/src
    ${PROJECT_BINARY_DIR}
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

if(MPI)
    set(EXEC_NAMES md_mpi repex_mpi eds_2box_mpi)
else()
    set(EXEC_NAMES md eds_2box)
endif()

# executables
foreach(LOOP_VAR ${EXEC_NAMES})
    add_executable(${LOOP_VAR} "${LOOP_VAR}.cc" ${GROMOS_OBJECTS})
    target_include_directories(${LOOP_VAR} PRIVATE ${PROGRAM_INCLUDES} ${EXTERNAL_INCLUDES})
    target_link_libraries(${LOOP_VAR} PRIVATE ${EXTERNAL_LIBRARIES})
endforeach()

if(COLVAR_TESTS)
    add_executable(colvar_extti_check
        "${PROJECT_SOURCE_DIR}/contrib/colvar_tests/extti_check.cc" ${GROMOS_OBJECTS})
    target_include_directories(colvar_extti_check PRIVATE ${PROGRAM_INCLUDES} ${EXTERNAL_INCLUDES})
    target_link_libraries(colvar_extti_check PRIVATE ${EXTERNAL_LIBRARIES})
    add_executable(colvar_le_check
        "${PROJECT_SOURCE_DIR}/contrib/colvar_tests/le_check.cc" ${GROMOS_OBJECTS})
    target_include_directories(colvar_le_check PRIVATE ${PROGRAM_INCLUDES} ${EXTERNAL_INCLUDES})
    target_link_libraries(colvar_le_check PRIVATE ${EXTERNAL_LIBRARIES})
    add_executable(colvar_eds_check
        "${PROJECT_SOURCE_DIR}/contrib/colvar_tests/eds_check.cc" ${GROMOS_OBJECTS})
    target_include_directories(colvar_eds_check PRIVATE ${PROGRAM_INCLUDES} ${EXTERNAL_INCLUDES})
    target_link_libraries(colvar_eds_check PRIVATE ${EXTERNAL_LIBRARIES})
    add_executable(colvar_force_check
        "${PROJECT_SOURCE_DIR}/contrib/colvar_tests/force_check.cc" ${GROMOS_OBJECTS})
    target_include_directories(colvar_force_check PRIVATE ${PROGRAM_INCLUDES} ${EXTERNAL_INCLUDES})
    target_link_libraries(colvar_force_check PRIVATE ${EXTERNAL_LIBRARIES})
    add_executable(colvar_angular_bias_check
        "${PROJECT_SOURCE_DIR}/contrib/colvar_tests/angular_bias_check.cc" ${GROMOS_OBJECTS})
    target_include_directories(colvar_angular_bias_check PRIVATE ${PROGRAM_INCLUDES} ${EXTERNAL_INCLUDES})
    target_link_libraries(colvar_angular_bias_check PRIVATE ${EXTERNAL_LIBRARIES})
    add_executable(colvar_average_check
        "${PROJECT_SOURCE_DIR}/contrib/colvar_tests/average_check.cc" ${GROMOS_OBJECTS})
    target_include_directories(colvar_average_check PRIVATE ${PROGRAM_INCLUDES} ${EXTERNAL_INCLUDES})
    target_link_libraries(colvar_average_check PRIVATE ${EXTERNAL_LIBRARIES})
    add_executable(colvar_energy_check
        "${PROJECT_SOURCE_DIR}/contrib/colvar_tests/energy_check.cc" ${GROMOS_OBJECTS})
    target_include_directories(colvar_energy_check PRIVATE ${PROGRAM_INCLUDES} ${EXTERNAL_INCLUDES})
    target_link_libraries(colvar_energy_check PRIVATE ${EXTERNAL_LIBRARIES})
endif()


# Full MD regression for CV geometry isolation and both SHAKE paths. Optional
# numerical checks already live under COLVAR_TESTS; Python is needed only here.
if(COLVAR_TESTS AND TARGET md)
    find_package(Python3 COMPONENTS Interpreter QUIET)
    if(Python3_Interpreter_FOUND)
        add_test(NAME colvar_angular_shake
            COMMAND ${Python3_EXECUTABLE}
                ${PROJECT_SOURCE_DIR}/contrib/colvar_tests/run_angular_shake.py
                --md $<TARGET_FILE:md>
                --reference ${PROJECT_SOURCE_DIR}/src/check/data)
        add_test(NAME colvar_extended_ti
            COMMAND ${Python3_EXECUTABLE}
                ${PROJECT_SOURCE_DIR}/contrib/colvar_tests/run_extti.py
                --md $<TARGET_FILE:md> --checker $<TARGET_FILE:colvar_extti_check>
                --reference ${PROJECT_SOURCE_DIR}/src/check/data)
        set_tests_properties(colvar_extended_ti PROPERTIES
            ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1"
            TIMEOUT 180)
        set_tests_properties(colvar_angular_shake PROPERTIES
            ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1"
            TIMEOUT 180)
    endif()
endif()

# install
install(TARGETS ${EXEC_NAMES}
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)
