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

# install
install(TARGETS ${EXEC_NAMES}
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)
