# include directories for executables
set(CMAKE_CXX_STANDARD 11)

set(PROGRAM_INCLUDES
    ${PROJECT_SOURCE_DIR}/src
    ${PROJECT_BINARY_DIR}
)

# gromos object files
set(GROMOS_OBJCECTS
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

if(MPI)
    set(EXEC_NAMES md_mpi)
else()
    set(EXEC_NAMES md)
endif()

# executables
foreach(LOOP_VAR ${EXEC_NAMES})
    add_executable(${LOOP_VAR} "${LOOP_VAR}.cc" ${GROMOS_OBJCECTS})
    target_include_directories(${LOOP_VAR} PRIVATE ${PROGRAM_INCLUDES})
    target_link_libraries(${LOOP_VAR} ${EXTERNAL_LIBRARIES})
    target_link_libraries(${LOOP_VAR} ${EXTERNAL_LIBRARIES})
endforeach()

# install
install(TARGETS ${EXEC_NAMES}
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)