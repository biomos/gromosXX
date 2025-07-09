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

if(USE_CUDA)
    find_package(CUDAToolkit REQUIRED)
    list(APPEND GROMOS_OBJECTS grocuda) # Link the static library instead of object files
endif()

if(MPI)
    set(EXEC_NAMES md_mpi repex_mpi eds_2box_mpi)
else()
    set(EXEC_NAMES md eds_2box)
endif()

# executables
foreach(LOOP_VAR ${EXEC_NAMES})
    add_executable(${LOOP_VAR} "${LOOP_VAR}.cc")
    target_link_libraries(${LOOP_VAR} PRIVATE ${GROMOS_OBJECTS} ${GROMOS_LIBS} ${EXTERNAL_LIBRARIES})
    target_include_directories(${LOOP_VAR} PRIVATE ${PROGRAM_INCLUDES} ${EXTERNAL_INCLUDES}
    ${PROJECT_SOURCE_DIR}/src)
endforeach()

# install
install(TARGETS ${EXEC_NAMES}
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)
