# Add GPU source files
set(GPU_SOURCES
    manager/cuda_manager.cc
)

# Add GPU header files (for clarity, not compiled directly)
set(GPU_HEADERS
    utils.h
    manager/cuda_memory_manager.h
    manager/cuda_device_worker.h
    manager/cuda_manager.h
    manager/cuda_device_manager.h
    kernels/reduction.h
    kernels/indexing.h
    memory/constants.h
    memory/cuvector.h
    memory/cumallocator.h
    memory/container.h
    algorithm/constraints.h
)

# If CUDA is enabled, include .cu files
if(USE_CUDA)
    enable_language(CUDA)
    find_package(CUDAToolkit REQUIRED)

    set(GPU_CUDA_SOURCES
        kernels/reduction.cu
        memory/container.cu
        manager/cuda_device_worker.cu
        manager/cuda_device_manager.cu
        manager/cuda_memory_manager.cu
        algorithm/constraints.cu
    )

    add_library(grocuda STATIC ${GPU_SOURCES} ${GPU_CUDA_SOURCES})

    set_target_properties(grocuda PROPERTIES
        CUDA_RUNTIME_LIBRARY Shared
        # CUDA_RUNTIME_LIBRARY Static ## use this for production
        CUDA_SEPARABLE_COMPILATION ON
        CUDA_ARCHITECTURES "60;62;70;75;80;86"
    )

    target_link_libraries(grocuda PRIVATE CUDA::cudart ${CUDAToolkit_LIBRARIES})
    target_compile_definitions(grocuda PRIVATE USE_CUDA)
else()
    # No CUDA, fallback to CPU implementation of constraints
    set(NO_GPU_SOURCES
        algorithm/constraints.cc
        manager/cuda_device_worker.cc
        manager/cuda_device_manager.cc
        manager/cuda_memory_manager.cc
    )

    add_library(grocuda STATIC ${GPU_SOURCES} ${NO_GPU_SOURCES})

endif()

set_target_properties(grocuda PROPERTIES
    CXX_STANDARD 17
    POSITION_INDEPENDENT_CODE ON
)

# Include paths (always the same)
target_include_directories(grocuda PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/algorithm
    ${CMAKE_CURRENT_SOURCE_DIR}/manager
    ${CMAKE_CURRENT_SOURCE_DIR}/kernels
    ${CMAKE_CURRENT_SOURCE_DIR}/memory
    ${PROJECT_BINARY_DIR}
    ${EXTERNAL_INCLUDES}
    ${PROJECT_SOURCE_DIR}/src
)
