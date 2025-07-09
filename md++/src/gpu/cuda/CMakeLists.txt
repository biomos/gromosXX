# Add GPU header files (not compiled directly)
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

if(USE_CUDA)

    set(GPU_SOURCES
        manager/cuda_manager.cc
        kernels/reduction.cu
        memory/container.cu
        manager/cuda_device_worker.cu
        manager/cuda_device_manager.cu
        manager/cuda_memory_manager.cu
        algorithm/constraints.cu
    )

    add_library(grocuda STATIC ${GPU_SOURCES})

    set_target_properties(grocuda PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON
        POSITION_INDEPENDENT_CODE ON
        CXX_STANDARD 17
        CUDA_ARCHITECTURES "60;62;70;75;80;86"
        CUDA_RUNTIME_LIBRARY Static  # Optional
    )

    target_link_libraries(grocuda PUBLIC CUDA::cudart ${CUDAToolkit_LIBRARIES})
    target_compile_definitions(grocuda PRIVATE USE_CUDA)

else()

    set(GPU_SOURCES
        algorithm/constraints.cc
        manager/cuda_manager.cc
        manager/cuda_device_worker.cc
        manager/cuda_device_manager.cc
        manager/cuda_memory_manager.cc
    )

    add_library(grocuda STATIC ${GPU_SOURCES})

    set_target_properties(grocuda PROPERTIES
        POSITION_INDEPENDENT_CODE ON
        CXX_STANDARD 17
    )
endif()

target_include_directories(grocuda PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/algorithm
    ${CMAKE_CURRENT_SOURCE_DIR}/manager
    ${CMAKE_CURRENT_SOURCE_DIR}/kernels
    ${CMAKE_CURRENT_SOURCE_DIR}/memory
    ${PROJECT_BINARY_DIR}
    ${EXTERNAL_INCLUDES}
    ${PROJECT_SOURCE_DIR}/src
)
