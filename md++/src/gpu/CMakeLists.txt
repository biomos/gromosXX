# Add GPU header files (not compiled directly)
set(GPU_HEADERS
    gpu.h
    cuda/cuheader.h
    cuda/manager/cuda_memory_manager.h
    cuda/manager/cuda_device_worker.h
    cuda/manager/cuda_manager.h
    cuda/manager/cuda_device_manager.h
    cuda/utils.h
    cuda/algorithm/constraints.h
    cuda/kernels/reduction.h
    cuda/kernels/indexing.h
    cuda/kernels/hello_world.h
    cuda/memory/float3.h
    cuda/memory/float9.h
    cuda/memory/types.h
    cuda/memory/float2.h
    cuda/memory/constants.h
    cuda/memory/float4.h
    cuda/memory/cupointer.h
    cuda/memory/cuvector.h
    cuda/memory/precision.h
    cuda/memory/array.h
    cuda/memory/cumallocator.h
    cuda/memory/pairlist/tile.h
    cuda/memory/pairlist/pairlist.h
    cuda/memory/configuration_struct.h
    cuda/memory/cuhallocator.h
    cuda/memory/topology_struct.h
    cuda/memory/container.h
    cuda/interaction/nonbonded/pairlist/cuda_pairlist_algorithm_impl_gpu.h
)

if(USE_CUDA)

    set(GPU_SOURCES
        gpu.cc
        cuda/utils.cu
        cuda/manager/cuda_manager.cu
        cuda/kernels/reduction.cu
        cuda/kernels/hello_world.cu
        cuda/memory/container.cu
        cuda/memory/topology_struct.cu
        cuda/memory/configuration_struct.cu
        cuda/manager/cuda_device_worker.cu
        cuda/manager/cuda_device_manager.cu
        cuda/manager/cuda_memory_manager.cu
        # cuda/algorithm/constraints.cu
        cuda/interaction/nonbonded/pairlist/cuda_pairlist_algorithm_impl_gpu.cu
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
        gpu.cc
        cuda/manager/cuda_manager.cc
        # cuda/algorithm/constraints.cc
    )

    add_library(grocuda STATIC ${GPU_SOURCES})

    set_target_properties(grocuda PROPERTIES
        POSITION_INDEPENDENT_CODE ON
        CXX_STANDARD 17
    )
endif()

target_include_directories(grocuda PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/algorithm
    ${CMAKE_CURRENT_SOURCE_DIR}/kernels
    ${CMAKE_CURRENT_SOURCE_DIR}/manager
    ${CMAKE_CURRENT_SOURCE_DIR}/memory
    ${PROJECT_BINARY_DIR}
    ${EXTERNAL_INCLUDES}
    ${PROJECT_SOURCE_DIR}/src
    ${PROJECT_SOURCE_DIR}/src/gpu
    ${PROJECT_SOURCE_DIR}/src/gpu/cuda
)
