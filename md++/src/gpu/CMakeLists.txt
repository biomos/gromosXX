# Add GPU source files
set(GPU_SOURCES
    manager/cuda_device_worker.cc
    manager/cuda_device_manager.cc
)

# Add GPU header files (for clarity, not compiled directly)
set(GPU_HEADERS
    manager/cuda_memory_manager.h
    manager/cuda_device_worker.h
    manager/cuda_device_worker.tcc
    manager/cuda_manager.h
    manager/cuda_device_manager.h
    manager/cuda_manager.tcc
    manager/cuda_memory_manager.tcc
    utils.h
    kernels/reduction.h
    kernels/indexing.h
    memory/constants.h
    memory/cuvector.h
    memory/cumallocator.h
    memory/container.h
)

# Add GPU CUDA source files
set(GPU_CUDA_SOURCES
    kernels/reduction.cu
    memory/container.cu
)

# Create an object library for GPU code
add_library(grocuda OBJECT ${GPU_SOURCES} ${GPU_CUDA_SOURCES} ${GPU_HEADERS})

# Set CUDA-specific properties
set_target_properties(grocuda PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    POSITION_INDEPENDENT_CODE ON
)

# Include directories for GPU code
target_include_directories(grocuda PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/manager
    ${CMAKE_CURRENT_SOURCE_DIR}/utils
    ${CMAKE_CURRENT_SOURCE_DIR}/kernels
    ${CMAKE_CURRENT_SOURCE_DIR}/memory
    ${PROJECT_BINARY_DIR}
)

# Link CUDA libraries
target_link_libraries(grocuda PRIVATE cuda cudart)