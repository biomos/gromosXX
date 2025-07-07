# Add GPU source files
set(GPU_SOURCES
    manager/cuda_device_worker.cc
    manager/cuda_device_manager.cc
)

# Add GPU header files (for clarity, not compiled directly)
set(GPU_HEADERS
    utils.h
    manager/cuda_memory_manager.h
    manager/cuda_device_worker.h
    manager/cuda_device_worker.tcc
    manager/cuda_manager.h
    manager/cuda_device_manager.h
    manager/cuda_manager.tcc
    manager/cuda_memory_manager.tcc
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
    CUDA_ARCHITECTURES "60;70;75;80;86;90" # Explicitly set CUDA architectures for this target
)

# Include directories for GPU code
target_include_directories(grocuda PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/manager
    ${CMAKE_CURRENT_SOURCE_DIR}/kernels
    ${CMAKE_CURRENT_SOURCE_DIR}/memory
    ${PROJECT_BINARY_DIR}
    ${EXTERNAL_INCLUDES}
    /usr/local/cuda/include # Explicitly add the CUDA include directory
)

# Link CUDA libraries
target_link_libraries(grocuda PRIVATE cuda cudart)