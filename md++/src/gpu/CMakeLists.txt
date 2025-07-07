# CUDA source files
set(CUDA_SOURCES
    kernels/md_kernels.cu
)

# Header files
set(CUDA_HEADERS
    cuda_manager/cuda_manager.h
    cuda_manager/cuda_manager.tpp
    cuda_manager/cuda_device_manager.h
    cuda_manager/cuda_device_worker.h
    cuda_manager/cuda_memory_manager.h
    memory/cuvector.h
    memory/cuallocator.h
    kernels/md_kernels.h
)

# Add a STATIC library for grocuda
add_library(grocuda STATIC ${CUDA_SOURCES} ${CUDA_HEADERS})