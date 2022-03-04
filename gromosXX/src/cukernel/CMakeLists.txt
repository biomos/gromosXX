set(CUDA_CC
    cudaKernel.cu  
    lib/utils.cu 
)

add_library(cukernel OBJECT ${CUDA_CC})

target_include_directories(cukernel PUBLIC
    ${PROJECT_BINARY_DIR}
    .
    ..
    lib
)
