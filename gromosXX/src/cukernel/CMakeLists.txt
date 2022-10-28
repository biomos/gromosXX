set(CUDA_CC
    cudaKernel.cu  
    lib/utils.cu 
)

add_library(libgrocukernel OBJECT ${CUDA_CC})

target_include_directories(libgrocukernel PUBLIC
    ${PROJECT_BINARY_DIR}
    .
    ..
    lib
)
