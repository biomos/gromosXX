set(CUDA_CC
    constraints.cu
    cudaKernel.cu
    cudaShake.cu  
    interaction.cu
    io.cu
    pairlist.cu
    lib/utils.cu 
)

add_library(grocukernel OBJECT ${CUDA_CC})

target_include_directories(grocukernel PUBLIC
    ${PROJECT_BINARY_DIR}
    ${EXTERNAL_INCLUDES}
    .
    ..
    lib
)
