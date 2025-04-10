set(CUDA_CC
  io.cu
  pairlist.cu
  interaction.cu
  constraints.cu
  cudaShake.cu
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
target_link_libraries(grocukernel PRIVATE ${EXTERNAL_LIBRARIES})
