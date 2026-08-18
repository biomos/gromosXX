# Add GPU header files (not compiled directly)
set(GPU_HEADERS
    gpu.h
    cuda/cuheader.h
    cuda/cuhostdevice.h
    cuda/manager/cuda_memory_manager.h
    cuda/manager/cuda_device_worker.h
    cuda/manager/cuda_manager.h
    cuda/manager/cuda_device_manager.h
    cuda/utils.h
    cuda/kernels/reduction.h
    cuda/kernels/indexing.h
    cuda/kernels/hello_world.h
    cuda/kernels/periodicity.h
    cuda/math/box.h
    cuda/math/periodicity.h
    cuda/memory/float3.h
    cuda/memory/float9.h
    cuda/memory/types.h
    cuda/memory/float2.h
    cuda/memory/constants.h
    cuda/memory/float4.h
    cuda/memory/cupointer.h
    cuda/memory/cuvector.h
    cuda/memory/precision.h
    cuda/memory/cumallocator.h
    cuda/memory/pairlist/tile.h
    # cuda/memory/pairlist/pairlist.h
    cuda/memory/configuration_struct.h
    cuda/memory/cuhallocator.h
    cuda/memory/topology_struct.h
    cuda/memory/container.h
    cuda/interaction/nonbonded/pairlist/cuda_pairlist_algorithm_impl.h
    cuda/interaction/nonbonded/pairlist/block_pairlist.h
    cuda/interaction/nonbonded/cuda_lj_params.h
    cuda/interaction/nonbonded/cuda_nb_sim_params.h
    cuda/interaction/nonbonded/kernels/lj_crf_tiles.h
    cuda/interaction/nonbonded/kernels/rf_excluded_kernels.h
    cuda/interaction/nonbonded/kernels/one_four_kernels.h
    cuda/interaction/nonbonded/kernels/displacement.h
    cuda/algorithm/integration/leap_frog_kernels.h
    cuda/algorithm/constraints/remove_com_motion_kernels.h
    cuda/algorithm/temperature/temperature_kernels.h
    cuda/interaction/bonded/quartic_bond_kernels.h
    cuda/interaction/bonded/angle_kernels.h
    cuda/interaction/bonded/improper_dihedral_kernels.h
    cuda/interaction/bonded/dihedral_kernels.h
    cuda/algorithm/constraints/shake_kernels.h
    cuda/algorithm/constraints/settle_kernels.h
    cuda/algorithm/constraints/lincs_kernels.h
    cuda/memory/vec3_convert.h
    cuda/algorithm/temperature/atom_bath_arrays.h
    cuda/algorithm/temperature/thermostat_velocity_scale.h
    cuda/interaction/bonded/sparse_force_accumulate.h
)

if(USE_CUDA)

    set(GPU_SOURCES
        gpu.cc
        cuda/utils.cu
        cuda/manager/cuda_manager.cu
        cuda/kernels/reduction.cu
        cuda/kernels/hello_world.cu
        cuda/kernels/periodicity.cu
        cuda/math/periodicity.cu
        cuda/memory/container.cu
        cuda/memory/topology_struct.cu
        cuda/memory/configuration_struct.cu
        cuda/manager/cuda_device_worker.cu
        cuda/manager/cuda_device_manager.cu
        cuda/manager/cuda_memory_manager.cu
        cuda/interaction/nonbonded/pairlist/cuda_pairlist_algorithm_impl.cu
        cuda/interaction/nonbonded/pairlist/block_pairlist.cu
        cuda/memory/pairlist/tile.cu
        cuda/interaction/nonbonded/cuda_lj_params.cu
        cuda/interaction/nonbonded/kernels/lj_crf_tiles.cu
        cuda/interaction/nonbonded/kernels/rf_excluded_kernels.cu
        cuda/interaction/nonbonded/kernels/one_four_kernels.cu
        cuda/interaction/nonbonded/kernels/displacement.cu
        cuda/algorithm/integration/leap_frog_kernels.cu
        cuda/algorithm/constraints/remove_com_motion_kernels.cu
        cuda/algorithm/temperature/temperature_kernels.cu
        cuda/interaction/bonded/quartic_bond_kernels.cu
        cuda/interaction/bonded/angle_kernels.cu
        cuda/interaction/bonded/improper_dihedral_kernels.cu
        cuda/interaction/bonded/dihedral_kernels.cu
        cuda/algorithm/constraints/shake_kernels.cu
        cuda/algorithm/constraints/settle_kernels.cu
        cuda/algorithm/constraints/lincs_kernels.cu
    )

    add_library(grocuda STATIC ${GPU_SOURCES})

    # Needed since block_pairlist.cu calls TileVecT<Interaction_Tile>::
    # push_back (a __device__ function defined in tile.cu) -- resolving a
    # device-side function call across .cu translation units requires
    # relocatable device code, unlike a host-launched __global__ kernel
    # defined in another .cu file (which already worked without this,
    # e.g. prepare_cog_kernel).
    set_target_properties(grocuda PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON
    )

    target_link_libraries(grocuda PUBLIC CUDA::cudart ${CUDAToolkit_LIBRARIES})
    target_compile_definitions(grocuda PRIVATE USE_CUDA)

else()

    set(GPU_SOURCES
        gpu.cc
        cuda/manager/cuda_manager.cc
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
    ${CMAKE_CURRENT_SOURCE_DIR}/math
    ${CMAKE_CURRENT_SOURCE_DIR}/memory
    ${PROJECT_BINARY_DIR}
    ${EXTERNAL_INCLUDES}
    ${PROJECT_SOURCE_DIR}/src
    ${PROJECT_SOURCE_DIR}/src/gpu
    ${PROJECT_SOURCE_DIR}/src/gpu/cuda
)
