/**
 * @file hello_world.h
 * @author poliak
 * example kernel function
 */

#pragma once


namespace gpu {
    // struct topology_struct;
    // struct configuration_state_struct;
    // __global__ void hello_world(gpu::topology_struct& topo, gpu::configuration_state_struct& conf);
    __global__ void hello_world(float* a, float* b);
}