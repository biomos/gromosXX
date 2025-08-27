/**
 * @file hello_world.h
 * @author poliak
 * example kernel function
 */

#pragma once


namespace gpu {
    __global__ void hello_world(gpu::Topology topo,
                                gpu::Configuration::View conf);
    // __global__ void hello_world(float* a, float* b);
}