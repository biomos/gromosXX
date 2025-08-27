/**
 * @file hello_world.h
 * @author poliak
 * example kernel function
 */

#pragma once


namespace gpu {
    __global__ void put_chargegroups_into_box(gpu::Topology topo,
                                            gpu::Configuration::View conf);
}