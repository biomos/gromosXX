/**
 * @file periodicity.h
 * @author poliak
 * Kernel functions for periodic boundary conditions implementation
 */

#pragma once


namespace gpu {
    __global__ void put_chargegroups_into_box(Topology topo,
                                              Configuration::View conf);

    template <math::boundary_enum BOUNDARY>
    __global__ void prepare_cog_kernel(Topology topo,
                                       Configuration::View conf,
                                       Periodicity<BOUNDARY> periodicity);
}
