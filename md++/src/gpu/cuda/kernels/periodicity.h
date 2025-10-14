/**
 * @file periodicity.h
 * @author poliak
 * Kernel functions for periodic boundary conditions implementation
 */

#pragma once


namespace gpu {
    template <math::boundary_enum BOUNDARY>
    __global__ void put_chargegroups_into_box_kernel(Topology::View topo,
                                                    Configuration::View conf,
                                                    Periodicity<BOUNDARY> periodicity);

    template <math::boundary_enum BOUNDARY>
    __global__ void prepare_cog_kernel(Topology::View topo,
                                       Configuration::View conf,
                                       Periodicity<BOUNDARY> periodicity,
                                       math::CuVArray::View cg_cog,
                                       gpu::cuvector<ushort4>::View cg_cells
                                    );
}
