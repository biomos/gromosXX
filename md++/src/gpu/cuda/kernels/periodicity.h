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

    /**
     * @param sort_key flat output of cg_cells[cg_i].w (the Morton cell
     * index), one entry per chargegroup -- a plain unsigned array is
     * easier to feed straight into a Thrust sort-by-key than pulling .w
     * out of the ushort4 cg_cells array (TILE_PAIRLIST_DESIGN.md §3 step 2).
     */
    template <math::boundary_enum BOUNDARY>
    __global__ void prepare_cog_kernel(Topology::View topo,
                                       Configuration::View conf,
                                       Periodicity<BOUNDARY> periodicity,
                                       math::CuVArray::View cg_cog,
                                       gpu::cuvector<ushort4>::View cg_cells,
                                       gpu::cuvector<unsigned>::View sort_key
                                    );
}
