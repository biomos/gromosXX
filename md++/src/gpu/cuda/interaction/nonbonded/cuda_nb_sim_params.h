/**
 * @file cuda_nb_sim_params.h
 * Plain-old-data struct passed by value to GPU force kernels,
 * holding all simulation constants needed for LJ-CRF nonbonded forces.
 */

#pragma once

#include "gpu/cuda/cuhostdevice.h"
#include "gpu/cuda/memory/precision.h"

namespace gpu {

/**
 * @brief Simulation constants for nonbonded LJ-CRF interaction.
 *
 * All values are in GROMOS internal units.
 * Passed by value into every GPU force kernel (fits in registers).
 */
struct NbSimParams {
    /** 1/(4*pi*eps0*eps_r) in kJ mol^-1 nm e^-2 */
    FPL_TYPE four_pi_eps_i;
    /** Reaction-field constant krf = crf / (2 * rc^3)  [nm^-3] */
    FPL_TYPE crf_2cut3i;
    /** Reaction-field energy constant Crf = (1 - crf/2) / rc  [nm^-1] */
    FPL_TYPE crf_cut;
    /** Short-range cutoff squared [nm^2] */
    FPL_TYPE cutoff_short_sq;
    /** Long-range cutoff squared [nm^2] */
    FPL_TYPE cutoff_long_sq;
    /** Number of energy groups (topo.energy_groups().size()) -- the tile
     * kernel's energy reduction buckets by [eg_i * num_energy_groups +
     * eg_j], matching configuration::Energy::lj_energy/crf_energy's
     * [num_groups][num_groups] shape. */
    unsigned num_energy_groups;
    /** 1,4-pair Coulomb scaling factor (Nonbonded_Parameter::
     * get_coulomb_scaling(), gpu::launch_one_four) -- 1.0 unless
     * param.amber.amber, matching create_nonbonded.cc's
     * set_coulomb_scaling call. */
    FPL_TYPE coulomb_scaling;
};

} // namespace gpu
