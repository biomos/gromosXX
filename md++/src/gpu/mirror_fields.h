/*
 * This file is part of GROMOS.
 *
 * Copyright (c) 2011, 2012, 2016, 2018, 2021, 2023 Biomos b.v.
 * See <https://www.gromos.net> for details.
 *
 * GROMOS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

/**
 * @file mirror_fields.h
 * Field-granularity flags for the CudaManager Configuration-mirror
 * freshness tracker (data-level cache coherence between
 * configuration::Configuration and its GPU mirror, gpu::Configuration).
 * Plain integer constants, zero CUDA dependency -- must be includable
 * from both CPU-only and CUDA builds (algorithm::Algorithm's
 * gpu_mirror_touches() hook needs this even when USE_CUDA is off).
 */

#pragma once

namespace gpu {

  /**
   * @brief Bits identify which part of a Configuration mirror
   * (current()+old() together, matching how copy_to_device()/
   * copy_pos_vel_to_device() already move both in one call) a piece of
   * code just made trustworthy on the GPU side, or might have
   * invalidated by writing the CPU-side configuration::Configuration
   * directly.
   */
  enum MirrorField : unsigned {
    MIRROR_POS             = 1u << 0,
    MIRROR_VEL             = 1u << 1,
    MIRROR_FORCE           = 1u << 2,
    MIRROR_BOX             = 1u << 3,
    /**
     * conf.{current,old}().constraint_force -- written by the
     * constraint algorithms (SHAKE/M-SHAKE/LINCS/SETTLE), each
     * accumulating into the shared GPU-resident buffer instead of a
     * private per-algorithm one, so consecutive constraint algorithms
     * in the same step never round-trip through the CPU between them.
     */
    MIRROR_CONSTRAINT_FORCE = 1u << 4,
    /**
     * conf.{current,old}().virial_tensor -- same rationale as
     * MIRROR_CONSTRAINT_FORCE. Deliberately not bundled with
     * MIRROR_FORCE: a plain force-field algorithm touches FORCE but
     * not necessarily VIRIAL (and vice versa for e.g.
     * Molecular_Virial_Interaction), so keeping them independent bits
     * avoids one causing an unnecessary resync of the other.
     */
    MIRROR_VIRIAL           = 1u << 5,
    /**
     * conf.special().lattice_shifts -- unlike every other bit here,
     * this doesn't cycle with current()/old() (it's a single
     * persistent per-atom array, not part of either state) -- see
     * gpu::Configuration's own `lattice_shifts` member (a sibling of
     * `current`/`old`, not inside either). Written by
     * Lattice_Shift_Tracker<gpuBackend> only, and read by nothing else
     * anywhere in the codebase except the final trajectory write
     * (io::Out_Configuration::_print_lattice_shifts(), form==final
     * only). Deliberately NOT part of MIRROR_ALL: unlike CONSTRAINT_
     * FORCE/VIRIAL (which other algorithms genuinely read), including
     * it there would mean every ordinary default-MIRROR_ALL algorithm
     * invalidates it every step -- since freshness invalidation forces
     * a real O(num_atoms) re-upload the next time it's requested
     * (unlike a true no-op), that turned Lattice_Shift_Tracker's own
     * apply() into a real-transfer-every-step cost again (measured:
     * 0.8s -> 17.5s over 10000 steps), defeating the entire point of
     * this port. Left fully self-managed instead, same as POS in
     * Lattice_Shift_Tracker<gpuBackend>'s own gpu_mirror_touches()==0.
     */
    MIRROR_LATTICE_SHIFT    = 1u << 6,
    MIRROR_ALL   = MIRROR_POS | MIRROR_VEL | MIRROR_FORCE | MIRROR_BOX |
                   MIRROR_CONSTRAINT_FORCE | MIRROR_VIRIAL,
  };

}
