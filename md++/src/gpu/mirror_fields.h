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
    MIRROR_POS   = 1u << 0,
    MIRROR_VEL   = 1u << 1,
    MIRROR_FORCE = 1u << 2,
    MIRROR_BOX   = 1u << 3,
    MIRROR_ALL   = MIRROR_POS | MIRROR_VEL | MIRROR_FORCE | MIRROR_BOX,
  };

}
