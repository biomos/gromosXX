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
 * @file constraint_error_slots.h
 * Named slots for the deferred constraint-error-flags buffer
 * (CudaManager). Plain integer constants + a pure translation
 * function, zero CUDA dependency -- includable from CPU-only builds
 * (Algorithm_Sequence::run() needs this even when USE_CUDA is off).
 *
 * Design: every GPU constraint algorithm writes its status into its
 * own reserved slot of one small, persistent, GPU-resident int array
 * instead of checking (and cudaDeviceSynchronize()-ing on) its own
 * private flag immediately after its own kernels. The whole array is
 * zeroed once per step and checked once, at the very end of
 * Algorithm_Sequence::run() -- one sync instead of one per algorithm.
 * Safe because a constraint failure is already a fatal, whole-step-
 * discarding event in GROMOS (program/md.cc aborts the run on any
 * nonzero Algorithm_Sequence::run() return) -- running the rest of
 * that step's kernels on top of already-corrupted positions before
 * the failure is noticed costs nothing extra, since the entire step
 * is being thrown away either way.
 *
 * Two different semantics share this one mechanism:
 *  - "Fatal" slots (SHAKE/M-SHAKE/SETTLE): 0 = ok, nonzero = a real
 *    error code (matching each algorithm's existing convention: 1 =
 *    orthogonal/singular, 2 = too many iterations) that must abort the
 *    run, same as today's immediate check did.
 *  - "Informational counter" slots (LINCS's rotation-count
 *    diagnostic): a plain count, atomicAdd'd, purely advisory --
 *    matches the CPU/pre-existing GPU behaviour of printing a notice,
 *    never a fatal return.
 */

#pragma once

#include <string>

namespace gpu {

  enum ConstraintErrorSlot : unsigned {
    ERR_SLOT_SHAKE_SOLUTE   = 0,
    ERR_SLOT_SHAKE_SOLVENT  = 1,
    ERR_SLOT_M_SHAKE        = 2,
    ERR_SLOT_LINCS_SOLUTE   = 3,
    ERR_SLOT_LINCS_SOLVENT  = 4,
    ERR_SLOT_SETTLE         = 5,
    NUM_CONSTRAINT_ERROR_SLOTS = 6,
  };

  /**
   * Whether a nonzero value in this slot is fatal (should abort the
   * run) or purely informational (LINCS's rotation counter).
   */
  inline bool constraint_error_slot_is_fatal(unsigned slot) {
    return slot != ERR_SLOT_LINCS_SOLUTE && slot != ERR_SLOT_LINCS_SOLVENT;
  }

  /**
   * Human-readable message for a nonzero slot value, matching each
   * algorithm's pre-existing immediate-check message text as closely
   * as possible so deferred detection doesn't change what the user
   * sees when a run fails.
   */
  inline std::string describe_constraint_error(unsigned slot, int code) {
    switch (slot) {
      case ERR_SLOT_SHAKE_SOLUTE:
        return code == 1 ? "SHAKE error. vectors orthogonal (solute)"
                          : "SHAKE error. too many iterations (solute)";
      case ERR_SLOT_SHAKE_SOLVENT:
        return code == 1 ? "SHAKE error. vectors orthogonal (solvent)"
                          : "SHAKE error. too many iterations (solvent)";
      case ERR_SLOT_M_SHAKE:
        return code == 1 ? "M_SHAKE error. vectors orthogonal (solvent)"
                          : "M_SHAKE error. too many iterations (solvent)";
      case ERR_SLOT_LINCS_SOLUTE:
        return "LINCS: too much rotation in " + std::to_string(code) + " case(s) (solute)";
      case ERR_SLOT_LINCS_SOLVENT:
        return "LINCS: too much rotation in " + std::to_string(code) + " case(s) (solvent)";
      case ERR_SLOT_SETTLE:
        return "SETTLE error";
      default:
        return "unknown constraint error slot " + std::to_string(slot);
    }
  }

} // namespace gpu
