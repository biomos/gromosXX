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
 * @file identity_token.h
 * Process-wide monotonic identity tokens (PLAN.md §3.2's "cache key
 * lifetime" resolution): used by topology::Topology and configuration::
 * Configuration to give every *object* (not just every distinct address)
 * a unique id, so CudaManager's identity-keyed GPU-mirror cache can
 * detect "this address was reused by an unrelated object" as a cache
 * miss instead of silently handing back a stale mirror for the wrong
 * object -- a correctness bug (silent wrong data), not just a crash.
 *
 * Deliberately has no GPU/CUDA dependency -- topology.h/configuration.h
 * need this in both USE_CUDA and non-CUDA builds (an object still needs
 * an identity even when nothing ever looks it up in a GPU cache), and
 * PLAN.md §3.2 explicitly calls out removing topology/configuration's
 * gpu/cuda/... header coupling as one of the points of this refactor --
 * adding a new one here would defeat that.
 */

#pragma once

#include <atomic>
#include <cstddef>

namespace util {

  /**
   * Returns a fresh id on every call, starting at 1 (0 is reserved as
   * "no id assigned" / "not a real object", e.g. for a default-
   * constructed cache-miss sentinel if one is ever needed) and never
   * repeating for the lifetime of the process. Thread-safe (the counter
   * itself, not any cache keyed on the result -- CudaManager's cache is
   * explicitly single-threaded, see PLAN.md §3.2's thread-safety note;
   * this only needs to be safe because object construction can
   * legitimately happen from more than one thread, e.g. MPI/replica
   * setups constructing independent Topology/Configuration objects).
   */
  inline std::size_t next_identity_token() {
    static std::atomic<std::size_t> counter{1};
    return counter.fetch_add(1, std::memory_order_relaxed);
  }

} // namespace util
