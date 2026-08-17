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
 * @file feature_checker_gpu.t.cc
 * Regression test for io::check_features()'s "gpu" feature-
 * compatibility matrix (io/parameter/check_parameter.cc) -- that
 * matrix predates the real CUDA implementation (PLAN.md §10 steps
 * 8-20) and, before this test's fix, still hard-errored on a normal
 * accelerator=cuda run using bonded forces (angle/dihedral/improper),
 * the standard LJ+CRF nonbonded interaction, COM motion removal, a
 * grid pairlist selector, chargegroup-based cutoff, and multiple
 * energy groups -- every one of which has a real, tested GPU
 * implementation elsewhere in this suite (quartic_bond_gpu.t.cc,
 * angle_gpu.t.cc, dihedral_gpu.t.cc, improper_dihedral_gpu.t.cc,
 * cuda_nonbonded_interaction.t.cc, remove_com_motion_gpu.t.cc,
 * cuda_atomic_cutoff_toggle.t.cc). This test drives
 * io::check_features() directly (no topology/trajectory needed, it's
 * a pure parameter cross-check) with exactly that combination and
 * checks it no longer fails. Runs on every build (this file has
 * nothing CUDA-specific in it -- it's checking parameter-level
 * bookkeeping, not running a kernel), so it also guards against a
 * *new* accelerator=cuda combination regressing this matrix again in
 * the future, CPU-only builds included.
 */

#include "../stdheader.h"

#include "../topology/topology.h"
#include "../simulation/simulation.h"
#include "../io/message.h"
#include "../io/parameter/check_parameter.h"

#include "check.h"

#ifdef XXMPI
  #include <mpi.h>
#endif

namespace {

  bool flush_messages(const char * label) {
    const io::message::severity_enum sev = io::messages.display(std::cerr);
    io::messages.clear();
    return sev >= io::message::error;
  }

  int run_case() {
    simulation::Simulation sim;
    simulation::Parameter & param = sim.param();

    param.gpu.accelerator = simulation::gpu_cuda;

    // Rectangular PBC (matching a normal periodic run, e.g. aladip) --
    // param.boundary.boundary defaults to vacuum, which trips the
    // unrelated, pre-existing "vacuum simulation does not work with
    // grid-based pairlist algorithm" lock (nothing to do with GPU).
    param.boundary.boundary = math::rectangular;

    // rf_excluded defaults to true ("new standard", parameter.h) but is
    // NOT implemented on GPU (see check_parameter.cc's comment) --
    // explicitly off here, since this test is about the locks that
    // *should* be lifted, not about that one, which correctly stays
    // locked.
    param.nonbonded.rf_excluded = false;

    // Bonded terms actually GPU-ported this session.
    param.force.bond = 1;
    param.force.angle = 1;
    param.force.dihedral = 1;
    param.force.improper = 1;

    // Standard (non-shifted) LJ + reaction-field nonbonded, the only
    // formula CUDA_Nonbonded_Interaction implements.
    param.force.nonbonded_crf = 1;
    param.force.nonbonded_vdw = 1;

    // COM motion removal (Remove_COM_Motion<gpuBackend>).
    param.centreofmass.remove_rot = true;
    param.centreofmass.remove_trans = true;

    // The pairlist.grid selector and pairlist.atomic_cutoff toggle are
    // both irrelevant/covered under accelerator=cuda (see the comments
    // in check_parameter.cc) -- exercise the "grid" + chargegroup-cutoff
    // combination specifically, since that's what actually triggered
    // this bug report.
    param.pairlist.grid = 1;
    param.pairlist.atomic_cutoff = false;

    // Multiple energy groups (bucketed by every GPU kernel this session
    // touched).
    param.force.energy_group.assign({11u, 71u});

    const int rc = io::check_features(sim);
    const bool had_errors = flush_messages("check_features");

    if (rc != 0 || had_errors) {
      std::cerr << "feature_checker_gpu: FAILED -- io::check_features() "
                << "rejected a GPU-accelerated bonded+nonbonded+COM-"
                << "removal+multi-energy-group run "
                << "(rc=" << rc << ")" << std::endl;
      return 1;
    }

    std::cout << "feature_checker_gpu: OK" << std::endl;
    return 0;
  }

} // namespace

int main(int argc, char* argv[]) {
#ifdef XXMPI
  MPI_Init(&argc, &argv);
#endif
  return run_case();
}
