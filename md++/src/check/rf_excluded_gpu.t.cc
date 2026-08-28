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
 * @file rf_excluded_gpu.t.cc
 * Correctness test for `gpu::launch_rf_excluded` (rf_excluded_kernels.cu),
 * wired into CUDA_Pairlist_Algorithm_Impl::compute_forces_energies()
 * whenever `param.nonbonded.rf_excluded` is set. Structured the same way
 * as cuda_nonbonded_interaction.t.cc (real CUDA_Pairlist_Algorithm +
 * CUDA_Nonbonded_Interaction pair, compared against a direct CPU
 * reference), but the CPU reference here is built from
 * Nonbonded_Term::rf_interaction over topo.exclusion(i) (solute
 * self-term + excluded pairs) plus the rigid-solvent excluded-pair energy
 * formula, exactly mirroring RF_excluded_outerloop's two innerloops
 * (nonbonded_innerloop.cc) rather than the regular pairlist-driven
 * LJ/CRF path cuda_nonbonded_interaction.t.cc already covers.
 *
 * aladip's own solute has real bonded exclusions (1-2/1-3/1-4) and its
 * solvent is SPC water (3-atom chargegroups) -- both branches of
 * gpu::launch_rf_excluded are genuinely exercised without needing a
 * synthetic topology.
 */

#include "../stdheader.h"

#include "../algorithm/algorithm.h"
#include "../topology/topology.h"
#include "../simulation/simulation.h"
#include "../configuration/configuration.h"

#include "../algorithm/algorithm/algorithm_sequence.h"

#include "../io/argument.h"
#include "../util/parse_verbosity.h"
#include "../util/usage.h"
#include "../io/topology/in_topology.h"
#include "../io/message.h"

#include "../util/create_simulation.h"

#include "../math/periodicity.h"
#include "../interaction/nonbonded/interaction/nonbonded_parameter.h"
#include "../interaction/nonbonded/interaction/nonbonded_term.h"
#include "../interaction/nonbonded/pairlist/pairlist.h"
#include "../interaction/nonbonded/pairlist/pairlist_algorithm.h"
#include "../interaction/nonbonded/pairlist/cuda_pairlist_algorithm.h"
#include "../interaction/nonbonded/interaction/cuda_nonbonded_interaction.h"

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

  struct Reference {
    std::vector<math::Vec> force;
    std::vector<std::vector<double> > crf_energy;
  };

  /**
   * Direct CPU reference for the RF-excluded contribution only, built
   * exactly like RF_excluded_outerloop's two innerloops
   * (nonbonded_innerloop.cc): solute self-term + topo.exclusion(i) pairs
   * via Nonbonded_Term::rf_interaction, and solvent chargegroups' rigid
   * internal-pair energy (no force, no self-term -- see that innerloop's
   * comment on why the distance-independent part is left out).
   */
  template <math::boundary_enum B>
  Reference compute_reference(topology::Topology & topo,
                               configuration::Configuration & conf,
                               simulation::Simulation & sim) {
    interaction::Nonbonded_Term term;
    term.init(sim);
    math::Periodicity<B> periodicity(conf.current().box);

    const unsigned num_groups = static_cast<unsigned>(topo.energy_groups().size());
    const unsigned num_solute_atoms = topo.num_solute_atoms();

    Reference ref;
    ref.force.resize(topo.num_atoms(), math::Vec(0.0, 0.0, 0.0));
    ref.crf_energy.assign(num_groups, std::vector<double>(num_groups, 0.0));

    // Solute: self-term + excluded pairs.
    for (unsigned i = 0; i < num_solute_atoms; ++i) {
      const double qi = topo.charge(i);
      math::Vec r(0.0, 0.0, 0.0), f;
      double e_crf = 0.0;
      term.rf_interaction(r, qi * qi, f, e_crf);
      const unsigned eg_i = topo.atom_energy_group(i);
      ref.crf_energy[eg_i][eg_i] += 0.5 * e_crf;

      for (topology::excl_cont_t::value_type::const_iterator
               it = topo.exclusion(i).begin(), to = topo.exclusion(i).end();
           it != to; ++it) {
        const unsigned j = *it;
        periodicity.nearest_image(conf.current().pos(i), conf.current().pos(j), r);
        term.rf_interaction(r, qi * topo.charge(j), f, e_crf);
        ref.force[i] += f;
        ref.force[j] -= f;
        const unsigned eg_j = topo.atom_energy_group(j);
        ref.crf_energy[eg_i][eg_j] += e_crf;
      }
    }

    // Solvent: rigid intramolecular pairs, energy only, no self-term.
    for (topology::Chargegroup_Iterator cg_it = topo.chargegroup_it(topo.num_solute_chargegroups()),
                                         cg_to = topo.chargegroup_end();
         cg_it != cg_to; ++cg_it) {
      for (topology::Atom_Iterator at_it = cg_it.begin(), at_end = cg_it.end();
           at_it != at_end; ++at_it) {
        for (topology::Atom_Iterator at2_it = at_it + 1; at2_it != at_end; ++at2_it) {
          math::Vec r;
          periodicity.nearest_image(conf.current().pos(*at_it), conf.current().pos(*at2_it), r);
          const double e_crf = -(topo.charge(*at_it) * topo.charge(*at2_it)) *
                                math::four_pi_eps_i * term.crf_2cut3i() * abs2(r);
          const unsigned eg1 = topo.atom_energy_group(*at_it);
          const unsigned eg2 = topo.atom_energy_group(*at2_it);
          ref.crf_energy[eg1][eg2] += e_crf;
        }
      }
    }

    return ref;
  }

  int run_case(const std::string & stopo, const std::string & sconf,
               const std::string & sinput, math::boundary_enum force_boundary,
               const char * label, bool quiet) {
    util::simulation_struct s;
    io::In_Topology in_topo;
    in_topo.quiet = quiet;

    if (util::create_simulation(stopo, "", sconf, sinput, s, in_topo,
                                 "", "", "", "", "", "", "", "", quiet) != 0) {
      std::cerr << label << ": creating simulation failed" << std::endl;
      return 1;
    }
    io::messages.display(std::cout);
    io::messages.clear();

    s.sim.param().pairlist.skin = 0.0;
    s.sim.param().nonbonded.rf_excluded = true;

    s.conf.boundary_type = force_boundary;
    s.sim.param().boundary.boundary = force_boundary;
    if (force_boundary == math::rectangular) {
      s.conf.current().box = math::Box(math::Vec(4.0, 0.0, 0.0),
                                        math::Vec(0.0, 4.0, 0.0),
                                        math::Vec(0.0, 0.0, 4.0));
    }

    const unsigned num_atoms = static_cast<unsigned>(s.topo.num_atoms());

    s.sim.param().pcouple.virial = math::no_virial;

    interaction::Nonbonded_Parameter param;
    in_topo.read_lj_parameter(param.lj_parameter(), std::cout);
    if (flush_messages("lj parameter read")) {
      std::cerr << label << ": reading LJ parameters failed" << std::endl;
      return 1;
    }

    Reference ref;
    if (force_boundary == math::vacuum) {
      ref = compute_reference<math::vacuum>(s.topo, s.conf, s.sim);
    } else {
      ref = compute_reference<math::rectangular>(s.topo, s.conf, s.sim);
    }

    interaction::Pairlist_Algorithm * pa = new interaction::CUDA_Pairlist_Algorithm();
    interaction::CUDA_Nonbonded_Interaction * ni =
        new interaction::CUDA_Nonbonded_Interaction(pa);

    in_topo.read_lj_parameter(ni->parameter().lj_parameter(), std::cout);
    if (flush_messages("lj parameter read (gpu)")) {
      std::cerr << label << ": reading LJ parameters failed (gpu)" << std::endl;
      delete ni;
      return 1;
    }

    if (ni->init(s.topo, s.conf, s.sim, std::cout, quiet) != 0) {
      flush_messages("gpu init");
      std::cerr << label << ": CUDA_Nonbonded_Interaction::init() failed" << std::endl;
      delete ni;
      return 1;
    }

    s.conf.current().force = 0.0;
    s.conf.current().energies.zero();

    if (ni->calculate_interactions(s.topo, s.conf, s.sim) != 0) {
      flush_messages("gpu calculate_interactions");
      std::cerr << label << ": CUDA_Nonbonded_Interaction::calculate_interactions() failed"
                << std::endl;
      delete ni;
      return 1;
    }
    if (flush_messages("gpu calculate_interactions")) {
      delete ni;
      return 1;
    }

    // CUDA_Nonbonded_Interaction writes force directly into the GPU
    // mirror and mark_gpu_dirty()s it (see angle_gpu.t.cc's comment for
    // why this standalone test needs its own explicit publish).
    s.sim.cuda().flush_gpu_dirty(s.conf, gpu::MIRROR_FORCE);
    // Energy is published to conf.old() (not .current()), matching
    // Pressure_Calculation's own convention -- see angle_gpu.t.cc's
    // identical comment.
    s.sim.cuda().flush_gpu_dirty(s.conf, gpu::MIRROR_ENERGY);

    // Isolate the RF-excluded contribution: subtract out the regular
    // pairlist-driven LJ/CRF result (computed the same way
    // cuda_nonbonded_interaction.t.cc does, rf_excluded left off) so this
    // test only checks the new kernel's contribution, not the whole
    // nonbonded force/energy (already covered elsewhere).
    {
      util::simulation_struct s2;
      if (util::create_simulation(stopo, "", sconf, sinput, s2, in_topo,
                                   "", "", "", "", "", "", "", "", true) != 0) {
        std::cerr << label << ": creating baseline simulation failed" << std::endl;
        delete ni;
        return 1;
      }
      io::messages.clear();
      s2.sim.param().pairlist.skin = 0.0;
      // rf_excluded defaults to true (parameter.h's "new standard") --
      // explicitly off here so this baseline run is genuinely the
      // "everything except the new contribution" run being subtracted.
      s2.sim.param().nonbonded.rf_excluded = false;
      s2.conf.boundary_type = force_boundary;
      s2.sim.param().boundary.boundary = force_boundary;
      if (force_boundary == math::rectangular) {
        s2.conf.current().box = s.conf.current().box;
      }
      s2.sim.param().pcouple.virial = math::no_virial;

      interaction::Pairlist_Algorithm * pa2 = new interaction::CUDA_Pairlist_Algorithm();
      interaction::CUDA_Nonbonded_Interaction * ni2 =
          new interaction::CUDA_Nonbonded_Interaction(pa2);
      in_topo.read_lj_parameter(ni2->parameter().lj_parameter(), std::cout);
      io::messages.clear();
      if (ni2->init(s2.topo, s2.conf, s2.sim, std::cout, true) != 0 ||
          ni2->calculate_interactions(s2.topo, s2.conf, s2.sim) != 0) {
        std::cerr << label << ": baseline (rf_excluded=off) GPU run failed" << std::endl;
        io::messages.clear();
        delete ni2;
        delete ni;
        return 1;
      }
      io::messages.clear();
      s2.sim.cuda().flush_gpu_dirty(s2.conf, gpu::MIRROR_FORCE);
      s2.sim.cuda().flush_gpu_dirty(s2.conf, gpu::MIRROR_ENERGY);

      for (unsigned i = 0; i < num_atoms; ++i) {
        s.conf.current().force(i) -= s2.conf.current().force(i);
      }
      const unsigned num_groups = static_cast<unsigned>(s.topo.energy_groups().size());
      for (unsigned gi = 0; gi < num_groups; ++gi) {
        for (unsigned gj = 0; gj < num_groups; ++gj) {
          s.conf.old().energies.crf_energy[gi][gj] -=
              s2.conf.old().energies.crf_energy[gi][gj];
        }
      }
      delete ni2;
    }

    const unsigned num_groups = static_cast<unsigned>(s.topo.energy_groups().size());

    int errors = 0;
    const double tol = 1e-4;

    for (unsigned gi = 0; gi < num_groups; ++gi) {
      for (unsigned gj = 0; gj < num_groups; ++gj) {
        const double gpu_e_crf = s.conf.old().energies.crf_energy[gi][gj];
        const double ref_e_crf = ref.crf_energy[gi][gj];
        if (std::abs(gpu_e_crf - ref_e_crf) > tol * std::max(1.0, std::abs(ref_e_crf))) {
          std::cerr << label << ": e_crf mismatch at group (" << gi << "," << gj
                    << "): gpu=" << gpu_e_crf << " cpu=" << ref_e_crf << std::endl;
          ++errors;
        }
      }
    }
    for (unsigned i = 0; i < num_atoms; ++i) {
      const math::Vec diff = s.conf.current().force(i) - ref.force[i];
      const double scale = std::max(1.0, math::abs(ref.force[i]));
      if (math::abs(diff) > tol * scale) {
        std::cerr << label << ": force mismatch at atom " << i
                  << ": gpu=" << math::v2s(s.conf.current().force(i))
                  << " cpu=" << math::v2s(ref.force[i]) << std::endl;
        ++errors;
      }
    }

    if (errors) {
      std::cerr << label << ": FAILED (" << errors << " mismatch(es))" << std::endl;
    } else {
      std::cout << label << ": OK (" << num_atoms << " atoms, " << num_groups
                << " energy group(s))" << std::endl;
    }

    delete ni;
    return errors;
  }

} // namespace

int main(int argc, char* argv[]) {
#ifdef XXMPI
  MPI_Init(&argc, &argv);
#endif

  util::Known knowns;
  knowns << "verb";

  std::string usage = argv[0];
  usage += "\n\t[@verb   <[module:][submodule:]level>]\n";

  io::Argument args;
  if (args.parse(argc, argv, knowns, true)) {
    std::cerr << usage << std::endl;
    return 1;
  }
  util::parse_verbosity(args);
  const bool quiet = (args.count("verb") == -1);

  std::string stopo, sconf, sinput;
  GETFILEPATH(stopo, "aladip.topo", "src/check/data/");
  GETFILEPATH(sconf, "aladip.conf", "src/check/data/");
  GETFILEPATH(sinput, "aladip_unperturbed.in", "src/check/data/");

  int total = 0;
  total += run_case(stopo, sconf, sinput, math::vacuum, "vacuum", quiet);
  total += run_case(stopo, sconf, sinput, math::rectangular, "rectangular", quiet);

  return total;
}
