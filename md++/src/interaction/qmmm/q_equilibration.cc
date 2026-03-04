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

#include "../../../stdheader.h"

#include "../../../algorithm/algorithm.h"      // <-- add this

#include "../../../topology/topology.h"
#include "../../../simulation/simulation.h"
#include "../../../configuration/configuration.h"

#include "../../../io/blockinput.h"
#include "../../../util/timing.h"
#include "../../../util/system_call.h"
#include "../../../util/debug.h"

#include "qm_zone.h"
#include "qm_atom.h"
#include "qm_link.h"
#include "mm_atom.h"
#include "q_equilibration.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE qmmm

namespace interaction {
namespace qeq {

// --- QEq parameters: simple element-based defaults (PLACEHOLDERS) ---
// Units: must be consistent with your electrostatics energy units.
// For now: treat chi and Jii as dimensionless in the same "energy per unit charge" system
// that matches phi. You will replace/tune these later.
static inline double qeq_chi_default(const unsigned Z)
{
  // Very rough placeholders (do NOT treat as physical)
  // You will later fit these to your reference charges.
  switch (Z) {
    case 1:  return 436.886;  // H
    case 6:  return 515.521;  // C
    case 7:  return 665.652;  // N
    case 8:  return 843.378;  // O
    case 16: return 668.450;  // S
    default: return 0.00;  // generic fallback
  }
}

static inline double qeq_Jii_default(const unsigned Z)
{
  // Must be > 0. Larger => stiffer charges (smaller response).
  // Again placeholders; you will later fit these.
  switch (Z) {
    case 1:  return 1340.220;  // H
    case 6:  return 977.010;  // C
    case 7:  return 1134.668;  // N
    case 8:  return 1289.430;  // O
    case 16: return 865.666;  // S
    default: return 1.20;  // generic fallback
  }
}

// --- helper: pair contribution (q/r) including COS, matching mopac_worker.cc logic ---
static inline double pair_potential_q_over_r(const math::Vec& qm_pos,
                                            const interaction::MM_Atom& mm_atom)
{
  if (mm_atom.is_polarisable) {
    // MM charge part
    double pot = (mm_atom.charge - mm_atom.cos_charge) / math::abs(qm_pos - mm_atom.pos);
    // COS charge part
    pot += mm_atom.cos_charge / math::abs(qm_pos - mm_atom.pos - mm_atom.cosV);
    return pot;
  } else {
    return mm_atom.charge / math::abs(qm_pos - mm_atom.pos);
  }
}

// --- helper: build excluded MM atoms for a linked QM atom (mode 1/2) ---
static inline void build_excluded_mm_atoms(const topology::Topology& topo,
                                          const simulation::Simulation& sim,
                                          const interaction::QM_Zone& qm_zone,
                                          const interaction::QM_Atom& qm_atom,
                                          std::set<unsigned>& excluded_mm)
{
  excluded_mm.clear();

  // Exclude directly linked MM atom(s)
  for (std::set<interaction::QM_Link>::const_iterator
         li_it = qm_zone.link.begin(), li_to = qm_zone.link.end();
       li_it != li_to; ++li_it) {
    if (li_it->qm_index == qm_atom.index) {
      excluded_mm.insert(li_it->mm_index);
      DEBUG(15, "QEq: excluding directly linked MM atom " << li_it->mm_index);
    }
  }

  // If mode 2: exclude the entire linked charge group(s)
  if (sim.param().qmmm.mopac.link_atom_mode == 2) {
    std::set<unsigned> excluded_cgs;

    for (std::set<unsigned>::const_iterator it = excluded_mm.begin(), to = excluded_mm.end();
         it != to; ++it) {

      int cg_idx = -1;
      for (unsigned cg = 0; cg < topo.num_chargegroups(); ++cg) {
        if (*it < unsigned(topo.chargegroup(cg + 1))) {
          cg_idx = int(cg);
          break;
        }
      }
      assert(cg_idx != -1);
      excluded_cgs.insert(unsigned(cg_idx));
    }

    // translate chargegroup indices into atom indices
    for (std::set<unsigned>::const_iterator it = excluded_cgs.begin(), to = excluded_cgs.end();
         it != to; ++it) {
      for (int a = topo.chargegroup(*it); a < topo.chargegroup(*it + 1); ++a) {
        excluded_mm.insert(unsigned(a));
      }
    }
  }
} 

// --- helper: total potential (sum q/r) at one QM atom, with optional link exclusions ---
static inline double total_potential_q_over_r(const topology::Topology& topo,
                                             const simulation::Simulation& sim,
                                             const interaction::QM_Zone& qm_zone,
                                             const interaction::QM_Atom& qm_atom)
{
  double pot = 0.0;

  // Same semantics as mopac_worker.cc:
  // 0: link atoms see no MM atoms
  // 1: exclude linked MM atom
  // 2: exclude linked charge group
  // 3: link atoms see all MM atoms
  if (!qm_atom.is_linked) {
    for (std::set<interaction::MM_Atom>::const_iterator
           mm_it = qm_zone.mm.begin(), mm_to = qm_zone.mm.end();
         mm_it != mm_to; ++mm_it) {
      pot += pair_potential_q_over_r(qm_atom.pos, *mm_it);
    }
    return pot;
  }

  if (sim.param().qmmm.mopac.link_atom_mode == 0) {
    return 0.0;
  }

  std::set<unsigned> excluded_mm;
  build_excluded_mm_atoms(topo, sim, qm_zone, qm_atom, excluded_mm);

  for (std::set<interaction::MM_Atom>::const_iterator
         mm_it = qm_zone.mm.begin(), mm_to = qm_zone.mm.end();
       mm_it != mm_to; ++mm_it) {
    if (excluded_mm.find(mm_it->index) == excluded_mm.end()) {
      pot += pair_potential_q_over_r(qm_atom.pos, *mm_it);
    }
  }

  return pot;
} // total_potential_q_over_r

int update_qm_charges_by_qeq(topology::Topology& topo,
                            const configuration::Configuration& conf,
                            const simulation::Simulation& sim,
                            interaction::QM_Zone& qm_zone)
{

  if (qm_zone.mm.empty()) {
    // This means get_mm_atoms() did not gather anything (e.g. cutoff too small,
    // or logic still gated somewhere).
    DEBUG(1, "QEq: WARNING: qm_zone.mm still empty after ensure_mm_atoms; phi^OR will be zero.");
  }

  // --- STEP 1: compute phi_i^OR at each QM atom i in INTERNAL/MM units ---
  // phi_i = (1/4*pi*eps0) * sum_j q_j / r_ij  (including COS as in MOPAC worker)
  const double prefactor = math::four_pi_eps_i;

  // Store potentials in a temporary map keyed by QM atom index
  std::map<unsigned, double> phi_on_qm;

  DEBUG(10, "QEq: qm_zone.mm size = " << qm_zone.mm.size()
          << " qm_zone.qm size = " << qm_zone.qm.size()
          << " link size = " << qm_zone.link.size());

  for (std::set<interaction::QM_Atom>::iterator
         qm_it = qm_zone.qm.begin(), qm_to = qm_zone.qm.end();
       qm_it != qm_to; ++qm_it) {

    const double sum_q_over_r = total_potential_q_over_r(topo, sim, qm_zone, *qm_it);
    const double phi = prefactor * sum_q_over_r;

    phi_on_qm[qm_it->index] = phi;

    DEBUG(8, "QEq: phi^OR on QM atom " << qm_it->index
              << " (Z=" << qm_it->atomic_number << ") = " << phi << " four_pi_eps_i = " << prefactor);
  }

  // --- STEP 2: diagonal-only QEq solve with charge constraint ---
  // Minimize:
  //   E(q) = sum_i (chi_i q_i + 1/2 J_ii q_i^2) + sum_i q_i * phi_i
  // subject to:
  //   sum_i q_i = Q_IR
  //
  // Stationarity gives:
  //   q_i = -(chi_i + phi_i + lambda) / J_ii
  // and lambda from the constraint:
  //   lambda = - ( Q_IR + sum_i (chi_i + phi_i)/J_ii ) / sum_i (1/J_ii)

  const double Q_IR = double(qm_zone.charge()); // net charge of the QM/IR zone

  double sum_invJ = 0.0;
  double sum_chi_phi_over_J = 0.0;
  unsigned n_qeq = 0;

  // First pass: accumulate sums over QM atoms that actually get QEq charges
  for (std::set<interaction::QM_Atom>::iterator
         qm_it = qm_zone.qm.begin(), qm_to = qm_zone.qm.end();
       qm_it != qm_to; ++qm_it) {

    // Only assign to actual QM atoms (matches your existing application logic)
    //if (!topo.is_qm(qm_it->index)) continue;

    const unsigned Z = qm_it->atomic_number;
    const double chi = qeq_chi_default(Z);
    const double Jii = qeq_Jii_default(Z);

    const double phi = phi_on_qm[qm_it->index];

    sum_invJ += 1.0 / Jii;
    sum_chi_phi_over_J += (chi + phi) / Jii;
    ++n_qeq;
  }

  if (n_qeq == 0) {
    DEBUG(10, "QEq: no QM atoms found for QEq charge assignment (n_qeq=0).");
    return 0;
  }

  const double lambda = -(Q_IR + sum_chi_phi_over_J) / sum_invJ;

  DEBUG(10, "QEq: solved lambda = " << lambda
            << " with Q_IR=" << Q_IR
            << " n_qeq=" << n_qeq);

  // Second pass: compute charges and store into qm_charge
  double q_sum_check = 0.0;

  for (std::set<interaction::QM_Atom>::iterator
         qm_it = qm_zone.qm.begin(), qm_to = qm_zone.qm.end();
       qm_it != qm_to; ++qm_it) {

    //if (!topo.is_qm(qm_it->index)) continue;

    const unsigned Z = qm_it->atomic_number;
    const double chi = qeq_chi_default(Z);
    const double Jii = qeq_Jii_default(Z);
    const double phi = phi_on_qm[qm_it->index];

    const double q = -(chi + phi + lambda) / Jii;

    qm_it->qm_charge = q;
    q_sum_check += q;

    DEBUG(10, "QEq: QM atom " << qm_it->index
              << " Z=" << Z
              << " phi=" << phi
              << " chi=" << chi
              << " Jii=" << Jii
              << " q=" << q);
  }

  DEBUG(10, "QEq: charge constraint check: sum(q)=" << q_sum_check
            << " target Q_IR=" << Q_IR);


  return 0;
}

} // namespace qeq
} // namespace interaction
