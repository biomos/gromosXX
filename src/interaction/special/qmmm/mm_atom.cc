/**
 * @file mm_atom.cc implements inclusion of MM atoms
 */

#include "../../../stdheader.h"

#include "../../../algorithm/algorithm.h"
#include "../../../topology/topology.h"
#include "../../../simulation/simulation.h"
#include "../../../configuration/configuration.h"
#include "../../../interaction/interaction.h"

#include "../../../math/periodicity.h"
#include "../../../util/template_split.h"

// special interactions
#include "../../../interaction/special/qmmm/mm_atom.h"

template<math::boundary_enum b>
void determine_mm_atoms_by_cutoff(topology::Topology& topo, 
        configuration::Configuration& conf, 
        simulation::Simulation& sim, 
        std::vector<interaction::MM_Atom>& mm_atoms);

void interaction::determine_mm_atoms(topology::Topology& topo, 
        configuration::Configuration& conf, 
        simulation::Simulation& sim, 
        std::vector<MM_Atom>& mm_atoms) {
  mm_atoms.clear();
  if (sim.param().qmmm.cutoff == 0.0) {
    // include all atoms
    for(unsigned int i = 0; i < topo.num_atoms(); ++i) {
      if (topo.in_qm_zone(i)) {
        continue;
      }
      mm_atoms.push_back(interaction::MM_Atom(i, conf.current().pos(i), topo.charge(i)));
    }
  } else {
    SPLIT_BOUNDARY(determine_mm_atoms_by_cutoff, topo, conf, sim, mm_atoms);
  }
}

template<math::boundary_enum b>
void determine_mm_atoms_by_cutoff(topology::Topology& topo, 
        configuration::Configuration& conf, 
        simulation::Simulation& sim, 
        std::vector<interaction::MM_Atom>& mm_atoms) {
  const double cut2 = sim.param().qmmm.cutoff * sim.param().qmmm.cutoff;
  // get the COG for the solute CGs
  math::VArray const &pos = conf.current().pos;
  math::VArray cogs;
  {
    cogs.resize(topo.num_solute_chargegroups());
    // calculate solute center of geometries
    topology::Chargegroup_Iterator cg = topo.chargegroup_begin();
    unsigned int i, num_cg = topo.num_solute_chargegroups();
    for (i = 0; i < num_cg; ++cg, ++i) {
      cg.cog(pos, cogs(i));
    }
  }
  
  std::set<unsigned int> atom_indices;
  math::Periodicity<b> periodicity(conf.current().box);
  const unsigned int num_solute_cg = topo.num_solute_chargegroups(),
          num_cg = topo.num_chargegroups();
  math::Vec r;
  
  // get the atoms close to the QM zone
  for (std::set<topology::qm_atom_struct>::const_iterator
    it = topo.qm_zone().begin(), to = topo.qm_zone().end(); it != to; ++it) {
    const math::Vec & qm_pos = pos(it->index);
    // SOLUTE
    for(unsigned int cg = 0; cg < num_solute_cg; ++cg) {
      // get distance
      periodicity.nearest_image(qm_pos, cogs(cg), r);
      const double r2 = math::abs2(r);
      // in cutoff?
      if (r2 < cut2) {
        for (unsigned int a = topo.chargegroup(cg), a_to = topo.chargegroup(cg + 1);
                a < a_to; ++a) {
          if (topo.in_qm_zone(a)) continue;
          atom_indices.insert(a);
        } // for atoms in cg
      }
    } // for solute cgs
    // SOLVENT
    for(unsigned int cg = num_solute_cg; cg < num_cg; ++cg) {
      // get distance
      periodicity.nearest_image(qm_pos, pos(topo.chargegroup(cg)), r);
      const double r2 = math::abs2(r);
      // in cutoff?
      if (r2 < cut2) {
        for (unsigned int a = topo.chargegroup(cg), a_to = topo.chargegroup(cg + 1);
                a < a_to; ++a) {
          if (topo.in_qm_zone(a)) continue;
          atom_indices.insert(a);
        } // for atoms in cg
      }
    } // for solvent cgs
  } // for QM zone
  
  // create the MM atoms
  mm_atoms.reserve(atom_indices.size());
  for(std::set<unsigned int>::const_iterator it = atom_indices.begin(),
          to = atom_indices.end(); it != to; ++it) {
    mm_atoms.push_back(interaction::MM_Atom(*it, pos(*it), topo.charge(*it)));
  }
}


