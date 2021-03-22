/**
 * @file position_restraint_interaction.cc
 * template methods of Position_Restraint_Interaction
 */
#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"
#include "../../interaction/interaction.h"

#include "../../math/periodicity.h"

// special interactions
#include "../../interaction/interaction_types.h"

#include "../../interaction/special/symmetry_restraint_interaction.h"

#include "../../util/template_split.h"
#include "../../util/debug.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE special

template<math::boundary_enum b>
void _calculate_interactions(topology::Topology& topo,
        configuration::Configuration& conf, simulation::Simulation& sim) {
  math::Periodicity<b> periodicity(conf.current().box);
  
  const std::vector<std::pair<math::Matrix, math::Vec> > & sym = sim.param().symrest.symmetry_operations;
  std::vector<std::pair<math::Matrix, math::Vec> > invsym;
  for(unsigned int i = 0; i < sym.size(); i++) {
    math::Matrix rot(math::inverse(sym[i].first));
    math::Vec trans(math::product(rot, -sym[i].second));
    invsym.push_back(std::pair<math::Matrix, math::Vec>(rot, trans));
  }
  
  conf.current().energies.symrest_total = 0.0;
  
  DEBUG(6, "symmetry restraints");
  std::vector<unsigned int>::const_iterator
  atom_it = topo.sym_restraints().begin(),
          atom_to = topo.sym_restraints().end();

  for (; atom_it != atom_to; ++atom_it) {
    const unsigned int atom_p = *atom_it - topo.sym_asu()[0];
    DEBUG(6, "atom: " << *atom_it);

    switch (sim.param().symrest.symrest) {
      case simulation::xray_symrest_off: break;
      case simulation::xray_symrest_ind:
      {
        DEBUG(6, "symmetry individual atoms");
        for (unsigned int i = 0; i < topo.sym_asu().size() - 1; ++i) {
          const unsigned int atom_i = topo.sym_asu()[i] + atom_p;
          const math::Vec pos_i = math::product(invsym[i].first, conf.current().pos(atom_i)) + invsym[i].second;
          for (unsigned int j = i + 1; j < topo.sym_asu().size(); ++j) {
            const unsigned int atom_j = topo.sym_asu()[j] + atom_p;
            const math::Vec pos_j = math::product(invsym[j].first, conf.current().pos(atom_j)) + invsym[j].second;

            DEBUG(8, "i: " << i << " j: " << j);
            DEBUG(8, "pos i   : " << math::v2s(pos_i));
            DEBUG(8, "pos j   : " << math::v2s(pos_j));

            math::Vec dist;
            periodicity.nearest_image(pos_i, pos_j, dist);
            DEBUG(9, "dist    : " << math::v2s(dist));
            // do the rotation
            const math::Vec r_i = math::product(sym[i].first, dist);
            const math::Vec r_j = math::product(sym[j].first, dist);

            const math::Vec f_i(-sim.param().symrest.force_constant * r_i);
            const math::Vec f_j(sim.param().symrest.force_constant * r_j);

            DEBUG(8, "f i     : " << math::v2s(f_i));
            DEBUG(8, "f j     : " << math::v2s(f_j));

            conf.current().force(atom_i) += f_i;
            conf.current().force(atom_j) += f_j;
            const double V = 0.5 * sim.param().symrest.force_constant * math::abs2(dist);
            conf.current().energies.symrest_total += V;
            DEBUG(7, "energy : " << V);
          }
        } // loop over images
        break;
      }
      case simulation::xray_symrest_constr:
      {
        DEBUG(6, "symmetry constrain atoms");
        // loop over images
        for (unsigned int i = 0; i < topo.sym_asu().size(); ++i) {
          const unsigned int atom_img = topo.sym_asu()[i] + atom_p;
          // optain the image position
          math::Vec pos_img = math::product(sym[i].first, conf.current().pos(atom_img)) +
                  sym[i].second;
          DEBUG(8, "pos     : " << math::v2s(conf.current().pos(atom_img)));
          periodicity.put_into_positive_box(pos_img);
          DEBUG(8, "new pos : " << math::v2s(pos_img));
          conf.current().pos(atom_img) = pos_img;
        } // loop over images
        break;
      }
    } // method switch
  } // loop over restrained atoms

  return;
}

int interaction::Symmetry_Restraint_Interaction::calculate_interactions(topology::Topology& topo,
        configuration::Configuration& conf, simulation::Simulation& sim) {
  m_timer.start();
  SPLIT_BOUNDARY(_calculate_interactions, topo, conf, sim);
  m_timer.stop();
  return 0;
}

int interaction::Symmetry_Restraint_Interaction::init(topology::Topology& topo,
        configuration::Configuration& conf, simulation::Simulation& sim, std::ostream& os, bool quiet) {
  if (!quiet) {
    os << "SYMMETRY RESTRAINS" << std::endl;
    switch (sim.param().symrest.symrest) {
      case simulation::xray_symrest_off:
        os << " - disabled";
        break;
      case simulation::xray_symrest_ind:
        os.precision(4);
        os << " - harmonic on atom position" << std::endl
                << " - force constant: " << std::setw(15) << sim.param().symrest.force_constant;
        break;
      case simulation::xray_symrest_constr:
        os << " - constrained.";
        break;
      default:
        io::messages.add("Method not implemented", "Symmetry_Restraint_Interaction",
                io::message::error);
        return 1;

    }
    os << std::endl;
    os << " - " << sim.param().symrest.symmetry_operations.size() << " symmetry operations:" << std::endl;
    os.precision(4);
    for (unsigned int i = 0; i < sim.param().symrest.symmetry_operations.size(); ++i) {
      for (unsigned int row = 0; row < 3; ++row) {
        os << "      [";
        for (unsigned int col = 0; col < 3; ++col) {
          // looks bad due to numerics
          double val = sim.param().symrest.symmetry_operations[i].first(row, col);
          if (val < math::epsilon) val = 0.0;
          os << std::setw(8) << val;
        }
        os << "] ";
        os << ((row == 1) ? "* ri +" : "      ");
        double val = sim.param().symrest.symmetry_operations[i].second[row];
        if (val < math::epsilon) val = 0.0;
        os << " [" << std::setw(8) << val << "]" << std::endl;
      }
      os << std::endl;
    }
    os << " - restrained atoms (excluding images): " << std::endl;
    for (unsigned int i = 0; i < topo.sym_restraints().size();) {
      os << std::setw(5) << topo.sym_restraints()[i] + 1;
      if (++i % 10 == 0)
        os << std::endl;
    }
    os << std::endl;
    os << "END" << std::endl;
  }
  return 0;
}
