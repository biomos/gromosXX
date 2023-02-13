/**
 * @file nonbonded_interaction.cc
 * template methods of Nonbonded_Interaction.
 */
#include "../../../stdheader.h"

#include "../../../algorithm/algorithm.h"
#include "../../../topology/topology.h"
#include "../../../simulation/simulation.h"
#include "../../../configuration/configuration.h"

#include "../../../simulation/parameter.h"

#include "../../../interaction/interaction.h"
#include "../../../interaction/interaction_types.h"
#include "../../../interaction/nonbonded/interaction/nonbonded_parameter.h"

#include "../../../interaction/nonbonded/pairlist/pairlist.h"
#include "../../../interaction/nonbonded/pairlist/pairlist_algorithm.h"

#include "../../../interaction/nonbonded/interaction/storage.h"

#include "../../../interaction/nonbonded/interaction/nonbonded_outerloop.h"
#include "../../../interaction/nonbonded/interaction/nonbonded_set_interface.h"
#include "../../../interaction/nonbonded/interaction/nonbonded_set.h"

#include "../../../interaction/nonbonded/interaction/nonbonded_term.h"
#include "../../../interaction/nonbonded/interaction/perturbed_nonbonded_term.h"
#include "../../../interaction/nonbonded/interaction/eds_nonbonded_term.h"

#include "../../../interaction/nonbonded/interaction/perturbed_nonbonded_pair.h"
#include "../../../interaction/nonbonded/interaction/perturbed_nonbonded_outerloop.h"
#include "../../../interaction/nonbonded/interaction/eds_nonbonded_outerloop.h"

#include "../../../interaction/nonbonded/interaction/perturbed_nonbonded_set.h"
#include "../../../interaction/nonbonded/interaction/eds_nonbonded_set.h"

#include "../../../interaction/nonbonded/interaction/nonbonded_interaction.h"

#include "../../../util/debug.h"

#include "../../../math/periodicity.h"
#include "../../../math/boundary_checks.h"
#include "../../../util/template_split.h"

#ifdef OMP
#include <omp.h>
#endif

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE nonbonded

/**
 * Constructor.
 */
interaction::Nonbonded_Interaction::Nonbonded_Interaction(Pairlist_Algorithm *pa)
: Interaction("NonBonded"),
m_pairlist_algorithm(pa),
m_parameter(),
m_set_size(1),
m_exp_conf(NULL) {
  m_pairlist_algorithm->timer_pointer(&m_timer);
}

/**
 * Destructor.
 * @bug change destruction of nonbonded set to be standard - conform!
 */
interaction::Nonbonded_Interaction::~Nonbonded_Interaction() {
  DEBUG(7, "Nonbonded_Interaction::destructor");
  delete m_pairlist_algorithm;
  DEBUG(12, "pairlist algorithm destroyed");

  for (unsigned int i = 0; i < m_nonbonded_set.size(); ++i) {
    DEBUG(12, "deleting set " << i);
    delete m_nonbonded_set[i];
    m_nonbonded_set[i] = NULL;
  }

  if (m_exp_conf != NULL)
    delete m_exp_conf;
}

/**
 * calculate nonbonded forces and energies.
 */
int interaction::Nonbonded_Interaction::
calculate_interactions(topology::Topology & topo,
        configuration::Configuration & conf,
        simulation::Simulation & sim) {
  DEBUG(4, "Nonbonded_Interaction::calculate_interactions");

  m_timer.start(sim);

  // check if we want to calculate nonbonded
  // might not be necessary if multiple time-stepping is enabled

  int steps = sim.param().multistep.steps;
  if (steps == 0) steps = 1;

  // std::cerr << "Nonbonded: steps = " << steps << std::endl;
  configuration::Configuration *p_conf = &conf;
  topology::Topology *p_topo = &topo;
  if ((sim.steps() % steps) == 0) {

    // std::cerr << "\tMULTISTEP: full non-bonded calculation" << std::endl;

    ////////////////////////////////////////////////////
    // multiple unit cell
    ////////////////////////////////////////////////////
    if (sim.param().multicell.multicell) {
      DEBUG(6, "nonbonded: MULTICELL");
      p_conf = m_exp_conf;
      p_topo = &topo.multicell_topo();
      expand_configuration(topo, conf, sim, *p_conf);
      if (!math::boundary_check_cutoff(p_conf->current().box, p_conf->boundary_type,
          sim.param().pairlist.cutoff_long)) {
        io::messages.add("box is too small: not twice the cutoff!",
                "configuration", io::message::error);
        return 1;
      }
      DEBUG(6, "\tmulticell conf: pos.size()=" << p_conf->current().pos.size());
    }

    // shared memory do this only once
    if (m_pairlist_algorithm->prepare(*p_topo, *p_conf, sim))
      return 1;

    // have to do all from here (probably it's only one,
    // but then maybe it's clearer like it is...)
    for (int i = 0; i < m_set_size; ++i) {
      if(m_nonbonded_set[i]->calculate_interactions(*p_topo, *p_conf, sim))
	return 1;
    }

    ///////////////////////////////////////////////////
    // end of multiple time stepping: calculate
    ////////////////////////////////////////////////////
  } else {
    // std::cerr << "\tMULTISTEP: no non-bonded calculation" << std::endl;
  }

  DEBUG(6, "sets are done, adding things up...");
  store_set_data(*p_topo, *p_conf, sim);
  
  if (sim.param().multicell.multicell) {
    reduce_configuration(topo, conf, sim, *p_conf);
  }

  ////////////////////////////////////////////////////
  // printing pairlist
  ////////////////////////////////////////////////////
  if (sim.param().pairlist.print &&
      (!(sim.steps() % sim.param().pairlist.skip_step))) {
    DEBUG(7, "print pairlist...");
    std::cerr << "printing pairlist!" << std::endl;
    print_pairlist(*p_topo, *p_conf, sim);
  }

  DEBUG(6, "Nonbonded_Interaction::calculate_interactions done");
  m_timer.stop();

  return 0;
}

/**
 * calculate the interaction for a given atom pair.
 * SLOW! as it has to create the periodicity...
 */
int interaction::Nonbonded_Interaction::calculate_interaction
(
        topology::Topology & topo,
        configuration::Configuration & conf,
        simulation::Simulation & sim,
        unsigned int atom_i, unsigned int atom_j,
        math::Vec & force,
        double &e_lj, double &e_crf
        ) {
  assert(m_nonbonded_set.size() >= 1);
  return m_nonbonded_set[0]->calculate_interaction(topo, conf, sim,
          atom_i, atom_j,
          force, e_lj, e_crf);
}

/**
 * calculate the hessian for a given atom.
 */
int interaction::Nonbonded_Interaction::calculate_hessian
(
        topology::Topology & topo,
        configuration::Configuration & conf,
        simulation::Simulation & sim,
        unsigned int atom_i, unsigned int atom_j,
        math::Matrix & hessian
        ) {
  std::vector<Nonbonded_Set_Interface *>::iterator
  it = m_nonbonded_set.begin(),
          to = m_nonbonded_set.end();

  hessian = 0.0;
  math::Matrix h;

  for (; it != to; ++it) {
    (*it)->calculate_hessian(topo, conf, sim, atom_i, atom_j, h);

    for (unsigned int d1 = 0; d1 < 3; ++d1)
      for (unsigned int d2 = 0; d2 < 3; ++d2)
        hessian(d1, d2) += h(d1, d2);
  }
  return 0;
}

/**
 * initialize the arrays
 */
int interaction::Nonbonded_Interaction::init(topology::Topology & topo,
        configuration::Configuration & conf,
        simulation::Simulation & sim,
        std::ostream & os,
        bool quiet) {
  if (!quiet)
    os << "NONBONDED INTERACTION\n";

  configuration::Configuration * p_conf = &conf;
  topology::Topology * p_topo = &topo;
  if (sim.param().multicell.multicell) {
    DEBUG(6, "nonbonded init: MULTICELL");
    if (!quiet)
      os << "\tperforming MULTICELL.\n";
    p_topo = &topo.multicell_topo();
    m_exp_conf = p_conf = new configuration::Configuration();
    init_expand_configuration(topo, conf, sim, *p_conf);
    expand_configuration(topo, conf, sim, *p_conf);
    DEBUG(7, "\tmulticell conf: pos.size()=" << p_conf->current().pos.size());

    // do some important multicell tests:
    // No surface term!
    if (sim.param().nonbonded.ls_epsilon != 0.0) {
      io::messages.add("MULTICELL simulation does only work under tinfoil "
              "boundary conditions (LSEPS=0.0)", "mutlicell", io::message::error);
      return 1;
    }

    // exclusions to other molecules can cause problems because they are only between
    // single copies of the box. Because we gather by molecules this can in rare
    // cases be a problem
    for(topology::Molecule_Iterator m_it = topo.molecule_begin(); m_it != topo.molecule_end(); ++m_it) {
      topology::Atom_Iterator a_it = m_it.begin(), a_start = m_it.begin(), a_end = m_it.end();
      DEBUG(10, "molecule from " << *a_start << " to " << *a_end);
      for(; a_it != a_end; ++a_it) {
        topology::excl_cont_t::value_type::const_iterator ex_it = topo.exclusion(*a_it).begin(),
                ex_to = topo.exclusion(*a_it).end();
        for (; ex_it != ex_to; ++ex_it) {
          if (*ex_it < int(*a_start) || *ex_it >= int(*a_end)) {
            io::messages.add("Exclusions are not in same molecule. This is can cause errors in multicell simulations.",
                    "multicell", io::message::error);
            return 1;
          }
        } // for exclusions
        topology::excl_cont_t::value_type::const_iterator of_it = topo.one_four_pair(*a_it).begin(),
                of_to = topo.one_four_pair(*a_it).end();
        for (; of_it != of_to; ++of_it) {
          if (*of_it < int(*a_start) || *of_it >= int(*a_end)) {
            io::messages.add("1-4 pairs are not in same molecule. This is can cause errors in multicell simulations.",
                    "multicell", io::message::error);
            return 1;
          }
        } // for 1-4 pairs
      }// for atoms
    } // for molecles
  } // check multicell stuff

  if (!math::boundary_check_cutoff(p_conf->current().box, p_conf->boundary_type,
            sim.param().pairlist.cutoff_long)) {
    io::messages.add("box is too small: not twice the cutoff!",
            "configuration", io::message::error);
    return 1;
  }

  // initialise the pairlist...
  m_pairlist_algorithm->init(*p_topo, *p_conf, sim, os, quiet);

  if (sim.param().nonbonded.method != simulation::el_reaction_field) {
    if (!quiet)
      os << "\tlattice-sum electrostatics\n";
    p_conf->lattice_sum().init(*p_topo, sim);
  }

  DEBUG(15, "nonbonded_interaction::initialize");
  m_nonbonded_set.clear();

  if (sim.param().perturbation.perturbation) {
    for (int i = 0; i < m_set_size; ++i) {
      m_nonbonded_set.push_back(new Perturbed_Nonbonded_Set(*m_pairlist_algorithm,
            m_parameter, i, m_set_size));
    }
  } else if (sim.param().eds.eds) {
    DEBUG(16, "creating EDS nonbonded set");
    for (int i = 0; i < m_set_size; ++i) {
      m_nonbonded_set.push_back(new Eds_Nonbonded_Set(*m_pairlist_algorithm,
              m_parameter, i, m_set_size));
      DEBUG(16, "pushed back EDS nonbonded set");
    }
  } else {
    for (int i = 0; i < m_set_size; ++i) {
      m_nonbonded_set.push_back(new Nonbonded_Set(*m_pairlist_algorithm,
              m_parameter, i, m_set_size));
    }
  }

  std::vector<Nonbonded_Set_Interface *>::iterator
  it = m_nonbonded_set.begin(),
          to = m_nonbonded_set.end();

  if (!quiet)
    os << "\tcreated " << m_nonbonded_set.size() << " set"
        <<  (m_nonbonded_set.size() > 1 ? "s" : "")    << "\n";

  bool q = quiet;
  for (; it != to; ++it) {
      (*it)->init(*p_topo, *p_conf, sim, os, q);
    // only print first time...
    q = true;
  }

  if (check_special_loop(*p_topo, *p_conf, sim, os, quiet) != 0) {
    io::messages.add("special solvent loop check failed", "Nonbonded_Interaction",
            io::message::error);
    return 1;
  }
  if (!quiet)
    os << "END\n";
  DEBUG(9, "nonbonded init done");
  return 0;
}


//***************************************************************************
// helper functions 
//***************************************************************************

int interaction::Nonbonded_Interaction::check_special_loop
(
        topology::Topology const & topo,
        configuration::Configuration const & conf,
        simulation::Simulation & sim,
        std::ostream & os,
        bool quiet) {
  DEBUG(7, "checking for spc interaction loops");
  DEBUG(10, " param: special_loop = " << sim.param().force.special_loop);

  if (sim.param().innerloop.method == simulation::sla_off) {
    sim.param().force.special_loop = simulation::special_loop_off;
    DEBUG(8, "standard loops, user request");
    if (!quiet)
      os << "\tusing standard solvent loops (user request)\n";
    return 0;
  }

  if (sim.param().pairlist.atomic_cutoff) {
    if (!quiet)
      os << "\tusing standard solvent loops (atomic cutoff)\n";
    return 1;
  }

  DEBUG(10, "num_solvents = " << topo.num_solvents());
  if (topo.num_solvents() > 0) {
    DEBUG(10, "molecules = " << topo.num_solvent_molecules(0));
    DEBUG(10, "atoms = " << topo.num_solvent_atoms(0));
  }

  // check whether there is a solvent.
  if (topo.num_solvents() != 1 || topo.num_solvent_molecules(0) < 1) {

    DEBUG(10, "standard loops...");

    if (!quiet) {
      if (topo.num_solvents() > 0) {
        if (topo.num_solvent_molecules(0) == 0) {
          os << "\tusing standard solvent loops (no solvent present!)\n";
        }
      } else {
        os << "\tusing standard solvent loops (no solvent in topology!)\n";
      }
    }
    return 1;
  }

  // check energy groups
  DEBUG(10, "checking energy groups ...");

  if (topo.atom_energy_group(topo.num_solute_atoms()) !=
      topo.atom_energy_group(topo.num_atoms() - 2)) {
    DEBUG(10, "- incompatible.");
    if (!quiet)
      os << "\tusing standard solvent loops (energy group partitioning incompatible).\n"
            << "\tAll solvent atoms must be in one single energy group. ";
    return 1;
  }

  if (sim.param().innerloop.method == simulation::sla_generic) {
    // the generic loop works for now.
    DEBUG(8, "generic solvent loop, user request");
    sim.param().force.special_loop = simulation::special_loop_generic;
    if (!quiet)
      os << "\tusing generic solvent loops (user request)\n";
    return 0;
  }

  if (sim.param().innerloop.method == simulation::sla_cuda) {
    // the generic loop works for now.
    DEBUG(8, "CUDA GPU solvent loop, user request");
    if (!quiet)
      os << "\tusing CUDA GPU solvent loops (user request)\n";
    return 0;
  }

  // check whether the number of atoms match for other methods
  switch (sim.param().innerloop.solvent) {
    case simulation::sls_spc :
    {
      if (topo.num_solvent_atoms(0) / topo.num_solvent_molecules(0) != 3) {
        os << "\tusing standard solvent loops (num solvents doesn't match)\n"
                << "\t\tnum solvents: "
                << topo.num_solvents() << "\n"
                << "\t\tsolvent atoms: "
                << topo.num_solvent_atoms(0) / topo.num_solvent_molecules(0) << "\n"
                << "\t\tmolecules: " << topo.num_solvent_molecules(0) << "\n";
        return 1;
      }
      break;
    }
    default:
      os << "\tusing standard solvent loops (unknown solvent)\n";
      return 1;
  }

  // check the hardcoded parameters if the method applies
  if (sim.param().innerloop.method == simulation::sla_hardcode) {
    // check four_pi_eps_i
    DEBUG(10, "checking (4 pi eps0)^-1 ...");
    if (math::four_pi_eps_i != 138.9354) {
      DEBUG(10, " does not match, force standard loop");
      if (!quiet)
        os << "\tusing standard solvent loops ((4 pi eps0)^-1 does not match)\n";
      return 1;
    }
    // check solvent specific parameters
    switch (sim.param().innerloop.solvent) {
      case simulation::sls_spc :
      {
        DEBUG(10, "checking charges...");

        // check charges
        if (topo.charge()(topo.num_solute_atoms()) != -0.82 ||
            topo.charge()(topo.num_solute_atoms() + 1) != 0.41 ||
            topo.charge()(topo.num_solute_atoms() + 2) != 0.41) {

          DEBUG(10, "charges don't match, standard loops");
          if (!quiet)
            os << "\tusing standard solvent loops (charges don't match)\n"
                  << "\t\tO  : " << topo.charge()(topo.num_solute_atoms()) << "\n"
            << "\t\tH1 : " << topo.charge()(topo.num_solute_atoms() + 1) << "\n"
            << "\t\tH2 : " << topo.charge()(topo.num_solute_atoms() + 2) << "\n";

          return 1;
        }

        // check lj parameters
        DEBUG(10, "checking LJ parameter...");
        const lj_parameter_struct &lj_OO =
                m_parameter.lj_parameter(topo.iac(topo.num_solute_atoms()),
                topo.iac(topo.num_solute_atoms()));

        const lj_parameter_struct &lj_OH1 =
                m_parameter.lj_parameter(topo.iac(topo.num_solute_atoms()),
                topo.iac(topo.num_solute_atoms() + 1));

        const lj_parameter_struct &lj_OH2 =
                m_parameter.lj_parameter(topo.iac(topo.num_solute_atoms()),
                topo.iac(topo.num_solute_atoms() + 2));

        const lj_parameter_struct &lj_H1H2 =
                m_parameter.lj_parameter(topo.iac(topo.num_solute_atoms() + 1),
                topo.iac(topo.num_solute_atoms() + 2));

        if (lj_OO.c6 != 2.617346E-3 ||
            lj_OO.c12 != 2.634129E-6 ||
            lj_OH1.c6 != 0.0 ||
            lj_OH1.c12 != 0.0 ||
            lj_OH2.c6 != 0.0 ||
            lj_OH2.c12 != 0.0 ||
            lj_H1H2.c6 != 0.0 ||
            lj_H1H2.c12 != 0.0) {

          DEBUG(10, "don't match, force standard loop");
          if (!quiet)
            os << "\tusing standard solvent loops (van der Waals parameter don't match)\n";
          return 1;
        }

        DEBUG(10, "happy to use spc loops");
        sim.param().force.special_loop = simulation::special_loop_spc;
        if (!quiet) {
          os << "\tusing hardcoded spc solvent loops\n";
        }

        break;
      }
      default:
        os << "\tusing standard solvent loops (unknown solvent)\n";
        return 1;
    } // solvents
  } // hardcoded parameters

  // check the tabulated parameters if the method applies
  if (sim.param().innerloop.method == simulation::sla_table) {
    // check solvent specific parameters
    switch (sim.param().innerloop.solvent) {
      case simulation::sls_spc :
      {
        // load the table (check section)
#define CHECK_PARAM "check_param"
#include "../../../interaction/nonbonded/interaction/spc_table.h"
#undef CHECK_PARAM

        // check four_pi_eps_i
        DEBUG(10, "checking (4 pi eps0)^-1 ...");
        if (math::four_pi_eps_i != four_pi_eps_i) {
          DEBUG(10, " does not match, force standard loop");
          if (!quiet)
            os << "\tusing standard solvent loops ((4 pi eps0)^-1 does not match)\n";
          return 1;
        }

        DEBUG(10, "checking cutoffs and RF parameters...");
        if (shortrange_cutoff != sim.param().pairlist.cutoff_short ||
            longrange_cutoff != sim.param().pairlist.cutoff_long ||
            longrange_cutoff != sim.param().nonbonded.rf_cutoff ||
            solvent_permittivity != sim.param().nonbonded.rf_epsilon ||
            solvent_kappa != sim.param().nonbonded.rf_kappa) {
          DEBUG(10, " do not match, force standard loop");
          if (!quiet)
            os << "\tusing standard solvent loops (RF parameters or cutoffs do not match)\n";
          return 1;
        }

        DEBUG(10, "checking charges...");
        // check charges
        const double qO = topo.charge()(topo.num_solute_atoms());
        const double qH1 = topo.charge()(topo.num_solute_atoms() + 1);
        const double qH2 = topo.charge()(topo.num_solute_atoms() + 2);
        if (qH2 - qH1 > math::epsilon || qH1 * qO - qOH > math::epsilon ||
            qO * qO - qOO > math::epsilon || qH1 * qH1 - qHH > math::epsilon) {
          DEBUG(10, "charges don't match, standard loops");
          if (!quiet)
            os << "\tusing standard solvent loops (charges don't match)\n"
                  << "\t\tO  : " << topo.charge()(topo.num_solute_atoms()) << "\n"
            << "\t\tH1 : " << topo.charge()(topo.num_solute_atoms() + 1) << "\n"
            << "\t\tH2 : " << topo.charge()(topo.num_solute_atoms() + 2) << "\n";

          return 1;
        }

        // check lj parameters
        DEBUG(10, "checking LJ parameter...");
        const lj_parameter_struct &lj_OO =
                m_parameter.lj_parameter(topo.iac(topo.num_solute_atoms()),
                topo.iac(topo.num_solute_atoms()));

        const lj_parameter_struct &lj_OH1 =
                m_parameter.lj_parameter(topo.iac(topo.num_solute_atoms()),
                topo.iac(topo.num_solute_atoms() + 1));

        const lj_parameter_struct &lj_OH2 =
                m_parameter.lj_parameter(topo.iac(topo.num_solute_atoms()),
                topo.iac(topo.num_solute_atoms() + 2));

        const lj_parameter_struct &lj_H1H2 =
                m_parameter.lj_parameter(topo.iac(topo.num_solute_atoms() + 1),
                topo.iac(topo.num_solute_atoms() + 2));

        if (lj_OO.c6 - c6 > math::epsilon ||
            lj_OO.c12 - c12 > math::epsilon ||
            lj_OH1.c6 != 0.0 ||
            lj_OH1.c12 != 0.0 ||
            lj_OH2.c6 != 0.0 ||
            lj_OH2.c12 != 0.0 ||
            lj_H1H2.c6 != 0.0 ||
            lj_H1H2.c12 != 0.0) {

          DEBUG(10, "don't match, force standard loop");
          if (!quiet)
            os << "\tusing standard solvent loops (van der Waals parameter don't match)\n";
          return 1;
        }

        // maybe one want's to check the solvent diameter
        DEBUG(10, "happy to use spc loops");
        sim.param().force.special_loop = simulation::special_loop_spc_table;
        if (!quiet) {
          os << "\tusing tabulated spc solvent loops\n";
          os << "\t\tshortrange table size: " << shortrange_table << "\n";
          os << "\t\tlongrange table size:  " << longrange_table << "\n";
        }

        break;
      }
      default:
        os << "\tusing standard solvent loops (unknown solvent)\n";
        return 1;
    } // solvents
  } // tabulated forces

  return 0;
}

/**
 * store data from sets into the configuration
 */
void interaction::Nonbonded_Interaction::store_set_data
(
        topology::Topology const & topo,
        configuration::Configuration & conf,
        simulation::Simulation const & sim
        ) {
  m_timer.start_subtimer("set data summation");
  std::vector<Nonbonded_Set_Interface *>::iterator
  it = m_nonbonded_set.begin(),
          to = m_nonbonded_set.end();

  // add the forces, energies, virial...
  for (; it != to; ++it) {
    DEBUG(7, "adding forces from set " << it - m_nonbonded_set.begin());
    (*it)->update_configuration(topo, conf, sim);
  }
  m_timer.stop_subtimer("set data summation");
}

void interaction::Nonbonded_Interaction::init_expand_configuration
(
        topology::Topology const & topo,
        configuration::Configuration & conf,
        simulation::Simulation & sim,
        configuration::Configuration & exp_conf) {
  const int mul = sim.param().multicell.x * sim.param().multicell.y * sim.param().multicell.z;
  DEBUG(8, "\tmul=" << mul);

  assert(topo.multicell_topo().num_atoms() == topo.num_atoms() * mul);
  // resize the configuration
  exp_conf.resize(topo.multicell_topo().num_atoms());
  
  exp_conf.boundary_type = conf.boundary_type;
  exp_conf.current().box(0) = sim.param().multicell.x * math::Vec(conf.current().box(0));
  exp_conf.current().box(1) = sim.param().multicell.y * math::Vec(conf.current().box(1));
  exp_conf.current().box(2) = sim.param().multicell.z * math::Vec(conf.current().box(2));

  exp_conf.init(topo.multicell_topo(), sim.param(), false);

  DEBUG(10, "\texp_conf initialised");
}

void interaction::Nonbonded_Interaction::expand_configuration
(
        topology::Topology const & topo,
        configuration::Configuration & conf,
        simulation::Simulation & sim,
        configuration::Configuration & exp_conf) {
  SPLIT_BOUNDARY(_expand_configuration, topo, conf, sim, exp_conf);
}

/**
 * expand a configuration for
 * multiple unit cell
 * simulations
 */
template<math::boundary_enum b>
void interaction::Nonbonded_Interaction::_expand_configuration
(
        topology::Topology const & topo,
        configuration::Configuration & conf,
        simulation::Simulation & sim,
        configuration::Configuration & exp_conf
        ) {

  DEBUG(7, "expanding configuration");

  const int mul = sim.param().multicell.x * sim.param().multicell.y * sim.param().multicell.z;
  DEBUG(8, "\tmul=" << mul);
  
  math::Periodicity<b> periodicity(conf.current().box);
  periodicity.gather_molecules_into_box(conf, topo);

  exp_conf.boundary_type = conf.boundary_type;
  exp_conf.current().box(0) = sim.param().multicell.x * math::Vec(conf.current().box(0));
  exp_conf.current().box(1) = sim.param().multicell.y * math::Vec(conf.current().box(1));
  exp_conf.current().box(2) = sim.param().multicell.z * math::Vec(conf.current().box(2));

  DEBUG(10, "\texp_conf initialised");

  math::Vec shift(0.0);
  unsigned int exp_i = 0;

  DEBUG(10, "Multicell box:" << math::v2s(exp_conf.current().box(0)) << std::endl
          << math::v2s(exp_conf.current().box(1)) << std::endl
          << math::v2s(exp_conf.current().box(2)));

  // SOLUTE
  DEBUG(10, "\tsolute");
  for (int z = 0; z < sim.param().multicell.z; ++z) {
    for (int y = 0; y < sim.param().multicell.y; ++y) {
      for (int x = 0; x < sim.param().multicell.x; ++x) {

        shift = x * conf.current().box(0) +
                y * conf.current().box(1) +
                z * conf.current().box(2);

        // this should be the NORMAL topo!
        for (unsigned int i = 0; i < topo.num_solute_atoms(); ++i, ++exp_i) {

          assert(exp_conf.current().pos.size() > unsigned(exp_i));
          assert(conf.current().pos.size() > i);

          exp_conf.current().pos(exp_i) = conf.current().pos(i) + shift;
          exp_conf.current().posV(exp_i) = conf.current().posV(i);
          DEBUG(10, "i: " << exp_i << " pos: " << math::v2s(exp_conf.current().pos(exp_i)));
          // exp_conf.old().pos(exp_i) = conf.old().pos(i) + shift;
        }
      }
    }
  }

  // SOLVENT
  for (int z = 0; z < sim.param().multicell.z; ++z) {
    for (int y = 0; y < sim.param().multicell.y; ++y) {
      for (int x = 0; x < sim.param().multicell.x; ++x) {

        shift = x * conf.current().box(0) +
                y * conf.current().box(1) +
                z * conf.current().box(2);

        for (unsigned int i = topo.num_solute_atoms(); i < topo.num_atoms(); ++i, ++exp_i) {

          assert(exp_conf.current().pos.size() > unsigned(exp_i));
          assert(conf.current().pos.size() > i);

          exp_conf.current().pos(exp_i) = conf.current().pos(i) + shift;
          exp_conf.current().posV(exp_i) = conf.current().posV(i);

          // exp_conf.old().pos(exp_i) = conf.old().pos(i) + shift;
          DEBUG(10, "i: " << exp_i << " pos: " << math::v2s(exp_conf.current().pos(exp_i)));
        }
      }
    }
  }

  assert(topo.multicell_topo().num_atoms() == exp_i);

  exp_conf.current().force = 0.0;
  exp_conf.current().energies.zero();
  exp_conf.current().perturbed_energy_derivatives.zero();
  exp_conf.current().virial_tensor = 0.0;
}

/**
 * reduce a configuration for
 * multiple unit cell
 * simulations
 */
void interaction::Nonbonded_Interaction::reduce_configuration
(
        topology::Topology const & topo,
        configuration::Configuration & conf,
        simulation::Simulation & sim,
        configuration::Configuration & exp_conf
        ) {

  // add one-four, rf excluded etc... all those things that go directly into 
  // the configuration and not into storages of the sets

  // reduce the forces
  const unsigned int cells = (sim.param().multicell.x * sim.param().multicell.y * sim.param().multicell.z);
  for (unsigned int i = 0; i < topo.num_solute_atoms(); ++i) {
    DEBUG(6, "i: " << i << " nb f: " << math::v2s(exp_conf.current().force(i)));
    conf.current().force(i) += exp_conf.current().force(i);
    DEBUG(6, "i: " << i << " f: " << math::v2s(conf.current().force(i)));
    conf.current().posV(i) = exp_conf.current().posV(i);
  }

  // one cell is is already contained in i!! -> cells - 1
  const unsigned int offset = topo.num_solute_atoms() * (cells - 1);
  for (unsigned int i = topo.num_solute_atoms(); i < topo.num_atoms(); ++i) {
    DEBUG(6, "i: " << i << " nb f: " << math::v2s(exp_conf.current().force(offset + i)));
    conf.current().force(i) += exp_conf.current().force(offset + i);
    DEBUG(6, "i: " << i << " f: " << math::v2s(conf.current().force(i)));
    conf.current().posV(i) = exp_conf.current().posV(offset + i);
  }

  const int ljs = conf.current().energies.lj_energy.size();
  configuration::Energy & e = conf.current().energies;
  configuration::Energy & exp_e = exp_conf.current().energies;
  configuration::Energy & pe = conf.current().perturbed_energy_derivatives;
  configuration::Energy & exp_pe = exp_conf.current().perturbed_energy_derivatives;

  const double cells_i = 1.0 / cells;

  // reduce the energies
  for (int i = 0; i < ljs; ++i) {
    for (int j = 0; j < ljs; ++j) {
      e.lj_energy[i][j] += exp_e.lj_energy[i][j] * cells_i;
      e.crf_energy[i][j] += exp_e.crf_energy[i][j] * cells_i;
      e.shift_extra_orig[i][j] += exp_e.shift_extra_orig[i][j] * cells_i;
      e.shift_extra_phys[i][j] += exp_e.shift_extra_phys[i][j] * cells_i;
      e.ls_real_energy[i][j] += exp_e.ls_real_energy[i][j] * cells_i;
      pe.lj_energy[i][j] += exp_pe.lj_energy[i][j] * cells_i;
      pe.crf_energy[i][j] += exp_pe.crf_energy[i][j] * cells_i;
    }
    e.self_energy[i] += exp_e.self_energy[i] * cells_i;
    pe.self_energy[i] += exp_pe.self_energy[i] * cells_i;
  }

  e.ls_a_term_total += exp_e.ls_a_term_total;
  e.ls_kspace_total += exp_e.ls_kspace_total * cells_i;
  e.ls_pair_total += exp_e.ls_pair_total;
  e.ls_self_total += exp_e.ls_self_total;
  e.ls_surface_total += exp_e.ls_surface_total;
  
  // reduce the virial
  if (sim.param().pcouple.virial) {
    DEBUG(7, "\tadd set virial");

    for (unsigned int i = 0; i < 3; ++i) {
      for (unsigned int j = 0; j < 3; ++j) {

        conf.current().virial_tensor(i, j) +=
                exp_conf.current().virial_tensor(i, j) * cells_i;
      }
    }
  }
}

/**
 * print the pairlist
 */
int interaction::Nonbonded_Interaction::print_pairlist
(
        topology::Topology & topo,
        configuration::Configuration & conf,
        simulation::Simulation & sim,
        std::ostream & os
        ) {
  DEBUG(4, "Nonbonded_Interaction::print_pairlist");

  Pairlist temp_solute, temp_solvent;
  temp_solute.resize(topo.num_atoms());
  temp_solvent.resize(topo.num_atoms());

  for (unsigned int atom_i = 0; atom_i < topo.num_atoms(); ++atom_i) {

    for (int i = 0; i < m_set_size; ++i) {

      assert(m_nonbonded_set.size() > unsigned(i));
      assert(m_nonbonded_set[i]->pairlist().solute_short.size() > atom_i);
      assert(m_nonbonded_set[i]->pairlist().solvent_short.size() > atom_i);

      for (unsigned int atom_j = 0;
              atom_j < m_nonbonded_set[i]->pairlist().solute_short[atom_i].size();
              ++atom_j) {

        assert(temp_solute.size() > atom_i);
        assert(temp_solute.size() > m_nonbonded_set[i]->pairlist().solute_short[atom_i][atom_j]);

        if (m_nonbonded_set[i]->pairlist().solute_short[atom_i][atom_j] < atom_i)
          temp_solute[m_nonbonded_set[i]->pairlist().solute_short[atom_i][atom_j]].push_back(atom_i);
        else
          temp_solute[atom_i].push_back(m_nonbonded_set[i]->pairlist().solute_short[atom_i][atom_j]);
      }
      for (unsigned int atom_j = 0;
              atom_j < m_nonbonded_set[i]->pairlist().solvent_short[atom_i].size();
              ++atom_j) {

        assert(temp_solvent.size() > atom_i);
        assert(temp_solvent.size() > m_nonbonded_set[i]->pairlist().solvent_short[atom_i][atom_j]);

        if (m_nonbonded_set[i]->pairlist().solvent_short[atom_i][atom_j] < atom_i)
          temp_solvent[m_nonbonded_set[i]->pairlist().solvent_short[atom_i][atom_j]].push_back(atom_i);
        else
          temp_solvent[atom_i].push_back(m_nonbonded_set[i]->pairlist().solvent_short[atom_i][atom_j]);
      }
    }
  }

  os << temp_solute << std::endl
          << temp_solvent << std::endl;
  return 0;
}
