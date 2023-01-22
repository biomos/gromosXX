/**
 * @file qmmm_nonbonded_outerloop.cc
 * template methods of QMMM_Nonbonded_Outerloop.
 */


#include "../../../stdheader.h"

#include "../../../algorithm/algorithm.h"
#include "../../../topology/topology.h"
#include "../../../simulation/simulation.h"
#include "../../../configuration/configuration.h"
#include "../../../interaction/interaction.h"

#include "../../../math/periodicity.h"
#include "../../../math/volume.h"

#include "../../../interaction/interaction_types.h"
#include "../../../interaction/nonbonded/innerloop_template.h"
#include "../../../interaction/nonbonded/pairlist/pairlist.h"
#include "../../../interaction/nonbonded/interaction/storage.h"
#include "../../../interaction/nonbonded/interaction/nonbonded_parameter.h"

#include "../../../interaction/nonbonded/interaction/nonbonded_term.h"
#include "../../../interaction/nonbonded/interaction/nonbonded_innerloop.h"

#include "../../../interaction/nonbonded/interaction_spec.h"

#include "../../../util/debug.h"

#include "../../../simulation/parameter.h"

#include "qmmm_nonbonded_outerloop.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE qmmm

/**
 * Constructor.
 */
interaction::QMMM_Nonbonded_Outerloop
::QMMM_Nonbonded_Outerloop(Nonbonded_Parameter &nbp)
: Nonbonded_Outerloop(nbp) {
}

//==================================================
// interaction loops
//==================================================

void interaction::QMMM_Nonbonded_Outerloop
::lj_outerloop(topology::Topology & topo,
        configuration::Configuration & conf,
        simulation::Simulation & sim,
        Pairlist const & pairlist_solute,
        Pairlist const & pairlist_solvent,
        Storage & storage,
        bool longrange, util::Algorithm_Timer & timer, bool master) {
  SPLIT_MY_INNERLOOP(_lj_outerloop, simulation::lj_func, topo, conf, sim,
          pairlist_solute, pairlist_solvent, storage, longrange, timer, master);
}

template<typename t_interaction_spec>
void interaction::QMMM_Nonbonded_Outerloop
::_lj_outerloop(topology::Topology & topo,
                configuration::Configuration & conf,
                simulation::Simulation & sim,
                Pairlist const & pairlist_solute,
                Pairlist const & pairlist_solvent,
                Storage & storage,
                bool longrange, util::Algorithm_Timer & timer, bool master) {
  
  math::Periodicity<t_interaction_spec::boundary_type> periodicity(conf.current().box);
  Nonbonded_Innerloop<t_interaction_spec> innerloop(m_param);

  innerloop.init(sim,simulation::lj_func);

  const unsigned end_solute = topo.num_solute_atoms();
  for (unsigned i = 0; i < end_solute; ++i) {
    for (std::vector<unsigned int>::const_iterator
          j_it = pairlist_solute[i].begin()
        , j_to = pairlist_solute[i].end()
        ; j_it != j_to; ++j_it) {
      DEBUG(15, "\tQMMM LJ nonbonded_interaction: i " << i << " j " << *j_it);

      innerloop.lj_innerloop(topo, conf, i, *j_it, storage, periodicity);
    }
  }
  const unsigned end_solvent = topo.num_atoms();
  for (unsigned i = end_solute; i < end_solvent; ++i) {
    for (std::vector<unsigned int>::const_iterator
          j_it = pairlist_solvent[i].begin()
        , j_to = pairlist_solvent[i].end()
        ; j_it != j_to; ++j_it) {
      DEBUG(15, "\tQMMM LJ nonbonded_interaction: i " << i << " j " << *j_it);

      innerloop.lj_innerloop(topo, conf, i, *j_it, storage, periodicity);
    }
  }
}

/**
 * helper function to calculate the forces and energies from the
 * QMMM 1,4 interactions.
 */
void interaction::QMMM_Nonbonded_Outerloop
::one_four_outerloop(topology::Topology & topo,
        configuration::Configuration & conf,
        simulation::Simulation & sim,
        Storage & storage,
        int rank, int size) {
  SPLIT_MY_INNERLOOP(_one_four_outerloop, simulation::lj_func
                      , topo, conf, sim, storage, rank, size);
}

/**
 * helper function to calculate the forces and energies from the
 * QMMM 1,4 interactions.
 */
template<typename t_interaction_spec>
void interaction::QMMM_Nonbonded_Outerloop
::_one_four_outerloop(topology::Topology & topo,
        configuration::Configuration & conf,
        simulation::Simulation & sim,
        Storage & storage,
        int rank, int size) {
  DEBUG(7, "\tcalculating QMMM 1,4-interactions");

  math::Periodicity<t_interaction_spec::boundary_type> periodicity(conf.current().box);
  Nonbonded_Innerloop<t_interaction_spec> innerloop(m_param);
  innerloop.init(sim);

  topology::excl_cont_t::value_type::const_iterator it, to;
  unsigned int const num_solute_atoms = topo.num_solute_atoms();

  for (unsigned int i = rank; i < num_solute_atoms; i += size) {
    it = topo.qm_one_four_pair(i).begin();
    to = topo.qm_one_four_pair(i).end();
    for (; it != to; ++it) {
      DEBUG(10, "\tQMMM 1,4-interaction:" << i << " " << *it);
      innerloop.one_four_interaction_innerloop(topo, conf, i, *it, storage, periodicity);
    } // loop over 1,4 pairs
  } // loop over solute atoms
}

/**
 * helper function to calculate the forces and energies from the
 * Lennard-Jones exception interaction
 */
void interaction::QMMM_Nonbonded_Outerloop
::lj_exception_outerloop(topology::Topology & topo,
		     configuration::Configuration & conf,
		     simulation::Simulation & sim,
		     Storage & storage,
                     int rank, int size)
{
  SPLIT_MY_INNERLOOP(_lj_exception_outerloop, simulation::lj_func, topo, conf, sim, storage, rank, size);
}

/**
 * helper function to calculate the forces and energies from the
 * Lennard-Jones exception interaction
 */
template<typename t_interaction_spec>
void interaction::QMMM_Nonbonded_Outerloop
::_lj_exception_outerloop(topology::Topology & topo,
		     configuration::Configuration & conf,
		     simulation::Simulation & sim,
		     Storage & storage,
                     int rank, int size)
{
  DEBUG(7, "\tcalculate Lennard-Jones-exception-interactions");

  math::Periodicity<t_interaction_spec::boundary_type> periodicity(conf.current().box);
  Nonbonded_Innerloop<t_interaction_spec> innerloop(m_param);
  innerloop.init(sim);
  unsigned int const num_lj_exceptions = topo.qm_lj_exceptions().size();

  for (unsigned int i = rank; i < num_lj_exceptions; i += size) {
    const topology::lj_exception_struct & ljex = topo.qm_lj_exceptions()[i];

    DEBUG(15, "\tQMMM LJ exception: " << ljex.i << " " << ljex.j);
    innerloop.lj_exception_innerloop(topo, conf, ljex, storage, periodicity);
  } // loop over LJ exceptions
}

/**
 * calculate the interaction for a given atom pair.
 * SLOW! as it has to create the periodicity...
 */
int interaction::QMMM_Nonbonded_Outerloop::calculate_interaction
(
        topology::Topology & topo,
        configuration::Configuration & conf,
        simulation::Simulation & sim,
        unsigned int atom_i, unsigned int atom_j,
        math::Vec & force,
        double &e_lj, double &e_crf
        ) {
  SPLIT_MY_INNERLOOP(_calculate_interaction, simulation::lj_func
                      , topo, conf, sim, atom_i, atom_j, force, e_lj, e_crf);
  return 0;
}

template<typename t_interaction_spec>
int interaction::QMMM_Nonbonded_Outerloop
::_calculate_interaction(topology::Topology & topo,
        configuration::Configuration & conf,
        simulation::Simulation & sim,
        unsigned int atom_i, unsigned int atom_j,
        math::Vec & force,
        double & e_lj, double & e_crf) {
  math::Vec r;
  math::Periodicity<t_interaction_spec::boundary_type> periodicity(conf.current().box);

  Nonbonded_Term term;
  term.init(sim, t_interaction_spec::interaction_func);

  const lj_parameter_struct &lj =
          m_param.lj_parameter(topo.iac(atom_i),
          topo.iac(atom_j));

  periodicity.nearest_image(conf.current().pos(atom_i), conf.current().pos(atom_j), r);
  DEBUG(10, "\tni i " << conf.current().pos(atom_i)(0) << " / "
          << conf.current().pos(atom_i)(1) << " / "
          << conf.current().pos(atom_i)(2));
  DEBUG(10, "\tni j " << conf.current().pos(atom_j)(0) << " / "
          << conf.current().pos(atom_j)(1) << " / "
          << conf.current().pos(atom_j)(2));
  DEBUG(10, "\tni r " << r(0) << " / " << r(1) << " / " << r(2));
  double f = 0.0;
  term.lj_interaction(r, lj.c6, lj.c12,
          f, e_lj);
  force = f * r;

  return 0;
}

/**
 * calculate the hessian for a given atom.
 * this will be VERY SLOW !
 */
int interaction::QMMM_Nonbonded_Outerloop
::calculate_hessian(topology::Topology & topo,
        configuration::Configuration & conf,
        simulation::Simulation & sim,
        unsigned int atom_i, unsigned int atom_j,
        math::Matrix & hessian,
        PairlistContainer const & pairlist) {

  SPLIT_MY_INNERLOOP(_calculate_hessian, simulation::lj_func, topo, conf, sim, atom_i, atom_j, hessian, pairlist);
  return 0;
}

template<typename t_interaction_spec>
int interaction::QMMM_Nonbonded_Outerloop
::_calculate_hessian(topology::Topology & topo,
        configuration::Configuration & conf,
        simulation::Simulation & sim,
        unsigned int atom_i, unsigned int atom_j,
        math::Matrix & hessian,
        PairlistContainer const & pairlist) {

  hessian = 0.0;

  // loop over the pairlist

  //*************************
  // standard implementation
  //*************************

  // check whether the pair is in one of the shortrange pairlists
  bool calculate_pair =
          std::find(pairlist.solute_short[atom_i].begin(),
          pairlist.solute_short[atom_i].end(),
          atom_j) != pairlist.solute_short[atom_i].end() || // i-j in solute
          std::find(pairlist.solute_short[atom_j].begin(),
          pairlist.solute_short[atom_j].end(),
          atom_i) != pairlist.solute_short[atom_j].end() || // j-i in solute
          std::find(pairlist.solvent_short[atom_i].begin(),
          pairlist.solvent_short[atom_i].end(),
          atom_j) != pairlist.solvent_short[atom_i].end() || // i-j in solvent
          std::find(pairlist.solvent_short[atom_j].begin(),
          pairlist.solvent_short[atom_j].end(),
          atom_i) != pairlist.solvent_short[atom_j].end(); // j-i in solvent

  if (calculate_pair) {
    math::Vec r;
    math::Matrix h;
    math::Periodicity<t_interaction_spec::boundary_type> periodicity(conf.current().box);

    Nonbonded_Term term;
    term.init(sim, t_interaction_spec::interaction_func);

    DEBUG(8, "\thessian pair in pairlist: " << atom_i << " - " << atom_j);

    periodicity.nearest_image(conf.current().pos(atom_i)
                            , conf.current().pos(atom_j), r);
    const lj_parameter_struct &lj = m_param.lj_parameter(topo.iac(atom_i)
                                                       , topo.iac(atom_j));

    term.lj_crf_hessian(r, lj.c6, lj.c12, 0, h);

    for (unsigned int d1 = 0; d1 < 3; ++d1)
      for (unsigned int d2 = 0; d2 < 3; ++d2)
        hessian(d1, d2) += h(d1, d2);
  }

  return 0;
}
