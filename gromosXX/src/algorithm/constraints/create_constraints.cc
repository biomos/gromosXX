/**
 * @file create_md_sequence.cc
 */

#include <util/stdheader.h>
#include <fstream>

#include <topology/core/core.h>
#include <topology/topology.h>
#include <simulation/multibath.h>
#include <simulation/parameter.h>
#include <simulation/simulation.h>
#include <configuration/energy.h>
#include <configuration/energy_average.h>
#include <configuration/configuration.h>
#include <interaction/interaction.h>
#include <interaction/interaction_types.h>

#include <io/argument.h>
#include <io/blockinput.h>
#include <io/instream.h>
#include <io/topology/in_topology.h>

#include <algorithm/algorithm.h>
#include <algorithm/algorithm_sequence.h>

#include <math/periodicity.h>
#include <algorithm/constraints/shake.h>
#include <algorithm/constraints/perturbed_shake.h>
#include <algorithm/constraints/lincs.h>
#include <algorithm/constraints/flexible_constraint.h>

#include <io/print_block.h>

#include "create_constraints.h"

#undef MODULE
#undef SUBMODULE
#define MODULE algorithm
#define SUBMODULE constraints


template<math::virial_enum do_virial>
static int _create_constraints(algorithm::Algorithm_Sequence &md_seq,
			       topology::Topology &topo,
			       configuration::Configuration & conf,
			       simulation::Simulation & sim,
			       io::In_Topology &it)
{
  // CONSTRAINTS
  DEBUG(7, "Constrain solute?");

  if (sim.param().start.shake_pos || sim.param().start.shake_vel){
    if (sim.param().constraint.solute.algorithm == simulation::constr_off)
      io::messages.add("no shaking of initial solute positions / velocities "
		       "without using constraints during the simulation.",
		       "create md sequence",
		       io::message::warning);

    if (sim.param().constraint.solvent.algorithm == simulation::constr_off)
      io::messages.add("no shaking of initial solvent positions / velocities "
		       "without using constraints during the simulation.",
		       "create md sequence",
		       io::message::warning);
  }

  switch(sim.param().constraint.solute.algorithm){
    case simulation::constr_shake:
      {
	// SHAKE
	algorithm::Shake<do_virial> * s = 
	  new algorithm::Shake<do_virial>
	  (sim.param().constraint.solute.shake_tolerance);
	it.read_harmonic_bonds(s->parameter());
	s->init(topo, conf, sim);
	md_seq.push_back(s);
	
	if (sim.param().perturbation.perturbation){
	  algorithm::Perturbed_Shake<do_virial> * ps =
	    new algorithm::Perturbed_Shake<do_virial>(*s);
	  ps->init(topo, conf, sim);
	  md_seq.push_back(ps);
	}
      }
    case simulation::constr_lincs:
      {
	algorithm::Lincs<do_virial> * s =
	  new algorithm::Lincs<do_virial>;
	it.read_harmonic_bonds(s->parameter());
	s->init(topo, conf, sim);
	md_seq.push_back(s);

	if (sim.param().perturbation.perturbation){
	  io::messages.add("no free energy derivatives for LINCS, so you better don't "
			   "change constrained bond lengths", "create_constraints",
			   io::message::warning);
	}
	break;
      }
    case simulation::constr_flexshake:
      {
	algorithm::Flexible_Constraint<do_virial> * fs = 
	  new algorithm::Flexible_Constraint<do_virial>
	  (sim.param().constraint.solute.shake_tolerance);
	it.read_harmonic_bonds(fs->parameter());
	fs->init(topo, conf, sim);
	md_seq.push_back(fs);

	if (sim.param().perturbation.perturbation){
	  io::messages.add("no free energy derivatives for flexible Constraints, "
			   "so you better don't change constrained bond lengths", 
			   "create_constraints",
			   io::message::warning);
	}
      }
    default:
      {
	// no constraints...
      }
      
  }

  // sovlent (if not the same as solute)
  if (sim.param().constraint.solute.algorithm == 
      sim.param().constraint.solvent.algorithm) return 0;

  switch(sim.param().constraint.solvent.algorithm){
    case simulation::constr_shake:
      {
	// SHAKE
	algorithm::Shake<do_virial> * s = 
	  new algorithm::Shake<do_virial>
	  (sim.param().constraint.solvent.shake_tolerance);
	it.read_harmonic_bonds(s->parameter());
	s->init(topo, conf, sim);
	md_seq.push_back(s);
	
	break;
      }
    case simulation::constr_lincs:
      {
	algorithm::Lincs<do_virial> * s =
	  new algorithm::Lincs<do_virial>;
	it.read_harmonic_bonds(s->parameter());
	s->init(topo, conf, sim);
	md_seq.push_back(s);
	
	break;
      }
    case simulation::constr_flexshake:
      {
	io::messages.add("Flexible Shake not implemented for solvent",
			 "create_constraints", io::message::error);
	break;
      }
    default:
      {
	// no constraints
	// should already be warned from In_Parameter
      }
  }

  return 0;
}


int algorithm::create_constraints(algorithm::Algorithm_Sequence & md_seq,
				  topology::Topology & topo,
				  configuration::Configuration & conf,
				  simulation::Simulation & sim,
				  io::In_Topology &it)
{

  DEBUG(7, "solute:  " << sim.param().constraint.solute.algorithm);
  DEBUG(7, "solvent: " << sim.param().constraint.solvent.algorithm);
  DEBUG(7, "\tNTC: " << sim.param().constraint.ntc);
  
  switch(sim.param().pcouple.virial){
    case math::no_virial:
    case math::molecular_virial:
      {
	DEBUG(8, "\twith no virial");
	
	return _create_constraints<math::no_virial>(md_seq, topo, conf, sim, it);
      }
    case math::atomic_virial:
      {
	DEBUG(8, "\twith atomic virial");

	if (sim.param().constraint.solute.algorithm == simulation::constr_lincs){
	  io::messages.add("atomic virial not implemented for lincs",
			   "create_constraints",
			   io::message::error);
	}
	
	return _create_constraints<math::atomic_virial>(md_seq, topo, conf, sim, it);
      }
    default:
      io::messages.add("wrong virial type", "create_constraints",
		       io::message::error);
      return -1;
  }
  
}
