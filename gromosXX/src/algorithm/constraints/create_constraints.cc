/**
 * @file create_constraints.cc
 * create the constraint algorithm
 * choosing from SHAKE, LINCS
 * and enabling or disabling perturbation
 */

#include <stdheader.h>
#include <fstream>

#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <simulation/simulation.h>
#include <configuration/configuration.h>

#include <configuration/state_properties.h>

#include <algorithm/algorithm/algorithm_sequence.h>

#include <interaction/interaction.h>
#include <interaction/interaction_types.h>

#include <interaction/forcefield/forcefield.h>

#include <math/periodicity.h>
#include <math/volume.h>

#include <io/argument.h>
#include <io/blockinput.h>
#include <io/instream.h>
#include <io/topology/in_topology.h>

#include <algorithm/constraints/shake.h>
#include <algorithm/constraints/perturbed_shake.h>
#include <algorithm/constraints/lincs.h>
#include <algorithm/constraints/flexible_constraint.h>
#include <algorithm/constraints/perturbed_flexible_constraint.h>

#include <algorithm/constraints/rottrans.h>

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
			       io::In_Topology &it,
			       bool quiet)
			       
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
	if (!sim.param().perturbation.perturbation){
	  
	  // SHAKE
	  algorithm::Shake<do_virial> * s = 
	    new algorithm::Shake<do_virial>
	    (sim.param().constraint.solute.shake_tolerance);
	  it.read_harmonic_bonds(s->parameter());
	  // s->init(topo, conf, sim, quiet);
	  md_seq.push_back(s);
	  
	}
	else{
	  // perturbed shake also calls normal shake...
	  algorithm::Perturbed_Shake<do_virial> * ps =
	    new algorithm::Perturbed_Shake<do_virial>
	    (sim.param().constraint.solute.shake_tolerance);
	  it.read_harmonic_bonds(ps->parameter());
	  // ps->init(topo, conf, sim, quiet);
	  md_seq.push_back(ps);

	}
	break;
      }
    case simulation::constr_lincs:
      {
	algorithm::Lincs<do_virial> * s =
	  new algorithm::Lincs<do_virial>;
	it.read_harmonic_bonds(s->parameter());
	// s->init(topo, conf, sim, quiet);
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
	// let's try to get the forcefield
	interaction::Forcefield * ff = NULL;

	for(size_t i=0; i < md_seq.size(); ++i){
	  if (md_seq[i]->name == "Forcefield"){
	    DEBUG(8, "flexible shake: forcefield found");
	    ff = dynamic_cast<interaction::Forcefield *>(md_seq[i]);
	    break;
	  }
	}

	if (!ff){
	  io::messages.add("no forcefield found", "create_constraints",
			   io::message::error);
	}

	algorithm::Flexible_Constraint<do_virial> * fs = 
	  new algorithm::Flexible_Constraint<do_virial>
	  (sim.param().constraint.solute.shake_tolerance, 1000, ff);

	it.read_harmonic_bonds(fs->parameter());
	// fs->init(topo, conf, sim, quiet);
	md_seq.push_back(fs);

	if (sim.param().perturbation.perturbation){

	  algorithm::Perturbed_Flexible_Constraint<do_virial> * pfc =
	    new algorithm::Perturbed_Flexible_Constraint<do_virial>(*fs);
	  // pfc->init(topo, conf, sim, quiet);
	  md_seq.push_back(pfc);

	}

	/*
	io::messages.add("flexible constraints disabled",
			 "create_constraints",
			 io::message::error);
	*/

	break;
      }
    default:
      {
	// no constraints...
      }
      
  }

  // sovlent (if not the same as solute)
  if (sim.param().constraint.solute.algorithm != 
      sim.param().constraint.solvent.algorithm){

    switch(sim.param().constraint.solvent.algorithm){
      case simulation::constr_shake:
	{
	  // SHAKE
	  algorithm::Shake<do_virial> * s = 
	    new algorithm::Shake<do_virial>
	    (sim.param().constraint.solvent.shake_tolerance);
	  it.read_harmonic_bonds(s->parameter());
	  // s->init(topo, conf, sim, quiet);
	  md_seq.push_back(s);
	  
	  break;
	}
      case simulation::constr_lincs:
	{
	  algorithm::Lincs<do_virial> * s =
	    new algorithm::Lincs<do_virial>;
	  it.read_harmonic_bonds(s->parameter());
	  // s->init(topo, conf, sim, quiet);
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
  }
  
  // roto-translational constraints
  if (sim.param().rottrans.rottrans){
    DEBUG(8, "creating roto-translational constraints");
    
    algorithm::Rottrans_Constraints<do_virial> * rtc =
      new algorithm::Rottrans_Constraints<do_virial>();
    md_seq.push_back(rtc);
  }

  return 0;
}


int algorithm::create_constraints(algorithm::Algorithm_Sequence & md_seq,
				  topology::Topology & topo,
				  configuration::Configuration & conf,
				  simulation::Simulation & sim,
				  io::In_Topology &it,
				  bool quiet)
{

  DEBUG(7, "solute:  " << sim.param().constraint.solute.algorithm);
  DEBUG(7, "solvent: " << sim.param().constraint.solvent.algorithm);
  DEBUG(7, "\tNTC: " << sim.param().constraint.ntc);
  DEBUG(7, "rottrans: " << sim.param().rottrans.rottrans);
  
  switch(sim.param().pcouple.virial){
    case math::no_virial:
    case math::molecular_virial:
      {
	DEBUG(8, "\twith no virial");
	
	return _create_constraints<math::no_virial>(md_seq, topo, conf, sim, it, quiet);
      }
    case math::atomic_virial:
      {
	DEBUG(8, "\twith atomic virial");

	if (sim.param().constraint.solute.algorithm == simulation::constr_lincs){
	  io::messages.add("atomic virial not implemented for lincs",
			   "create_constraints",
			   io::message::error);
	}
	
	return _create_constraints<math::atomic_virial>(md_seq, topo, conf, sim, it, quiet);
      }
    default:
      io::messages.add("wrong virial type", "create_constraints",
		       io::message::error);
      return -1;
  }
  
}
