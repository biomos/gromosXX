/**
 * @file parse_tcouple.cc
 * parse TCOUPLE into multibath
 */

#include <stdheader.h>

#include <topology/core/core.h>

#include <topology/solute.h>
#include <topology/solvent.h>
#include <topology/perturbed_atom.h>
#include <topology/perturbed_solute.h>

#include <topology/topology.h>

#include <simulation/multibath.h>
#include <simulation/parameter.h>

#include "parse_tcouple.h"

void util::parse_TCOUPLE(simulation::Parameter &param,
		   topology::Topology const & topo)
{
  // the last solute atom (and really the last index)
  int last_solute = topo.num_solute_atoms() - 1;
  // the last solvent atom
  int last_solvent = topo.num_atoms() - 1;

  // cases to handle
  // 0 0 0
  if (param.multibath.tcouple.ntt[0] == 0 && 
      param.multibath.tcouple.ntt[1] == 0 && 
      param.multibath.tcouple.ntt[2] == 0){
    // nothing
  }
  // 1 0 0
  else if (param.multibath.tcouple.ntt[0] == 1 && 
	   param.multibath.tcouple.ntt[1] == 0 && 
	   param.multibath.tcouple.ntt[2] == 0){

    param.multibath.couple = true;
    
    // the baths
    param.multibath.multibath.add_bath(param.multibath.tcouple.temp0[0], 
			     param.multibath.tcouple.tau[0]);
    param.multibath.multibath.add_bath(0, -1);
    // the atoms in the baths
    param.multibath.multibath.add_bath_index(last_solute, 0, 0, 1);
    // and an uncoupled one for the solvent...
    if (last_solvent != last_solute)
      param.multibath.multibath.add_bath_index(last_solvent, 0, 1, 1);

  }
  // 0 1 0
  else if (param.multibath.tcouple.ntt[0] == 0 && 
	   param.multibath.tcouple.ntt[1] == 1 && 
	   param.multibath.tcouple.ntt[2] == 0){
    param.multibath.couple = true;
    // the baths
    param.multibath.multibath.add_bath(param.multibath.tcouple.temp0[1], 
			     param.multibath.tcouple.tau[1]);
    param.multibath.multibath.add_bath(0, -1);
    // the atoms in the baths
    param.multibath.multibath.add_bath_index(last_solute, 0, 1, 0);
    // and an uncoupled one for the solvent...
    if (last_solvent != last_solute)
      param.multibath.multibath.add_bath_index(last_solvent, 0, 1, 1);

  }
  // 0 0 1
  else if (param.multibath.tcouple.ntt[0] == 0 && 
	   param.multibath.tcouple.ntt[1] == 0 && 
	   param.multibath.tcouple.ntt[2] == 1){
    param.multibath.couple = true;
    // the baths
    param.multibath.multibath.add_bath(param.multibath.tcouple.temp0[2], 
			     param.multibath.tcouple.tau[2]);
    param.multibath.multibath.add_bath(0, -1);
    // the atoms in the baths
    param.multibath.multibath.add_bath_index(last_solute, 0, 1, 1);
    param.multibath.multibath.add_bath_index(last_solvent, 0, 0, 0);
  }
  // 1 1 0
  else if (param.multibath.tcouple.ntt[0] == 1 && 
	   param.multibath.tcouple.ntt[1] == 1 && 
	   param.multibath.tcouple.ntt[2] == 0){
    param.multibath.couple = true;
    // the baths
    param.multibath.multibath.add_bath(param.multibath.tcouple.temp0[0], 
			     param.multibath.tcouple.tau[0]);
    param.multibath.multibath.add_bath(param.multibath.tcouple.temp0[1], 
			     param.multibath.tcouple.tau[1]);
    // the atoms in the baths
    param.multibath.multibath.add_bath_index(last_solute, 0, 0, 1);

    // and an uncoupled one for the solvent...
    if (last_solvent != last_solute){
      param.multibath.multibath.add_bath(0, -1);
      param.multibath.multibath.add_bath_index(last_solvent, 0, 2, 2);
    }
          
  }
  // 1 1 1
  else if (param.multibath.tcouple.ntt[0] == 1 && 
	   param.multibath.tcouple.ntt[1] == 1 && 
	   param.multibath.tcouple.ntt[2] == 1){
    param.multibath.couple = true;
    // the baths
    param.multibath.multibath.add_bath(param.multibath.tcouple.temp0[0], 
			     param.multibath.tcouple.tau[0]);
    param.multibath.multibath.add_bath(param.multibath.tcouple.temp0[1], 
			     param.multibath.tcouple.tau[1]);
    param.multibath.multibath.add_bath(param.multibath.tcouple.temp0[2], 
			     param.multibath.tcouple.tau[2]);
    // the atoms in the baths
    param.multibath.multibath.add_bath_index(last_solute, 0, 0, 1);
    param.multibath.multibath.add_bath_index(last_solvent, 0, 2, 2);
  }
  // 2 -2 0
  else if (param.multibath.tcouple.ntt[0] == 2 && 
	   param.multibath.tcouple.ntt[1] == -2 && 
	   param.multibath.tcouple.ntt[2] == 0){
    param.multibath.couple = true;
    // the bath
    param.multibath.multibath.add_bath(param.multibath.tcouple.temp0[0], 
			     param.multibath.tcouple.tau[0]);
    // the atoms in the bath
    param.multibath.multibath.add_bath_index(last_solute, 0, 0, 0);
    // and an uncoupled one for the solvent...
    if (last_solvent != last_solute){
      param.multibath.multibath.add_bath(0, -1);
      param.multibath.multibath.add_bath_index(last_solvent, 0, 1, 1);
    }

  }
  // 2 -2 1
  else if (param.multibath.tcouple.ntt[0] == 2 && 
	   param.multibath.tcouple.ntt[1] == -2 && 
	   param.multibath.tcouple.ntt[2] == 1){
    param.multibath.couple = true;
    // the baths
    param.multibath.multibath.add_bath(param.multibath.tcouple.temp0[0], 
			     param.multibath.tcouple.tau[0]);
    param.multibath.multibath.add_bath(param.multibath.tcouple.temp0[2], 
			     param.multibath.tcouple.tau[2]);
    // the atoms in the baths
    param.multibath.multibath.add_bath_index(last_solute, 0, 0, 0);
    param.multibath.multibath.add_bath_index(last_solvent, 0, 1, 1);
  }
  // 3 3 3
  else if (param.multibath.tcouple.ntt[0] == 3 && 
	   param.multibath.tcouple.ntt[1] == -3 && 
	   param.multibath.tcouple.ntt[2] == -3){
    param.multibath.couple = true;
    // the bath
    param.multibath.multibath.add_bath(param.multibath.tcouple.temp0[0], 
			     param.multibath.tcouple.tau[0]);
    // the atoms in the bath
    param.multibath.multibath.add_bath_index(last_solvent, 0, 0, 0);
  }
  // rest is not handled!
  else{
    io::messages.add("TCOUPLE param.multibath.tcouple.ntt combination not handled",
		     "InInput", io::message::error);
  }

}
