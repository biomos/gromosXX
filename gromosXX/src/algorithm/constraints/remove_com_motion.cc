/**
 * @file remove_com_motion.cc
 * remove com motion.
 */

#include <stdheader.h>

#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <simulation/simulation.h>
#include <configuration/configuration.h>

#include "remove_com_motion.h"

#include <io/print_block.h>

#undef MODULE
#undef SUBMODULE
#define MODULE algorithm
#define SUBMODULE constraints

/**
 * apply the COM removal.
 */
int algorithm::Remove_COM_Motion
::apply(topology::Topology & topo,
	configuration::Configuration & conf,
	simulation::Simulation & sim)
{
  bool remove_it = false;
  bool print_it = false;
  
  const double start = util::now();

  // check if nothing to do
  if (sim.steps() == 0){
    if (sim.param().start.remove_com) remove_it = true;
    print_it = true;
  }
  else{
    if (sim.param().centreofmass.skip_step &&
	(sim.steps() % sim.param().centreofmass.skip_step) == 0)
      remove_it = true;
    if (sim.param().print.centreofmass &&
	(sim.steps() % sim.param().print.centreofmass) == 0)
      print_it = true;
  }

  if (!print_it && !remove_it) return 0;

  math::Vec com_v (0.0), com_r(0.0);
  double com_mass = 0.0;
  
  for(unsigned int i = 0; i < topo.num_atoms(); ++i){

    com_mass += topo.mass()(i);
    com_v += topo.mass()(i) * conf.current().vel(i);
    // positions should be at same time than velocities
    com_r += topo.mass()(i) * conf.current().pos(i) - 
      0.5 * topo.mass()(i) * conf.current().vel(i) * sim.time_step_size();
  }

  com_v /= com_mass;
  com_r /= com_mass;

  double ekin_trans = 0.5*com_mass*abs2(com_v);

  DEBUG(7, "totmass " << com_mass);
  DEBUG(7, "com_v " << math::v2s(com_v));
  DEBUG(7, "com_r " << math::v2s(com_r));
  DEBUG(7, "com_Ekin " << ekin_trans);

  math::Vec com_L(0.0);
  math::Matrix com_I;
  com_I = 0.0;
  
  for(unsigned int i = 0; i < topo.num_atoms(); ++i){
    math::Vec r = conf.current().pos(i) - 
      0.5 * sim.time_step_size() * conf.current().vel(i) - com_r;

    DEBUG(15, "pos  " << i << " " 
	  << math::v2s(conf.current().pos(i) - 
	  0.5 * sim.time_step_size() * conf.current().vel(i)));
    DEBUG(15, "posp " << i << " " << math::v2s(r));
    
    // should this be pos or r???
    com_L += topo.mass()(i) * 
      math::cross(conf.current().pos(i), conf.current().vel(i));

    // inertia tensor
    // double r2 = abs2(r);
    com_I(0,0) += topo.mass()(i) * (r(1)*r(1)+r(2)*r(2));
    com_I(1,1) += topo.mass()(i) * (r(0)*r(0)+r(2)*r(2));
    com_I(2,2) += topo.mass()(i) * (r(0)*r(0)+r(1)*r(1));
    com_I(1,0) += topo.mass()(i) * (-r(0)*r(1));
    com_I(0,1) += topo.mass()(i) * (-r(0)*r(1));
    com_I(2,0) += topo.mass()(i) * (-r(0)*r(2));
    com_I(0,2) += topo.mass()(i) * (-r(0)*r(2));
    com_I(2,1) += topo.mass()(i) * (-r(1)*r(2));
    com_I(1,2) += topo.mass()(i) * (-r(1)*r(2));
  }

  com_L -= com_mass * math::cross(com_r, com_v);
  
  DEBUG(7, "Angular momentum " << math::v2s(com_L));
  
  // invert the inertia tensor
  math::Matrix com_II;
  const double denom = -com_I(2,0)*com_I(2,0)*com_I(1,1)
    + 2 * com_I(0,1) * com_I(0,2) * com_I(1,2)
    - com_I(0, 0) * com_I(1,2) * com_I(1,2)
    - com_I(0,1) * com_I(0,1) * com_I(2,2)
    + com_I(0,0) * com_I(1,1) * com_I(2,2);
  
  com_II(0,0) = (-com_I(1,2)*com_I(1,2) + com_I(1,1) * com_I(2,2));
  com_II(1,0) = com_II(0,1) = (com_I(0,2) * com_I(1,2)
    - com_I(0,1) * com_I(2,2));
  com_II(0,2) = com_II(2,0) = (-com_I(0,2)*com_I(1,1)
    + com_I(0,1)*com_I(1,2));

  com_II(1,1) = (-com_I(0,2)*com_I(0,2) + com_I(0,0) * com_I(2,2));
  com_II(1,2) = com_II(2,1) = (com_I(0,1)*com_I(0,2)
    - com_I(0,0) * com_I(1,2));

  com_II(2,2) = (-com_I(0,1)*com_I(0,1) + com_I(0,0)*com_I(1,1));
  
  DEBUG(7, "inertia tensor:\n"<< math::m2s(com_I));
  DEBUG(7, "determinant : " << denom);
  DEBUG(7, "inverted tens :\n" << math::m2s(com_II));
  
  // get the angular velocity around the COM
  math::Vec com_O;
  if (denom < math::epsilon)
    com_O = math::Vec(0.0, 0.0, 0.0);
  else
    com_O = math::product(com_II, com_L) / denom;
  
  DEBUG(7, " angular velocity " << math::v2s(com_O));

  double ekin_rot =  0.5 * dot(com_O, com_L);
  DEBUG(7, " com_Ekin_rot " << ekin_rot);
   
  if (print_it){
    io::print_CENTREOFMASS(os, ekin_trans, ekin_rot);
  }

  // remove if necessary
  // angular momentum
  if(remove_it && sim.param().centreofmass.remove_rot){

    os << "removing center of mass angular momentum\n";
    
    for(unsigned int i=0; i<topo.num_atoms(); ++i){
      math::Vec r = conf.current().pos(i) - 
	0.5 * sim.time_step_size() * conf.current().vel(i) - com_r;
      conf.current().vel(i) -= math::cross(com_O, r); 
    }
  }

  // translational momentum
  if(remove_it && sim.param().centreofmass.remove_trans){
    
    os << "removing center of mass translational momentum\n";
    
    // get the corrected velocities
    for(unsigned int i=0; i<topo.num_atoms(); ++i){
      conf.current().vel(i) -= com_v;
    }
  }

  m_timing += util::now() - start;

  // return success!
  return 0;
		   
}

