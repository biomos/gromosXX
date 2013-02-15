/**
 * @file remove_com_motion.cc
 * remove com motion.
 */

#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"

#include "remove_com_motion.h"

#include "../../io/print_block.h"

#undef MODULE
#undef SUBMODULE
#define MODULE algorithm
#define SUBMODULE constraints

int algorithm::Remove_COM_Motion::init
(
 topology::Topology &topo, 
 configuration::Configuration &conf,
 simulation::Simulation &sim,
 std::ostream &os,
 bool quiet
)
{
  if (quiet) return 0;
  
  os << "CENTRE OF MASS MOTION\n";

  if (sim.param().centreofmass.skip_step){
    if (sim.param().centreofmass.skip_step > 1)
      os << "\tremoving centre of mass motion every " 
	 << sim.param().centreofmass.skip_step
	 << " steps\n";
    else
      os << "\tremoving centre of mass motion every step\n";

    if (sim.param().centreofmass.remove_rot){
      os << "\tremoving centre of mass rotation" << std::endl;
    }
    if (sim.param().centreofmass.remove_trans){
      os << "\tremoving centre of mass translation" << std::endl;
    }
    os << "\n";
  }
  if (sim.param().print.centreofmass > 1){
    os << "\tprinting centre of mass motion every "
       << sim.param().print.centreofmass
       << " steps\n";
  }
  if (sim.param().print.centreofmass == 1){
    os << "\tprinting centre of mass motion every step\n";
  }

  if (sim.param().start.remove_com_translation)
    os << "\n\tremoving initial centre of mass translation\n";
  if (sim.param().start.remove_com_rotation)
    os << "\n\tremoving initial centre of mass rotation\n";
  
  os << "END\n";
  
  return 0;
};

double algorithm::Remove_COM_Motion
::remove_com_translation
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation & sim,
 bool remove_trans
 )
{
  math::Vec com_v (0.0);
  double com_mass = 0.0;
  
  for(unsigned int i = 0; i < topo.num_atoms(); ++i){
    com_mass += topo.mass()(i);
    com_v += topo.mass()(i) * conf.current().vel(i);
  }

  com_v /= com_mass;
  double ekin_trans = 0.5*com_mass*abs2(com_v);

  // remove if necessary
  if(remove_trans){

    // os << "removing center of mass translational momentum\n";
    for(unsigned int i=0; i<topo.num_atoms(); ++i){
      conf.current().vel(i) -= com_v;
    }
  }
  return ekin_trans;
}

double algorithm::Remove_COM_Motion
::remove_com_rotation
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation & sim,
 bool remove_rot
 )
{
  math::Vec com_v (0.0);
  math::Vec com_r(0.0);
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

  DEBUG(7, "com_r " << math::v2s(com_r));
  DEBUG(7, "com_v " << math::v2s(com_v));
  
  math::Vec com_L(0.0);
  math::Matrix com_I;
  com_I = 0.0;
  
  for(unsigned int i = 0; i < topo.num_atoms(); ++i){
    math::Vec r = conf.current().pos(i) - 
      0.5 * sim.time_step_size() * conf.current().vel(i) - com_r;

    com_L += topo.mass()(i) * 
      math::cross(r, (conf.current().vel(i) - com_v));

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
   
  // remove if necessary
  if(remove_rot){
    // os << "removing center of mass angular momentum\n";
    
    for(unsigned int i=0; i<topo.num_atoms(); ++i){
      math::Vec r = conf.current().pos(i) - 
	0.5 * sim.time_step_size() * conf.current().vel(i) - com_r;
      conf.current().vel(i) -= math::cross(com_O, r); 
    }
  }

  return ekin_rot;
}


/**
 * apply the COM removal.
 */
int algorithm::Remove_COM_Motion
::apply(topology::Topology & topo,
	configuration::Configuration & conf,
	simulation::Simulation & sim)
{
  bool remove_rot = false;
  bool remove_trans = false;
  bool print_it = false;
  
  m_timer.start();

  // check if nothing to do
  if (sim.steps() == 0){
    remove_rot = sim.param().start.remove_com_rotation;
    remove_trans = sim.param().start.remove_com_translation;
    if (sim.param().print.centreofmass) print_it = true;
  }
  else{
    if (sim.param().centreofmass.skip_step &&
	(sim.steps() % sim.param().centreofmass.skip_step) == 0)
      remove_rot = remove_trans = true;
    if (sim.param().print.centreofmass &&
	(sim.steps() % sim.param().print.centreofmass) == 0)
      print_it = true;
  }

  DEBUG(9, "centre of mass: print " << print_it << " remove " <<
           (remove_trans || remove_rot) );
  if (!print_it && !remove_trans && !remove_rot) return 0;
  
  if (sim.steps() != 0){
    remove_rot = remove_rot && sim.param().centreofmass.remove_rot;
    remove_trans = remove_trans && sim.param().centreofmass.remove_trans;
  }
  
  DEBUG(9, "centre of mass: trans " << remove_trans << " rot " << remove_rot);
  
  double ekin_trans = 0.0, ekin_rot = 0.0;

  if (print_it || remove_trans){
    ekin_trans = remove_com_translation(topo, conf, sim, remove_trans);
  }
  if (print_it || remove_rot){
    ekin_rot = remove_com_rotation(topo, conf, sim, remove_rot);
  }

  if (print_it){
    io::print_CENTREOFMASS(os, ekin_trans, ekin_rot);
  }

  m_timer.stop();

  return 0;		   
}

double algorithm::Remove_COM_Motion
::add_com_rotation
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation & sim,
 math::Vec com_L
 )
{
  // assumption is that com_L of the velocities is zero right now

  math::Vec com_v (0.0);
  math::Vec com_r(0.0);
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

  math::Matrix com_I;
  com_I = 0.0;
  
  for(unsigned int i = 0; i < topo.num_atoms(); ++i){
    math::Vec r = conf.current().pos(i) - 
      0.5 * sim.time_step_size() * conf.current().vel(i) - com_r;

    // inertia tensor
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

  // com_L -= com_mass * math::cross(com_r, com_v);
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
   
  // and add it
    
  for(unsigned int i=0; i<topo.num_atoms(); ++i){
    math::Vec r = conf.current().pos(i) - 
      0.5 * sim.time_step_size() * conf.current().vel(i) - com_r;
    conf.current().vel(i) += math::cross(com_O, r); 
  }

  return ekin_rot;
}
