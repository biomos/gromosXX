/**
 * @file rottrans.cc
 * contains the template methods for
 * the class Rottrans_Constraints.
 */

#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../simulation/parameter.h"
#include "../../configuration/configuration.h"

#include "../../configuration/state_properties.h"
#include "../../math/periodicity.h"

#include "../../math/volume.h"

#include "../../algorithm/constraints/rottrans.h"

#include "../../util/template_split.h"
#include "../../util/debug.h"

#undef MODULE
#undef SUBMODULE
#define MODULE algorithm
#define SUBMODULE constraints

/**
 * Constructor.
 */
algorithm::Rottrans_Constraints
::Rottrans_Constraints(std::string const name)
  : Algorithm(name)
{
}

/**
 * Destructor.
 */
algorithm::Rottrans_Constraints
::~Rottrans_Constraints()
{
}

template<math::boundary_enum b, math::virial_enum do_virial>
static void _apply(topology::Topology & topo,
		  configuration::Configuration & conf,
		  simulation::Simulation & sim)
{
  // calculate c
  DEBUG(8, "Roto-translational constraints");
  
  math::Vec c_trans(0.0), c_rot(0.0), lambda_trans(0.0), lambda_rot(0.0);
  
  for(int i=0; i < sim.param().rottrans.last; ++i){

    const math::Vec diff = conf.current().pos(i) - conf.old().pos(i);
    
    c_trans(0) += topo.mass()(i) * diff(0);
    c_trans(1) += topo.mass()(i) * diff(1);
    c_trans(2) += topo.mass()(i) * diff(2);

    c_rot(0) += topo.mass()(i) * (conf.special().rottrans_constr.pos(i)(1) * diff(2) -
			      conf.special().rottrans_constr.pos(i)(2) * diff(1));

    c_rot(1) += topo.mass()(i) * (conf.special().rottrans_constr.pos(i)(2) * diff(0) -
			      conf.special().rottrans_constr.pos(i)(0) * diff(2));

    c_rot(2) += topo.mass()(i) * (conf.special().rottrans_constr.pos(i)(0) * diff(1) -
			      conf.special().rottrans_constr.pos(i)(1) * diff(0));
  }

  // lambda = - theta_inv * c
  for(int d1=0; d1<3; ++d1){
    lambda_trans(d1) -= 
      conf.special().rottrans_constr.theta_inv_trans(d1, d1) * c_trans(d1);
    DEBUG(10, "\tlambda[" << d1 << "] = " << lambda_trans(d1) << "\n");
  }
  
  for(int d1=0; d1 < 3; ++ d1){
    for(int d2=0; d2 < 3; ++d2){
      lambda_rot(d1) -=
	conf.special().rottrans_constr.theta_inv_rot(d1, d2) * c_rot(d2);
    }
    DEBUG(10, "\tlambda[" << d1+3 << "] = " << lambda_rot(d1) << "\n");
  }
  
  // update the positions
  for(int i=0; i < sim.param().rottrans.last; ++i){

    conf.current().pos(i)(0) +=
      ( lambda_trans(0) +
	lambda_rot(1) * conf.special().rottrans_constr.pos(i)(2) -
	lambda_rot(2) * conf.special().rottrans_constr.pos(i)(1)
	);

    conf.current().pos(i)(1) +=
      ( lambda_trans(1) -
	lambda_rot(0) * conf.special().rottrans_constr.pos(i)(2) +
	lambda_rot(2) * conf.special().rottrans_constr.pos(i)(0)
	);
    
    conf.current().pos(i)(2) +=
      ( lambda_trans(2) +
	lambda_rot(0) * conf.special().rottrans_constr.pos(i)(1) -
	lambda_rot(1) * conf.special().rottrans_constr.pos(i)(0)
	);
  }

  // update velocities
  for(unsigned int i=0; i<topo.num_atoms(); ++i)
    conf.current().vel(i) = (conf.current().pos(i) - conf.old().pos(i)) / 
      sim.time_step_size();

#ifndef NDEBUG
  //==================================================
  // CHECK
  //==================================================

  math::Vec v(0.0);
  for(int i=0; i < sim.param().rottrans.last; ++i){
    v += topo.mass()(i) * (math::cross(conf.special().rottrans_constr.pos(i),
				       conf.current().pos(i)));
  }
  DEBUG(10, "v: " << math::v2s(v));

#endif

}

/**
 * apply roto-translational constraints
 */
int algorithm::Rottrans_Constraints
::apply(topology::Topology & topo,
	configuration::Configuration & conf,
	simulation::Simulation & sim)
{
  DEBUG(7, "roto-constrational translaints");

  // the roto-constrational translaints
  if (sim.param().rottrans.rottrans){
    
    DEBUG(8, "\tapplying roto-constrational translaints");

    SPLIT_VIRIAL_BOUNDARY(_apply, topo, conf, sim);
    
  }
  
  return 0;		   
}


template<math::boundary_enum b, math::virial_enum do_virial>
static void _init(topology::Topology & topo,
		  configuration::Configuration & conf,
		  simulation::Simulation & sim,
		  std::ostream & os,
		  bool quiet)
{
  if (!quiet) {
    os << "ROTTRANS\n";
    os << "\tRoto-translational constraints\tON\n";
    os << "\tusing last atom:\t" << sim.param().rottrans.last << "\n\n";
    if (sim.param().start.read_rottrans)
      os << "\treading initial reference orientation and position form configuration";
    else
      os << "\treseting initial reference orientation and position.";
    os << "\nEND\n";
  }
  
  if (!sim.param().start.read_rottrans) {
    math::Periodicity<b> periodicity(conf.current().box);
    configuration::State_Properties sp(conf);

    math::Vec com_pos;
    math::Matrix com_e_kin;

    periodicity.put_chargegroups_into_box(conf, topo);
    sp.center_of_mass(0, sim.param().rottrans.last, topo.mass(), com_pos, com_e_kin);

    DEBUG(12, "com: " << math::v2s(com_pos));

    for (unsigned int i = 0; i < topo.num_atoms(); ++i)
      conf.current().pos(i) -= com_pos;
    periodicity.put_chargegroups_into_box(conf, topo);
    conf.exchange_state();
    for (unsigned int i = 0; i < topo.num_atoms(); ++i)
      conf.current().pos(i) -= com_pos;
    periodicity.put_chargegroups_into_box(conf, topo);
    conf.exchange_state();
    
    // if position restraints are used, the reference positions should be moved along
    if(sim.param().posrest.posrest != simulation::posrest_off){
      for(unsigned int i = 0; i < topo.num_atoms(); ++i){
	conf.special().reference_positions(i) -= com_pos;
	// and they should also be put in the box
	periodicity.put_into_box(conf.special().reference_positions(i));      
      }
    }
    
    // check
    // sp.center_of_mass(0, topo.num_solute_atoms(), topo.mass(), com_pos, com_e_kin);
    // DEBUG(10, "com (centered): " << math::v2s(com_pos));

    // store initial (reference) positions
    conf.special().rottrans_constr.pos.resize(sim.param().rottrans.last);

    for (int i = 0; i < sim.param().rottrans.last; ++i) {
      conf.special().rottrans_constr.pos(i) = conf.current().pos(i);
    }

    // now calculate the thetas
    std::vector<math::Vec> theta(6, math::Vec(0.0));
    const int X = 0;
    const int Y = 1;
    const int Z = 2;

    for (int i = 0; i < sim.param().rottrans.last; ++i) {

      theta[0](X) += topo.mass()(i);
      theta[1](Y) += topo.mass()(i);
      theta[2](Z) += topo.mass()(i);

      theta[3](X) += topo.mass()(i) * (conf.current().pos(i)(Z) * conf.current().pos(i)(Z) +
              conf.current().pos(i)(Y) * conf.current().pos(i)(Y));

      theta[3](Y) += topo.mass()(i) * (-conf.current().pos(i)(Y) * conf.current().pos(i)(X));

      theta[3](Z) += topo.mass()(i) * (-conf.current().pos(i)(X) * conf.current().pos(i)(Z));


      theta[4](X) += topo.mass()(i) * (-conf.current().pos(i)(X) * conf.current().pos(i)(Y));

      theta[4](Y) += topo.mass()(i) * (conf.current().pos(i)(X) * conf.current().pos(i)(X) +
              conf.current().pos(i)(Z) * conf.current().pos(i)(Z));

      theta[4](Z) += topo.mass()(i) * (-conf.current().pos(i)(Y) * conf.current().pos(i)(Z));


      theta[5](X) += topo.mass()(i) * (-conf.current().pos(i)(X) * conf.current().pos(i)(Z));


      theta[5](Y) += topo.mass()(i) * (-conf.current().pos(i)(Y) * conf.current().pos(i)(Z));

      theta[5](Z) += topo.mass()(i) * (conf.current().pos(i)(X) * conf.current().pos(i)(X) +
              conf.current().pos(i)(Y) * conf.current().pos(i)(Y));

    }

    DEBUG(10, "theta[4] " << math::v2s(theta[3]));
    DEBUG(10, "theta[5] " << math::v2s(theta[4]));
    DEBUG(10, "theta[6] " << math::v2s(theta[5]));

    const double d = math::dot(math::cross(theta[3], theta[4]), theta[5]);

    // inverse that
    math::Matrix & theta_inv_trans = conf.special().rottrans_constr.theta_inv_trans;
    math::Matrix & theta_inv_rot = conf.special().rottrans_constr.theta_inv_rot;

    theta_inv_trans = 0.0;
    theta_inv_rot = 0.0;

    theta_inv_trans(0, 0) = theta_inv_trans(1, 1) = theta_inv_trans(2, 2) =
            1.0 / theta[0](0);

    math::Vec h = math::cross(theta[4], theta[5]) / d;
    theta_inv_rot(0, 0) = h(0);
    theta_inv_rot(1, 0) = h(1);
    theta_inv_rot(2, 0) = h(2);

    h = math::cross(theta[3], theta[5]) / (-d);

    theta_inv_rot(0, 1) = h(0);
    theta_inv_rot(1, 1) = h(1);
    theta_inv_rot(2, 1) = h(2);

    h = math::cross(theta[3], theta[4]) / d;

    theta_inv_rot(0, 2) = h(0);
    theta_inv_rot(1, 2) = h(1);
    theta_inv_rot(2, 2) = h(2);
  }
}

int algorithm::Rottrans_Constraints
::init(topology::Topology & topo,
       configuration::Configuration & conf,
       simulation::Simulation & sim,
       std::ostream & os,
       bool quiet)
{

  if (sim.param().rottrans.rottrans){
    SPLIT_VIRIAL_BOUNDARY(_init, topo, conf, sim, os, quiet);
  }
  
  return 0;
}
