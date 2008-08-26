/**
 * @file out_configuration.cc
 * definition of the Out_Configuration methods.
 */

#include <stdheader.h>

#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <simulation/simulation.h>
#include <simulation/parameter.h>
#include <configuration/configuration.h>
#include <configuration/energy.h>

#include <math/periodicity.h>
#include <math/volume.h>

#include <io/print_block.h>
#include <io/argument.h>

#include <sstream>

#include "out_configuration.h"

#include <util/replica_data.h>
#include <util/template_split.h>
#include <util/debug.h>

#include <limits>

#undef MODULE
#undef SUBMODULE
#define MODULE io
#define SUBMODULE configuration

// declarations
static void _print_energyred_helper(std::ostream & os, configuration::Energy const &e);

static void _print_volumepressurered_helper(std::ostream &os,
					    double mass,
					    simulation::Multibath const & m,
					    std::vector<double> const & s,
					    configuration::Energy const & e,
					    math::Box const & b,
					    math::boundary_enum t,
					    math::Matrix const & p,
					    math::Matrix const & v,
					    math::Matrix const & k);


io::Out_Configuration::Out_Configuration(std::string title,
					 std::ostream & os)
  : m_pos_traj(NULL),
    m_final_conf(NULL),
    m_vel_traj(NULL),
    m_force_traj(NULL),
    m_energy_traj(NULL),
    m_has_replica_traj(false),
    m_replica_traj(NULL),
    m_free_energy_traj(NULL),
    m_blockaveraged_energy(NULL),
    m_blockaveraged_free_energy(NULL),
    m_ramd_traj(NULL),
    m_special_traj(NULL),
    m_output(os),
    m_final(false),
    m_replica(false),
    m_every_pos(0),
    m_every_vel(0),
    m_every_force(0),
    m_every_energy(0),
    m_every_free_energy(0),
    m_every_blockaverage(0),
    m_every_ramd(0),
    m_every_cos_pos(0),
    m_write_blockaverage_energy(false),
    m_write_blockaverage_free_energy(false),
    m_precision(9),
    m_force_precision(9),
    m_distance_restraint_precision(7),
    m_width(15),
    m_force_width(18),
    m_title(title),
    m_compressed(false),
    minimum_energy(std::numeric_limits<double>::max())
{
  _print_title(m_title, "output file", os);
}

io::Out_Configuration::~Out_Configuration()
{
  // std::cout << "out_configuration destructor" << std::endl;

  if (m_every_pos) {
    if (m_compressed) {
      dynamic_cast<io::ogzstream*> (m_pos_traj)->flush();
      dynamic_cast<io::ogzstream*> (m_pos_traj)->close();
    } else {
      dynamic_cast<std::ofstream*> (m_pos_traj)->flush();
      dynamic_cast<std::ofstream*> (m_pos_traj)->close();
    }
    delete m_pos_traj;
  }
  if (m_final){
    if (m_compressed) {
      dynamic_cast<io::ogzstream*> (m_final_conf)->flush();
      dynamic_cast<io::ogzstream*> (m_final_conf)->close();
    } else {
      dynamic_cast<std::ofstream*> (m_final_conf)->flush();
      dynamic_cast<std::ofstream*> (m_final_conf)->close();
    }
    delete m_final_conf;
  }

  if (m_every_vel){
    if (m_compressed) {
      dynamic_cast<io::ogzstream*> (m_vel_traj)->flush();
      dynamic_cast<io::ogzstream*> (m_vel_traj)->close();
    } else {
      dynamic_cast<std::ofstream*> (m_vel_traj)->flush();
      dynamic_cast<std::ofstream*> (m_vel_traj)->close();
    }
    delete m_vel_traj;
  }
    
  if (m_every_force){
    if (m_compressed) {
      dynamic_cast<io::ogzstream*> (m_force_traj)->flush();
      dynamic_cast<io::ogzstream*> (m_force_traj)->close();
    } else {
      dynamic_cast<std::ofstream*> (m_force_traj)->flush();
      dynamic_cast<std::ofstream*> (m_force_traj)->close();
    }
    delete m_force_traj;
  }
    
  if (m_every_energy){
    if (m_compressed) {
      dynamic_cast<io::ogzstream*> (m_energy_traj)->flush();
      dynamic_cast<io::ogzstream*> (m_energy_traj)->close();
    } else {
      dynamic_cast<std::ofstream*> (m_energy_traj)->flush();
      dynamic_cast<std::ofstream*> (m_energy_traj)->close();
    }
    delete m_energy_traj;
  }

  if (m_has_replica_traj){
    if (m_compressed) {
      dynamic_cast<io::ogzstream*> (m_replica_traj)->flush();
      dynamic_cast<io::ogzstream*> (m_replica_traj)->close();
    } else {
      dynamic_cast<std::ofstream*> (m_replica_traj)->flush();
      dynamic_cast<std::ofstream*> (m_replica_traj)->close();
    }
    delete m_replica_traj;
  }

  if (m_every_free_energy){
    if (m_compressed) {
      dynamic_cast<io::ogzstream*> (m_free_energy_traj)->flush();
      dynamic_cast<io::ogzstream*> (m_free_energy_traj)->close();
    } else {
      dynamic_cast<std::ofstream*> (m_free_energy_traj)->flush();
      dynamic_cast<std::ofstream*> (m_free_energy_traj)->close();
    }
    delete m_free_energy_traj;
  }

  if (m_write_blockaverage_energy){
    if (m_compressed) {
      dynamic_cast<io::ogzstream*> (m_blockaveraged_energy)->flush();
      dynamic_cast<io::ogzstream*> (m_blockaveraged_energy)->close();
    } else {
      dynamic_cast<std::ofstream*> (m_blockaveraged_energy)->flush();
      dynamic_cast<std::ofstream*> (m_blockaveraged_energy)->close();
    }
    delete m_blockaveraged_energy;
  }
  
  if (m_write_blockaverage_free_energy){
    if (m_compressed) {
      dynamic_cast<io::ogzstream*> (m_blockaveraged_free_energy)->flush();
      dynamic_cast<io::ogzstream*> (m_blockaveraged_free_energy)->close();
    } else {
      dynamic_cast<std::ofstream*> (m_blockaveraged_free_energy)->flush();
      dynamic_cast<std::ofstream*> (m_blockaveraged_free_energy)->close();
    }
    delete m_blockaveraged_free_energy;
  }
  
  if (m_every_ramd){
    if (m_compressed) {
      dynamic_cast<io::ogzstream*> (m_ramd_traj)->flush();
      dynamic_cast<io::ogzstream*> (m_ramd_traj)->close();
    } else {
      dynamic_cast<std::ofstream*> (m_ramd_traj)->flush();
      dynamic_cast<std::ofstream*> (m_ramd_traj)->close();
    }
    delete m_ramd_traj;
  }  
  
  if (m_every_cos_pos){ // add others if there are any
    if (m_compressed) {
      dynamic_cast<io::ogzstream*> (m_special_traj)->flush();
      dynamic_cast<io::ogzstream*> (m_special_traj)->close();
    } else {
      dynamic_cast<std::ofstream*> (m_special_traj)->flush();
      dynamic_cast<std::ofstream*> (m_special_traj)->close();
    }
    delete m_special_traj;
  }  
}

void io::Out_Configuration::_print_title(std::string title,
					 std::string name,
					 std::ostream &os)
{
  os << "TITLE\n\t"
     << title << "\n"
     << "\t" << name 
     << "\nEND\n";
}

void io::Out_Configuration::init(io::Argument & args,
				 simulation::Parameter const & param)
{
   // has to be the first so that all streams get initialised correctly
  if (args.count(argname_gzip) >= 0)
    compressed(true);  
  
  if (args.count(argname_fin) > 0)
    final_configuration(args[argname_fin]);
  else io::messages.add("argument "+argname_fin+" for final configuration required!",
			"Out_Configuration",
			io::message::error);

  if (args.count(argname_trj) > 0)
    trajectory(args[argname_trj], param.write.position);
  else if (param.write.position)
    io::messages.add("write trajectory but no "+argname_trj+" argument",
		     "Out_Configuration",
		     io::message::error);

  if (args.count(argname_trv) > 0)
    velocity_trajectory(args[argname_trv], param.write.velocity);
  else if (param.write.velocity)
    io::messages.add("write velocity trajectory but no trv argument",
		     "Out_Configuration",
		     io::message::error);

  if (args.count(argname_trf) > 0)
    force_trajectory(args[argname_trf], param.write.force);
  else if (param.write.force)
    io::messages.add("write force trajectory but no trf argument",
		     "Out_Configuration",
		     io::message::error);
  
  if (args.count(argname_trs) > 0)
    special_trajectory(args[argname_trs], param.polarize.write);
  else if (param.polarize.write) // check for other that also go to this traj.
    io::messages.add("write special trajectory but no trs argument",
		     "Out_Configuration",
		     io::message::error);

  if (args.count(argname_re) > 0)
    replica_trajectory(args[argname_re]);
  
  if (args.count(argname_tre) > 0)
    energy_trajectory(args[argname_tre], param.write.energy);
  else if (param.write.energy)
    io::messages.add("write energy trajectory but no "+argname_tre+" argument",
		     "Out_Configuration",
		     io::message::error);

  if (args.count(argname_trg) > 0)
    free_energy_trajectory(args[argname_trg], param.write.free_energy);
  else if (param.write.free_energy)
    io::messages.add("write free energy trajectory but no trg argument",
		     "Out_Configuration",
		     io::message::error);

  if (args.count(argname_bae) > 0)
    block_averaged_energy(args[argname_bae], param.write.block_average);
  else if (param.write.block_average && param.write.energy)
    io::messages.add("write block averaged energy but no bae argument",
		     "Out_Configuration",
		     io::message::error);

  if (param.perturbation.perturbation){
    if (args.count(argname_bag) > 0)
      block_averaged_free_energy(args[argname_bag], 
				 param.write.block_average);
    else if (param.write.block_average && param.write.free_energy)
      io::messages.add("write block averaged free energy "
			"but no bag argument",
		       "Out_Configuration",
		       io::message::error);
  }
  if (args.count(argname_tramd) > 0)
    ramd_trajectory(args[argname_tramd], param.ramd.every);
  else if (param.ramd.fc!=0.0 && param.ramd.every)
    io::messages.add("write RAMD trajectory but no tramd argument",
		     "Out_Configuration",
		     io::message::error);

  if (param.replica.num_T * param.replica.num_l){
    m_replica = true;
  }
  
}

void io::Out_Configuration::write(configuration::Configuration &conf,
				  topology::Topology const &topo,
				  simulation::Simulation const &sim,
				  output_format const form)
{ 
  // standard trajectories
  
  bool constraint_force = sim.param().constraint.solute.algorithm == simulation::constr_shake &&
    sim.param().constraint.solvent.algorithm == simulation::constr_shake;
  
  // check whether a new energy minimum was found
  bool minimum_found = false;
  if (sim.param().write.energy_index > 0) {
    double current_energy = conf.old().energies.get_energy_by_index(sim.param().write.energy_index);
    
    // found a new minimum?
    if (current_energy < minimum_energy) {
      minimum_found = true;
      minimum_energy = current_energy;
    }
  }
  
  if (form == reduced){
    /**
     * set this to true when you print the timestep to the special traj.
     * make sure you don't print it twice. 
     */
    bool special_timestep_printed = false; 

    if(m_every_pos && ((sim.steps() % m_every_pos) == 0 || minimum_found)){
      // don't write starting configuration if analyzing a trajectory
      if (sim.steps() || !sim.param().analyze.analyze){
	_print_timestep(sim, *m_pos_traj);
	
	if (sim.param().write.position_solute_only)
	  _print_positionred(conf, topo,  topo.num_solute_atoms(), *m_pos_traj);
	else
	  _print_positionred(conf, topo,  topo.num_atoms(), *m_pos_traj);
	
	if (conf.boundary_type != math::vacuum)
	  _print_box(conf, *m_pos_traj);
      }
      // a new block begins. let's reset the minimum
      minimum_energy = conf.old().energies.get_energy_by_index(sim.param().write.energy_index);
    }
    
    if (m_every_vel && (sim.steps() % m_every_vel) == 0){
      _print_timestep(sim, *m_vel_traj);
      if (sim.param().write.velocity_solute_only)
        _print_velocityred(conf, topo.num_solute_atoms(), *m_vel_traj);
      else
        _print_velocityred(conf, topo.num_atoms(), *m_vel_traj);
    }
    
    if(m_every_force && ((sim.steps()) % m_every_force) == 0){
      if(sim.steps()){
	_print_old_timestep(sim, *m_force_traj);
        if (sim.param().write.force_solute_only)
	  _print_forcered(conf, topo.num_solute_atoms(), *m_force_traj, constraint_force);
        else
          _print_forcered(conf, topo.num_atoms(), *m_force_traj, constraint_force);
      }
    }
    
    if(m_every_cos_pos && (sim.steps() % m_every_cos_pos) == 0){
      if (!special_timestep_printed) {
	_print_timestep(sim, *m_special_traj);
        special_timestep_printed = true;
      }
       _print_cos_position(conf, topo, *m_special_traj);
    }
    
    if(m_every_energy && ((sim.steps() % m_every_energy) == 0 || minimum_found)){
      if(sim.steps()){
	_print_old_timestep(sim, *m_energy_traj);
	_print_energyred(conf, *m_energy_traj);
	_print_volumepressurered(topo, conf, sim, *m_energy_traj);
      }
    }
    
    if(m_every_free_energy && (sim.steps() % m_every_free_energy) == 0){
      if(sim.steps()){
	_print_old_timestep(sim, *m_free_energy_traj);
	_print_free_energyred(conf, topo, *m_free_energy_traj);
      }
    }

    if (m_every_blockaverage && (sim.steps() % m_every_blockaverage) == 0){
      
      if(m_write_blockaverage_energy){
	if(sim.steps()){
	  _print_old_timestep(sim, *m_blockaveraged_energy);
	  _print_blockaveraged_energyred(conf, *m_blockaveraged_energy);
	  _print_blockaveraged_volumepressurered(conf, sim, *m_blockaveraged_energy);
	}
      }

      if(m_write_blockaverage_free_energy){
	if(sim.steps()){
	  _print_old_timestep(sim, *m_blockaveraged_free_energy);
	  _print_blockaveraged_free_energyred(conf, sim.param().perturbation.dlamt,
					      *m_blockaveraged_free_energy);
	}
      }
      conf.current().averages.block().zero();
    }
    
    if (m_every_ramd && (sim.steps() % m_every_ramd) == 0){
      _print_timestep(sim, *m_ramd_traj);
      _print_ramd(topo, conf, sim, *m_ramd_traj);
    }
    
  }
  else if(form == final && m_final){
    _print_timestep(sim, *m_final_conf);
    _print_position(conf, topo, *m_final_conf);
    _print_lattice_shifts(conf, topo, *m_final_conf);
    
    if (sim.param().polarize.cos)
      _print_cos_position(conf, topo, *m_final_conf);

    if(sim.param().minimise.ntem == 0)
      _print_velocity(conf, topo, *m_final_conf);

    _print_box(conf, *m_final_conf);

    if(sim.param().constraint.solute.algorithm
       == simulation::constr_flexshake){
      _print_flexv(conf, topo, *m_final_conf);
    }

    if(sim.param().stochastic.sd){
      _print_stochastic_integral(conf, topo, *m_final_conf);
    }
    
    if(sim.param().perturbation.perturbation){
      _print_pertdata(topo, *m_final_conf);
    }
    
    if(sim.param().distanceres.distanceres < 0){
      _print_distance_restraint_averages(conf, topo, *m_final_conf);
    }
    
    if (sim.param().posrest.posrest != simulation::posrest_off) {
      _print_position_restraints(sim, topo, *m_final_conf);
    }
    
    if(sim.param().jvalue.mode != simulation::restr_off){
      _print_jvalue(sim.param(), conf, topo, *m_final_conf);
    }
    
    if (sim.param().rottrans.rottrans) {
      _print_rottrans(conf, sim, *m_final_conf);
    }

    if(sim.param().pscale.jrest){
      _print_pscale_jrest(conf, topo, *m_final_conf);
    }
    
    if (sim.param().multibath.nosehoover > 1) {
      _print_nose_hoover_chain_variables(sim.multibath(), *m_final_conf);
    }
    
    // forces and energies still go to their trajectories
    if (m_every_force && ((sim.steps()) % m_every_force) == 0){
      _print_old_timestep(sim, *m_force_traj);
      if (sim.param().write.force_solute_only)
	_print_forcered(conf, topo.num_solute_atoms(), *m_force_traj, constraint_force);
      else
        _print_forcered(conf, topo.num_atoms(), *m_force_traj, constraint_force);
    }

    if(m_every_energy && (sim.steps() % m_every_energy) == 0){
      _print_old_timestep(sim, *m_energy_traj);
      _print_energyred(conf, *m_energy_traj);
      _print_volumepressurered(topo, conf, sim, *m_energy_traj);
    }

    if(m_every_free_energy && (sim.steps() % m_every_free_energy) == 0){
      _print_old_timestep(sim, *m_free_energy_traj);
      _print_free_energyred(conf, topo, *m_free_energy_traj);
    }

  }
  else{

    // not reduced or final (so: decorated)

    if(m_every_pos && (sim.steps() % m_every_pos) == 0){
      _print_timestep(sim, *m_pos_traj);
      _print_position(conf, topo, *m_pos_traj);
      if (conf.boundary_type != math::vacuum)
	_print_box(conf, *m_pos_traj);
    }
    
    if (m_every_vel && (sim.steps() % m_every_vel) == 0){
      _print_timestep(sim, *m_vel_traj);
      _print_velocity(conf, topo, *m_vel_traj);
    }
    
    if(m_every_force && (sim.steps() % m_every_force) == 0){
      if (sim.steps()){
	_print_timestep(sim, *m_force_traj);	
	_print_force(conf, topo, *m_force_traj, constraint_force);
      }
    }
  }

  // done writing!

}

/**
 * write out replicas
 */
void io::Out_Configuration::write_replica
(
 std::vector<util::Replica_Data> & replica_data,
 std::vector<configuration::Configuration> & conf,
 topology::Topology const & topo,
 simulation::Simulation const &sim,
 output_format const form)
{
  // standard trajectories
  if (form == reduced){

    if(m_every_pos && (sim.steps() % m_every_pos) == 0){
      _print_timestep(sim, *m_pos_traj);
      _print_replica_information(replica_data, *m_pos_traj);
      
      for(unsigned int i=0; i<conf.size(); ++i){
	_print_positionred(conf[i], topo,  topo.num_atoms(), *m_pos_traj);
	_print_velocityred(conf[0], topo.num_atoms(), *m_vel_traj);
	
	if (conf[i].boundary_type != math::vacuum)
	  _print_box(conf[i], *m_pos_traj);
      }
    }
  }
  else if(form == final && m_final){
    for(unsigned int i=0; i<conf.size(); ++i){
      
      *m_final_conf << "REPLICAFRAME\n"
		   << std::setw(12) << i + 1
		   << "\nEND\n";

      if (i==0){
	_print_timestep(sim, *m_final_conf);
	_print_replica_information(replica_data, *m_final_conf);
      }
      
      _print_position(conf[i], topo, *m_final_conf);
      if (sim.param().polarize.cos)
        _print_cos_position(conf[i], topo, *m_final_conf);
      _print_velocity(conf[i], topo, *m_final_conf);
      _print_lattice_shifts(conf[i], topo, *m_final_conf);
      _print_box(conf[i], *m_final_conf);
    }
    
    /*
    if(sim.param().jvalue.mode != simulation::restr_off){
      _print_jvalue(sim.param(), conf[0], topo, *m_final_conf);
    }
    if(sim.param().pscale.jrest){
      _print_pscale_jrest(conf[0], topo, *m_final_conf);
    }
    */

    if (sim.param().multibath.nosehoover > 1) {
      _print_nose_hoover_chain_variables(sim.multibath(), *m_final_conf);
    }
  }
  // done writing replicas!
}

void io::Out_Configuration
::final_configuration(std::string name)
{
  if (m_compressed)
    m_final_conf = new io::ogzstream(name.c_str());
  else
    m_final_conf = new std::ofstream(name.c_str());
  
  _print_title(m_title, "final configuration", *m_final_conf);
  m_final = true;
}

void io::Out_Configuration
::trajectory(std::string name, int every)
{
  if (m_compressed) 
    m_pos_traj = new io::ogzstream(name.c_str());
  else
    m_pos_traj = new std::ofstream(name.c_str());
  
  m_every_pos = every;
  _print_title(m_title, "position trajectory", *m_pos_traj);
}

void io::Out_Configuration
::velocity_trajectory(std::string name, int every)
{
  if (m_compressed)
    m_vel_traj = new io::ogzstream(name.c_str());
  else
    m_vel_traj = new std::ofstream(name.c_str());
  
  m_every_vel = every;
  _print_title(m_title, "velocity trajectory", *m_vel_traj);
}

void io::Out_Configuration
::force_trajectory(std::string name, int every)
{
  if (m_compressed)
    m_force_traj = new io::ogzstream(name.c_str());
  else
    m_force_traj = new std::ofstream(name.c_str());
  
  m_every_force = every;
  _print_title(m_title, "force trajectory", *m_force_traj);
}

void io::Out_Configuration
::special_trajectory(std::string name, int every_cos)
{
  if (m_compressed)
    m_special_traj = new io::ogzstream(name.c_str());
  else
    m_special_traj = new std::ofstream(name.c_str());
  
  m_every_cos_pos = every_cos;
  _print_title(m_title, "special trajectory", *m_special_traj);
}

void io::Out_Configuration
::energy_trajectory(std::string name, int every)
{
  if (m_compressed)
    m_energy_traj = new io::ogzstream(name.c_str());
  else
    m_energy_traj = new std::ofstream(name.c_str());
  
  m_every_energy = every;
  _print_title(m_title, "energy trajectory", *m_energy_traj);
}

void io::Out_Configuration
::replica_trajectory(std::string name)
{
  if (m_compressed)
    m_replica_traj = new io::ogzstream(name.c_str());
  else
    m_replica_traj = new std::ofstream(name.c_str());
  
  m_has_replica_traj = true;
  _print_title(m_title, "replica trajectory", *m_replica_traj);
}

void io::Out_Configuration
::free_energy_trajectory(std::string name, int every)
{
  if (m_compressed)
    m_free_energy_traj = new io::ogzstream(name.c_str());
  else
    m_free_energy_traj = new std::ofstream(name.c_str());
  
  m_every_free_energy = every;
  _print_title(m_title, "free energy trajectory", *m_free_energy_traj);
}

void io::Out_Configuration
::block_averaged_energy(std::string name, int every)
{
  if (m_compressed)
    m_blockaveraged_energy = new io::ogzstream(name.c_str());
  else
    m_blockaveraged_energy = new std::ofstream(name.c_str());
  
  if (m_every_blockaverage && m_every_blockaverage != every){
    io::messages.add("overwriting how often block averages are written out illegal",
		     "Out_Configuration",
		     io::message::error);
  }
  m_every_blockaverage = every;
  m_write_blockaverage_energy = true;
  _print_title(m_title, "block averaged energies", *m_blockaveraged_energy);
}

void io::Out_Configuration
::block_averaged_free_energy(std::string name, int every)
{
  if (m_compressed)
    m_blockaveraged_free_energy = new io::ogzstream(name.c_str());
  else
    m_blockaveraged_free_energy = new std::ofstream(name.c_str());
  
  if (m_every_blockaverage && m_every_blockaverage != every){
    io::messages.add("overwriting how often block averages are written out illegal",
		     "Out_Configuration",
		     io::message::error);
  }
  m_every_blockaverage = every;
  m_write_blockaverage_free_energy = true;
  _print_title(m_title, "block averaged free energies", *m_blockaveraged_free_energy);
}

void io::Out_Configuration
::ramd_trajectory(std::string name, int every)
{
  if (m_compressed)
    m_ramd_traj = new io::ogzstream(name.c_str());
  else
    m_ramd_traj = new std::ofstream(name.c_str());
  
  m_every_ramd = every;
  _print_title(m_title, "RAMD trajectory", *m_ramd_traj);
}


void io::Out_Configuration
::_print_timestep(simulation::Simulation const &sim, 
		  std::ostream &os)
{
  os.setf(std::ios::fixed, std::ios::floatfield);
  os.precision(m_precision);
  
  os << "TIMESTEP\n"
     << std::setw(m_width) << sim.steps()
     << " "
     << std::setw(m_width-1) << sim.time()
     << "\nEND\n";
  
}

void io::Out_Configuration
::_print_old_timestep(simulation::Simulation const &sim,
		      std::ostream &os)
{
  os.setf(std::ios::fixed, std::ios::floatfield);
  os.precision(m_precision);
  
  os << "TIMESTEP\n"
     << std::setw(m_width) << sim.steps()-1
     << " "
     << std::setw(m_width-1) << sim.time() - sim.time_step_size()
     << "\nEND\n";
  
}

template<math::boundary_enum b>
void _print_g96_position_bound(configuration::Configuration const &conf,
				      topology::Topology const &topo,
				      std::ostream &os, int width)
{
  math::Periodicity<b> periodicity(conf.current().box);
  math::VArray const &pos = conf.current().pos;
  topology::Solute const &solute = topo.solute();
  std::vector<std::string> const &residue_name = topo.residue_names();

  math::Vec v, v_box, trans, r;
  
  os << "# first 24 chars ignored\n";

  // put chargegroups into the box (on the fly)
  topology::Chargegroup_Iterator cg_it = topo.chargegroup_begin(),
    cg_to = topo.chargegroup_end();

  // solute chargegroups...
  unsigned int i = 0;
  for( ; i < topo.num_solute_chargegroups(); ++cg_it, ++i){
    // gather on first atom...
    v = pos(*cg_it.begin());
    v_box = v;
    periodicity.put_into_positive_box(v_box);
    trans = v_box - v;
    
    // atoms in a chargegroup
    topology::Atom_Iterator at_it = cg_it.begin(),
      at_to = cg_it.end();

    for( ; at_it != at_to; ++at_it){
      r = pos(*at_it) + trans;

      os << std::setw(5)  << solute.atom(*at_it).residue_nr+1 << " "
	 << std::setw(5)  << std::left 
	 << residue_name[solute.atom(*at_it).residue_nr] << " "
	 << std::setw(6)  << std::left << solute.atom(*at_it).name 
	 << std::right
	 << std::setw(6)  << *at_it+1
	 << std::setw(width) << r(0)
	 << std::setw(width) << r(1)
	 << std::setw(width) << r(2)
	 << "\n";
    } 
  }

  // solvent chargegroups
  unsigned int s = 0;
  unsigned int mol = 0;

  for( ; cg_it != cg_to; ++cg_it, ++mol){
    v = pos(**cg_it);
    v_box = v;
    periodicity.put_into_positive_box(v_box);
    trans = v_box - v;
    
    if (mol >= topo.num_solvent_molecules(s)) ++s;
    
    // loop over the atoms
    topology::Atom_Iterator at_it = cg_it.begin(),
      at_to = cg_it.end();
    // one chargegroup per solvent
    for(unsigned int atom=0; at_it != at_to; ++at_it, ++atom){
      r = pos(*at_it) + trans;
	
      os << std::setw(5)  << mol+1
	 << ' ' << std::setw(5)  << std::left
	 << residue_name[topo.solvent(s).atom(atom).residue_nr] << " "
	 << std::setw(6)  << std::left << topo.solvent(s).atom(atom).name
	 << std::right
	 << std::setw(6)  << *at_it + 1
	 << std::setw(width) << r(0)
	 << std::setw(width) << r(1)
	 << std::setw(width) << r(2)
	 << "\n";
    }
  }

}


/**
 * i need a specialized function to put the particles into the box.
 */
template<math::boundary_enum b>
void _print_position_bound(configuration::Configuration const &conf,
				  topology::Topology const &topo,
				  std::ostream &os, int width)
{
  math::Periodicity<b> periodicity(conf.current().box);
  math::VArray const &pos = conf.current().pos;
  topology::Solute const &solute = topo.solute();
  std::vector<std::string> const &residue_name = topo.residue_names();

  math::Vec v;
  
  os << "# first 24 chars ignored\n";
  
  for(int i=0,to = topo.num_solute_atoms(); i<to; ++i){

    v = pos(i);
    periodicity.put_into_box(v);

    os << std::setw(5)  << solute.atom(i).residue_nr+1 << " "
       << std::setw(5)  << std::left 
       << residue_name[solute.atom(i).residue_nr] << " "
       << std::setw(6)  << std::left << solute.atom(i).name << std::right
       << std::setw(6)  << i+1
       << std::setw(width) << v(0)
       << std::setw(width) << v(1)
       << std::setw(width) << v(2)
       << "\n";
  }

  int index = topo.num_solute_atoms();
  int res_nr = 1;

  for(unsigned int s=0; s < topo.num_solvents(); ++s){

    for(unsigned int m=0; m < topo.num_solvent_molecules(s); ++m, ++res_nr){
      
      for(unsigned int a=0; a < topo.solvent(s).num_atoms(); ++a, ++index){
	
	v = pos(index);
	periodicity.put_into_positive_box(v);
	
	os << std::setw(5)  << res_nr
	   << ' ' << std::setw(5)  << std::left
	   << residue_name[topo.solvent(s).atom(a).residue_nr] << " "
	   << std::setw(6)  << std::left 
	   << topo.solvent(s).atom(a).name << std::right
	   << std::setw(6)  << index + 1
	   << std::setw(width) << v(0)
	   << std::setw(width) << v(1)
	   << std::setw(width) << v(2)
	   << "\n";
      }
    }
  }
}


void io::Out_Configuration
::_print_position(configuration::Configuration const &conf,
		  topology::Topology const &topo, 
		  std::ostream &os)
{
  os.setf(std::ios::fixed, std::ios::floatfield);
  os.precision(m_precision);
  
  os << "POSITION\n";

  SPLIT_BOUNDARY(_print_g96_position_bound,
		 conf, topo, os, m_width);
  
  os << "END\n";
  
}

template<math::boundary_enum b>
void _print_g96_positionred_bound(configuration::Configuration const &conf,
				  topology::Topology const &topo,
				  int num,
				  std::ostream &os, int width)
{
  DEBUG(10, "g96 positionred");
  
  math::Periodicity<b> periodicity(conf.current().box);
  math::VArray const &pos = conf.current().pos;

  math::Vec v, v_box, trans, r;
  
  assert(num >= 0);
  
  // put chargegroups into the box (on the fly)
  topology::Chargegroup_Iterator cg_it = topo.chargegroup_begin(),
    cg_to = topo.chargegroup_end();
  DEBUG(10, "cg to : " <<  **cg_to << std::endl);
  
  // solute chargegroups...
  unsigned int i = 0, count = 0;
  for( ; i < topo.num_solute_chargegroups(); ++cg_it, ++i){
    DEBUG(10, "solute cg: " << i);
    // gather on first atom...
    v = pos(*cg_it.begin());
    v_box = v;
    periodicity.put_into_positive_box(v_box);
    trans = v_box - v;
    
    // atoms in a chargegroup
    topology::Atom_Iterator at_it = cg_it.begin(),
      at_to = cg_it.end();

    for( ; at_it != at_to; ++at_it, ++count){

      if (*at_it >= unsigned(num)) return;

      DEBUG(10, "atom: " << count);
      r = pos(*at_it) + trans;

      os << std::setw(width) << r(0)
	 << std::setw(width) << r(1)
	 << std::setw(width) << r(2)
	 << "\n";

      if ((count+1) % 10 == 0) os << '#' << std::setw(10) << count+1 << "\n";
    
    } 
  }

  DEBUG(10, "solvent");
  
  // solvent chargegroups
  unsigned int mol = 0;

  for( ; cg_it != cg_to; ++cg_it, ++mol){
    DEBUG(10, "solvent " << mol);
    
    v = pos(**cg_it);
    v_box = v;
    periodicity.put_into_positive_box(v_box);
    trans = v_box - v;
    
    // loop over the atoms
    topology::Atom_Iterator at_it = cg_it.begin(),
      at_to = cg_it.end();
    // one chargegroup per solvent
    for( ; at_it != at_to; ++at_it, ++count){
      DEBUG(10, "\tatom " << count);
      
      if (*at_it >= unsigned(num)) return;

      r = pos(*at_it) + trans;
	
      os << std::setw(width) << r(0)
	 << std::setw(width) << r(1)
	 << std::setw(width) << r(2)
	 << "\n";

      if ((count+1) % 10 == 0) os << '#' << std::setw(10) << count+1 << "\n";

    }
  }
}


template<math::boundary_enum b>
void
_print_positionred_bound(configuration::Configuration const &conf,
			 int num,
			 std::ostream &os, int width)
{
  math::Periodicity<b> periodicity(conf.current().box);
  
  math::VArray const &pos = conf.current().pos;
  math::Vec v;

  DEBUG(10, "writing POSITIONRED " << pos.size() );
  
  for(int i=0; i<num; ++i){

    v = pos(i);
    periodicity.put_into_box(v);

    os << std::setw(width) << v(0)
       << std::setw(width) << v(1)
       << std::setw(width) << v(2)
       << "\n";

    if((i+1)% 10 == 0) os << '#' << std::setw(10) << i+1 << "\n";
  }
  
}

inline void io::Out_Configuration
::_print_positionred(configuration::Configuration const &conf,
		     topology::Topology const &topo,
		     int num,
		     std::ostream &os)
{
  os.setf(std::ios::fixed, std::ios::floatfield);
  os.precision(m_precision);
  
  os << "POSITIONRED\n";
  DEBUG(7, "configuration boundary type :" << conf.boundary_type);

  SPLIT_BOUNDARY(_print_g96_positionred_bound, conf, topo, num, os, m_width);
  
  os << "END\n";

}

inline void io::Out_Configuration
::_print_cos_position(configuration::Configuration const &conf,
		      topology::Topology const &topo,
		      std::ostream &os)
{
  os.setf(std::ios::fixed, std::ios::floatfield);
  os.precision(m_precision);
  
  os << "COSDISPLACEMENTS\n";
  
  math::VArray const &posV = conf.current().posV;
  
  for(unsigned int i=0; i<posV.size(); ++i){

    os << std::setw(m_width) << posV(i)(0)
       << std::setw(m_width) << posV(i)(1)
       << std::setw(m_width) << posV(i)(2)
       << "\n";

    if((i+1)% 10 == 0) os << '#' << std::setw(10) << i+1 << "\n";
  }  
  
  os << "END\n";
}

void io::Out_Configuration
::_print_velocity(configuration::Configuration const &conf,
		  topology::Topology const &topo, 
		  std::ostream &os)
{
  os.setf(std::ios::fixed, std::ios::floatfield);
  os.precision(m_precision);
  
  os << "VELOCITY\n";
  
  math::VArray const &vel = conf.current().vel;
  topology::Solute const &solute = topo.solute();
  std::vector<std::string> const &residue_name = topo.residue_names();

  os << "# first 24 chars ignored\n";
  
  for(int i=0,to = topo.num_solute_atoms(); i<to; ++i){

    os << std::setw(5)  << solute.atom(i).residue_nr+1 << " "
       << std::setw(5)  << std::left << residue_name[solute.atom(i).residue_nr] << " "
       << std::setw(6)  << std::left << solute.atom(i).name << std::right
       << std::setw(6)  << i+1
       << std::setw(m_width) << vel(i)(0)
       << std::setw(m_width) << vel(i)(1)
       << std::setw(m_width) << vel(i)(2)
       << "\n";
  }
  
  int index = topo.num_solute_atoms();
  int res_num = 1;
  
  for(unsigned int s=0; s < topo.num_solvents(); ++s){

    for(unsigned int m=0; m < topo.num_solvent_molecules(s); ++m, ++res_num){
      
      for(unsigned int a=0; a < topo.solvent(s).num_atoms(); ++a, ++index){
	
	os << std::setw(5)  << res_num << " "
	   << std::setw(5)  << std::left 
	   << residue_name[topo.solvent(s).atom(a).residue_nr] << " "
	   << std::setw(6)  << std::left << topo.solvent(s).atom(a).name << std::right
	   << std::setw(6)  << index + 1
	   << std::setw(m_width) << vel(index)(0)
	   << std::setw(m_width) << vel(index)(1)
	   << std::setw(m_width) << vel(index)(2)
	   << "\n";
      }
    }
  }

  os << "END\n";

}

void io::Out_Configuration
::_print_lattice_shifts(configuration::Configuration const &conf,
		  topology::Topology const &topo, 
		  std::ostream &os)
{  
  os << "LATTICESHIFTS\n";
  math::VArray const &shift = conf.special().lattice_shifts;
  for(int i=0,to = topo.num_atoms(); i < to; ++i){
    os << std::setw(10) << int(rint(shift(i)(0)))
       << std::setw(10) << int(rint(shift(i)(1)))
       << std::setw(10) << int(rint(shift(i)(2)))
       << "\n";
  }
  os << "END\n";
}

void io::Out_Configuration
::_print_velocityred(configuration::Configuration const &conf,
	             int num,
		     std::ostream &os)
{
  os.setf(std::ios::fixed, std::ios::floatfield);
  os.precision(m_precision);
  
  os << "VELOCITYRED\n";
  
  math::VArray const &vel = conf.current().vel;
  
  assert(num <= int(vel.size()));
  
  for(int i=0; i<num; ++i){

    os << std::setw(m_width) << vel(i)(0)
       << std::setw(m_width) << vel(i)(1)
       << std::setw(m_width) << vel(i)(2)
       << "\n";

    if((i+1)% 10 == 0) os << '#' << std::setw(10) << i+1 << "\n";
  }
  
  os << "END\n";
  
}

void io::Out_Configuration
::_print_force(configuration::Configuration const &conf,
	       topology::Topology const &topo,
	       std::ostream &os, bool constraint_force)
{
  os.setf(std::ios::fixed, std::ios::floatfield);
  os.precision(m_force_precision);
  
  os << "FREEFORCE\n";
  
  math::VArray const & force = conf.current().force;
  topology::Solute const &solute = topo.solute();
  std::vector<std::string> const &residue_name = topo.residue_names();

  os << "# first 24 chars ignored\n";
  for(int i=0,to = topo.num_solute_atoms(); i<to; ++i){
    os << std::setw(6)  << solute.atom(i).residue_nr+1
       << std::setw(5)  << residue_name[solute.atom(i).residue_nr]
       << std::setw(6)  << solute.atom(i).name
       << std::setw(8)  << i+1
       << std::setw(m_force_width) << force(i)(0)
       << std::setw(m_force_width) << force(i)(1)
       << std::setw(m_force_width) << force(i)(2)
       << "\n";
  }
  int index = topo.num_solute_atoms();
  
  for(unsigned int s=0; s < topo.num_solvents(); ++s){
    for(unsigned int m=0; m < topo.num_solvent_molecules(s); ++m){
      for(unsigned int a=0; a < topo.solvent(s).num_atoms(); ++a, ++index){
	os << std::setw(6)  << topo.solvent(s).atom(a).residue_nr+1
	   << std::setw(5)  << residue_name[topo.solvent(s).atom(a).residue_nr]
	   << std::setw(6)  << topo.solvent(s).atom(a).name
	   << std::setw(8)  << index + 1
	   << std::setw(m_force_width) << force(index)(0)
	   << std::setw(m_force_width) << force(index)(1)
	   << std::setw(m_force_width) << force(index)(2)
	   << "\n";
      }
    }
  }
  os << "END\n";

  if (constraint_force) {
    os << "CONSFORCE\n";

    math::VArray const & cons_force = conf.current().constraint_force;

    os << "# first 24 chars ignored\n";
    for (int i = 0, to = topo.num_solute_atoms(); i < to; ++i) {
      os << std::setw(6) << solute.atom(i).residue_nr + 1
              << std::setw(5) << residue_name[solute.atom(i).residue_nr]
              << std::setw(6) << solute.atom(i).name
              << std::setw(8) << i + 1
              << std::setw(m_force_width) << cons_force(i)(0)
              << std::setw(m_force_width) << cons_force(i)(1)
              << std::setw(m_force_width) << cons_force(i)(2)
              << "\n";
    }
    index = topo.num_solute_atoms();

    for (unsigned int s = 0; s < topo.num_solvents(); ++s) {
      for (unsigned int m = 0; m < topo.num_solvent_molecules(s); ++m) {
        for (unsigned int a = 0; a < topo.solvent(s).num_atoms(); ++a, ++index) {
          os << std::setw(6) << topo.solvent(s).atom(a).residue_nr + 1
                  << std::setw(5) << residue_name[topo.solvent(s).atom(a).residue_nr]
                  << std::setw(6) << topo.solvent(s).atom(a).name
                  << std::setw(8) << index + 1
                  << std::setw(m_force_width) << cons_force(index)(0)
                  << std::setw(m_force_width) << cons_force(index)(1)
                  << std::setw(m_force_width) << cons_force(index)(2)
                  << "\n";
        }
      }
    }
    os << "END\n";
  }
}

void io::Out_Configuration
::_print_forcered(configuration::Configuration const &conf,
                  int num,
		  std::ostream &os,
                  bool constraint_force)
{
  os.setf(std::ios::fixed, std::ios::floatfield);
  os.precision(m_force_precision);
  
  os << "FREEFORCERED\n";
  
  const math::VArray &force = conf.old().force;
  assert(num <= int(force.size()));
  
  for(int i=0; i<num; ++i){
    os << std::setw(m_force_width) << force(i)(0)
       << std::setw(m_force_width) << force(i)(1)
       << std::setw(m_force_width) << force(i)(2)
       << "\n";

    if((i+1)% 10 == 0) os << '#' << std::setw(10) << i+1 << "\n";
  }
  
  os << "END\n";

  if (constraint_force) {
    os << "CONSFORCERED\n";

    const math::VArray & cons_force = conf.old().constraint_force;
    assert(num <= int(cons_force.size()));

    for (int i = 0; i < num; ++i) {
      os << std::setw(m_force_width) << cons_force(i)(0)
              << std::setw(m_force_width) << cons_force(i)(1)
              << std::setw(m_force_width) << cons_force(i)(2)
              << "\n";

      if ((i + 1) % 10 == 0) os << '#' << std::setw(10) << i + 1 << "\n";
    }

    os << "END\n";
  }
}

void io::Out_Configuration
::_print_energyred(configuration::Configuration const &conf,
		   std::ostream &os)
{
  os.setf(std::ios::scientific, std::ios::floatfield);
  os.precision(m_precision);
  
  os << "ENERGY03\n";
  _print_energyred_helper(os, conf.old().energies);
  os << "END\n";
  
}

void io::Out_Configuration
::write_replica_energy(util::Replica_Data const & replica_data,
		       simulation::Simulation const & sim,
		       configuration::Energy const & energy,
		       int reeval)
{
  std::ostream &replica_traj = *m_replica_traj;
  replica_traj.setf(std::ios::scientific, std::ios::floatfield);
  replica_traj.precision(m_precision);
  
  print_REMD(replica_traj, replica_data, sim.param(), reeval);
  _print_timestep(sim, replica_traj);

  replica_traj << "ENERGY03\n";
  _print_energyred_helper(replica_traj, energy);
  replica_traj << "END\n";
  
}

void io::Out_Configuration
::_print_free_energyred(configuration::Configuration const &conf,
			topology::Topology const &topo,
			std::ostream &os)
{
  // assert(m_free_energy_traj.is_open());
  
  os.setf(std::ios::scientific, std::ios::floatfield);
  os.precision(m_precision);

  os << "FREEENERDERIVS03\n"
     << "# lambda\n"
     << std::setw(18) << topo.old_lambda() << "\n";

  _print_energyred_helper(os, conf.old().perturbed_energy_derivatives);

  os << "END\n";

}

void io::Out_Configuration
::_print_volumepressurered(topology::Topology const & topo,
			   configuration::Configuration const &conf,
			   simulation::Simulation const &sim,
			   std::ostream &os)
{
  std::vector<double> const s;

  os.setf(std::ios::scientific, std::ios::floatfield);
  os.precision(m_precision);

  os << "VOLUMEPRESSURE03\n";

  _print_volumepressurered_helper(os,
				  math::sum(topo.mass()),
				  sim.multibath(),
				  s,
				  conf.old().energies,
				  conf.old().box,
				  conf.boundary_type,
				  conf.old().pressure_tensor,
				  conf.old().virial_tensor,
				  conf.old().kinetic_energy_tensor);
  
  os << "END\n";
  
}

template<math::boundary_enum b>
void _print_ramd_bound(topology::Topology const & topo,
		       configuration::Configuration const &conf,
		       simulation::Simulation const &sim,
		       std::ostream &os, int width)
{
  
  math::Periodicity<b> periodicity(conf.current().box);
  math::Vec com(0.0,0.0,0.0);
  math::Vec r;
  math::Vec f = conf.special().ramd.force_direction * sim.param().ramd.fc;
  std::set<unsigned int>::const_iterator
    it = sim.param().ramd.atom.begin(),
    i0 = sim.param().ramd.atom.begin(),
    to = sim.param().ramd.atom.end();

  for(; it!=to; ++it){
    periodicity.nearest_image(conf.current().pos(*it), 
			      conf.current().pos(*i0), r);
    com += topo.mass()(*it) * r;
  }
  
  com /= conf.special().ramd.total_mass;
  com += conf.current().pos(*i0);
  
  os << "# force\n";
  os << std::setw(width) << f(0)
     << std::setw(width) << f(1)
     << std::setw(width) << f(2)
     << "\n";
  os << "# com RAMD atoms\n";
  os << std::setw(width) << com(0)
     << std::setw(width) << com(1)
     << std::setw(width) << com(2)
     << "\n";
  
}

void io::Out_Configuration
::_print_ramd(topology::Topology const &topo,
	      configuration::Configuration const &conf,
	      simulation::Simulation const &sim,
	      std::ostream &os)
{
  os.setf(std::ios::fixed, std::ios::floatfield);
  os.precision(m_precision);
  os << "RAMD\n";
  
  SPLIT_BOUNDARY(_print_ramd_bound,
		 topo, conf, sim, os, m_width);
  os << "END\n";
  
}

  
void io::Out_Configuration
::_print_box(configuration::Configuration const &conf,
	     std::ostream &os)
{
  os.setf(std::ios::fixed, std::ios::floatfield);
  os.precision(m_precision);
  //change to GENBOX
  os << "GENBOX\n";
  
  math::Box const &box = conf.current().box;
  
  os << std::setw(5) << conf.boundary_type << "\n";
  
  long double a, b, c, alpha, beta, gamma, phi, theta, psi;
  
  a=math::abs(box(0));
  b=math::abs(box(1));
  c=math::abs(box(2));

  os << std::setw(m_width) << a
     << std::setw(m_width) << b
     << std::setw(m_width) << c
     << "\n";
  
  alpha = acos(math::costest(dot(box(1),box(2))/(abs(box(1))*abs(box(2))))); 
  beta  = acos(math::costest(dot(box(0),box(2))/(abs(box(0))*abs(box(2)))));
  gamma = acos(math::costest(dot(box(0),box(1))/(abs(box(0))*abs(box(1)))));
  
  os << std::setw(m_width) << alpha*180/math::Pi
     << std::setw(m_width) << beta*180/math::Pi
     << std::setw(m_width) << gamma*180/math::Pi
     << "\n";
  
  long double cosdelta=(cosl(alpha)-cosl(beta)*cosl(gamma))/(sinl(beta)*sinl(gamma)); 
  
  long double sindelta=sqrtl(1-cosdelta*cosdelta); 
  long double cotdelta=cosdelta/sindelta;
  
  long double cotgamma, cotbeta;
  //check for cot(gamma)
  if(gamma==90){
        cotgamma=0;
  }
  else{
        cotgamma=1/(tanl(gamma));
  }
  //check for cot(beta)
  if(beta==90){
        cotbeta=0;
  }
  else{
        cotbeta=1/(tanl(beta));
  }

  math::Vecl BSx(1/a,0.0,0.0);
  math::Vecl BSy(-cotdelta/a,
          1/(b*sinl(gamma)),0.0);
  math::Vecl BSz((cotdelta*cotgamma-cotbeta/sindelta)/a, 
          -cotdelta/(b*sinl(gamma)),
          1/(c*sinl(beta)*sindelta));
  
  
  math::Matrixl BSmat(BSx,BSy,BSz);
  math::Matrix boxmat(box(0),box(1),box(2));

  math::Matrixl Rmat=product(BSmat,boxmat);
 /* 
  os <<"boxmat\n" ;
  for (int i=0;i<3; i++){
  os << std::setw(m_width) << box(0)(i)
     << std::setw(m_width) << box(1)(i)
     << std::setw(m_width) << box(2)(i)
    << "\n";  
  }

    os <<"Bsmat\n" ;
  for (int i=0;i<3; i++){
  os << std::setw(m_width) << BSmat(0,i)
     << std::setw(m_width) << BSmat(1,i)
     << std::setw(m_width) << BSmat(2,i)
     << "\n";  
  }
  
  os <<"Rmat\n" ;
  for (int i=0;i<3; i++){
    os << std::setw(m_width) << Rmat(0,i)
     << std::setw(m_width) << Rmat(1,i)
     << std::setw(m_width) << Rmat(2,i)
     << "\n";  
  }
 */   
  long double R11R21 = sqrtl(Rmat(0,0)*Rmat(0,0)+Rmat(0,1)*Rmat(0,1));
  if(R11R21==0.0)
  {
      theta = -math::sign(Rmat(0,2))*M_PI/2;
      psi   = 0.0;
      phi   =-math::sign(Rmat(1,0))*acosl(math::costest(Rmat(1,1)));
  }
  else
  {
      theta = -math::sign(Rmat(0,2))*acosl(math::costest(R11R21));
      long double costheta=cosl(theta);
      psi   = math::sign(Rmat(1,2)/costheta)*acosl(math::costest(Rmat(2,2)/costheta));
      phi   = math::sign(Rmat(0,1)/costheta)*acosl(math::costest(Rmat(0,0)/costheta));

  }
  
 
  
  if(theta==-0) theta =0;
  if(psi==-0)   psi=0;
  if(phi==-0)   phi=0;
          
          
  os << std::setw(m_width) << phi*180/math::Pi
     << std::setw(m_width) << theta*180/math::Pi
     << std::setw(m_width) << psi*180/math::Pi
     << "\n";
  
  double origin=0.0;
  
  os << std::setw(m_width) << origin
     << std::setw(m_width) << origin
     << std::setw(m_width) << origin
     << "\n";

  
  os << "END\n";
  
}

void io::Out_Configuration
::precision(int prec, int add)
{
  m_precision = prec;
  m_width = prec + add;
}

void io::Out_Configuration
::force_precision(int prec, int add)
{
  m_force_precision = prec;
  m_force_width = prec + add;
}

int io::Out_Configuration
::precision()
{
  return m_precision;
}

int io::Out_Configuration
::force_precision()
{
  return m_force_precision;
}


void io::Out_Configuration
::print(topology::Topology const & topo,
	configuration::Configuration & conf,
	simulation::Simulation const & sim)
{
  if (sim.param().print.stepblock && (sim.steps() % sim.param().print.stepblock) == 0){

    m_output << "\n---------------------------------------------------"
	     << "-----------------------------\n";
    
    _print_timestep(sim, m_output);
    
    print_ENERGY(m_output, conf.old().energies, topo.energy_groups());
    
    if (sim.param().perturbation.perturbation){
      
      m_output << "lambda: " << topo.old_lambda() << "\n";
      
      print_ENERGY(m_output, conf.old().perturbed_energy_derivatives, 
		   topo.energy_groups(), "dE/dLAMBDA", "dE_");
    }
    
    print_MULTIBATH(m_output, sim.multibath(), conf.old().energies);

    // flexible shake kinetic energy
    if (sim.param().constraint.solute.algorithm == simulation::constr_flexshake){
      m_output << "FLEXSHAKE\n";
      m_output << "\tflex_ekin";
      for(unsigned int i=0; i < conf.special().flexible_constraint.flexible_ekin.size(); ++i)
	m_output << std::setw(12) << std::setprecision(4) << std::scientific
		 << conf.special().flexible_constraint.flexible_ekin[i];
      m_output << "\nEND\n";
    }
    
    if (sim.param().pcouple.calculate)
      print_PRESSURE(m_output, conf);
    
    m_output.flush();

  }
  if (sim.param().ramd.fc!=0.0 && 
      sim.param().ramd.every &&
      (sim.steps() % sim.param().ramd.every) == 0){
    print_RAMD(m_output, conf, topo.old_lambda());
  }
}

void io::Out_Configuration
::print_final(topology::Topology const & topo,
	      configuration::Configuration & conf,
	      simulation::Simulation const & sim)
{
  m_output << "\n============================================================\n";
  m_output << "FINAL DATA\n";
  m_output << "============================================================\n\n\n";
  
  m_output << "\tsimulation time  :" << std::setw(10) << sim.time() << "\n"
	   << "\tsimulation steps :" << std::setw(10) << sim.steps() << "\n\n";

  configuration::Energy e, ef;
  math::Matrix p, pf, v, vf, et, etf;
  
  if (sim.param().minimise.ntem){
    print_ENERGY(m_output, conf.current().energies, topo.energy_groups(), "MINIMIZED ENERGY",
		 "<EMIN>_");
  }

  // new averages
  conf.current().averages.simulation().energy_average(e, ef);
  conf.current().averages.simulation().pressure_average(p, pf, v, vf, et, etf);

  print_ENERGY(m_output, e, topo.energy_groups(), "ENERGY AVERAGES", "<E>_");
  m_output << "\n";
  print_ENERGY(m_output, ef, topo.energy_groups(), "ENERGY FLUCTUATIONS", "<<E>>_");
  m_output << "\n";

  m_output << "\n";
  print_MULTIBATH(m_output, sim.multibath(), e, "TEMPERATURE AVERAGES");
  m_output << "\n";
  print_MULTIBATH(m_output, sim.multibath(), ef, "TEMPERATURE FLUCTUATIONS");

  m_output << "\n\n";
  if (sim.param().pcouple.calculate){
    print_MATRIX(m_output, p, "PRESSURE AVERAGE");
    m_output << "\n";
    print_MATRIX(m_output, pf, "PRESSURE FLUCTUATION");
  }
  m_output << "\n\n";    

  if (sim.param().perturbation.perturbation){

    double lambda, lambda_fluct;
    conf.current().averages.simulation().
      energy_derivative_average(e, ef, lambda, lambda_fluct, sim.param().perturbation.dlamt);

    if (sim.param().perturbation.dlamt){
	
      print_ENERGY(m_output, e, topo.energy_groups(), "CUMULATIVE DG", "DG_");
      m_output << "\n";

      // what's that anyway...
      //print_ENERGY(m_output, ef, topo.energy_groups(), "DG FLUCTUATIONS", "<<DG>>_");
    }

    else{

      std::ostringstream ss, pre;
      ss << "dE/dLAMBDA ";
      pre << "dE/dl";
      
      print_ENERGY(m_output, e, topo.energy_groups(),
		   ss.str() + "AVERAGES", "<" + pre.str() + ">_");
      
      m_output << "\n";
      
      print_ENERGY(m_output, ef, topo.energy_groups(),
		   ss.str() + "FLUCTUATIONS", "<<" + pre.str() + ">>_");
      
      m_output << "\n";
    }
  }
  if(sim.param().ramd.fc!=0.0 && sim.param().ramd.every)
    print_RAMD(m_output, conf, topo.old_lambda());

}

void io::Out_Configuration
::_print_blockaveraged_energyred(configuration::Configuration const &conf,
				 std::ostream &os)
{
  os.setf(std::ios::scientific, std::ios::floatfield);
  os.precision(m_precision);
  
  configuration::Energy e, ef;
  // energies are in old(), but averages in current()!
  conf.current().averages.block().energy_average(e, ef);

  os << "BAENERGY03\n";
  _print_energyred_helper(os, e);
  os << "END\n";

  os << "BAEFLUCT03\n";
  _print_energyred_helper(os, ef);
  os << "END\n";
  
}

void io::Out_Configuration
::_print_blockaveraged_volumepressurered(configuration::Configuration const & conf,
					 simulation::Simulation const & sim,
					 std::ostream &os)
{
  double mass, massf, vol, volf;
  std::vector<double> s, sf;
  configuration::Energy e, ef;
  math::Box b, bf;
  math::Matrix p, pf, v, vf, k, kf;
  
  // needed again. do in 1 function (energyred && volumepressurered)
  // energies in old(), but averages in current()!
  conf.current().averages.block().energy_average(e, ef);
  conf.current().averages.block().pressure_average(p, pf, v, vf, k, kf);
  conf.current().averages.block().mbs_average(mass, massf, vol, volf, b, bf, s, sf);

  os.setf(std::ios::scientific, std::ios::floatfield);
  os.precision(m_precision);

  os << "BAVOLUMEPRESSURE03\n";

  _print_volumepressurered_helper(os,
				  mass,
				  sim.multibath(),
				  s,
				  e,
				  b,
				  conf.boundary_type,
				  p,
				  v,
				  k);
  
  os << "END\n";

  os << "BAVPFLUCT03\n";

  _print_volumepressurered_helper(os,
				  massf,
				  sim.multibath(),
				  sf,
				  ef,
				  bf,
				  conf.boundary_type,
				  pf,
				  vf,
				  kf);
  
  os << "END\n";

}

void io::Out_Configuration
::_print_blockaveraged_free_energyred(configuration::Configuration const &conf,
				      double dlamt,
				      std::ostream &os)
{
  os.setf(std::ios::scientific, std::ios::floatfield);
  os.precision(m_precision);

  configuration::Energy e, ef;
  double lambda, lambda_fluct;
  
  // energies in old(), but averages in current()!
  conf.current().averages.block().energy_derivative_average(e, ef, lambda, lambda_fluct, dlamt);

  os << "BAFREEENERDERIVS03\n"
     << "# lambda\n"
     << std::setw(18) << lambda << "\n";

  _print_energyred_helper(os, e);

  os << "END\n";

  os << "BAFREEEFLUCTS03\n"
     << "# lambda\n"
     << std::setw(18) << lambda_fluct << "\n";

  _print_energyred_helper(os, ef);

  os << "END\n";

}

void io::Out_Configuration::_print_flexv(configuration::Configuration const &conf,
					 topology::Topology const &topo,
					 std::ostream &os)
{
  DEBUG(10, "FLEXV");
  
  unsigned int k = 0;
  
  std::vector<topology::two_body_term_struct>::const_iterator
    constr_it = topo.solute().distance_constraints().begin(),
    constr_to = topo.solute().distance_constraints().end();
  
  os << "FLEXV\n";
  os << "#\tflexible constraints (" 
     << topo.solute().distance_constraints().size()
     << ")\n";

  for( ; constr_it != constr_to; ++constr_it, ++k){
    
    assert(conf.special().flexible_constraint.flex_len.size() > k);
    assert(conf.special().flexible_constraint.flexible_vel.size() > k);
    
    os << std::setw(15) << constr_it->i+1
       << std::setw(10) << constr_it->j+1
       << std::setw(20) << conf.special().flexible_constraint.flex_len[k]
       << std::setw(20) << conf.special().flexible_constraint.flexible_vel[k]
       << "\n";
  }

  std::vector<topology::perturbed_two_body_term_struct>::const_iterator
    pconstr_it = topo.perturbed_solute().distance_constraints().begin(),
    pconstr_to = topo.perturbed_solute().distance_constraints().end();

  os << "#\tperturbed flexible constraints (" 
     << topo.perturbed_solute().distance_constraints().size()
     << " of "
     << conf.special().flexible_constraint.flex_len.size()
     << ")\n";

  for( ; pconstr_it != pconstr_to; ++pconstr_it, ++k){

    assert(conf.special().flexible_constraint.flex_len.size() > k);
    assert(conf.special().flexible_constraint.flexible_vel.size() > k);
    
    os << std::setw(15) << pconstr_it->i+1
       << std::setw(10) << pconstr_it->j+1
       << std::setw(20) << conf.special().flexible_constraint.flex_len[k]
       << std::setw(20) << conf.special().flexible_constraint.flexible_vel[k]
       << "\n";
  }
  
  os << "END\n";

}

void io::Out_Configuration::_print_stochastic_integral(configuration::Configuration const &conf,
						      topology::Topology const &topo,
						      std::ostream &os)
{
  topology::Solute const &solute = topo.solute();
  std::vector<std::string> const &residue_name = topo.residue_names();

  os.setf(std::ios::fixed, std::ios::floatfield);
  os.precision(m_precision - 2); // because it's scientific

  os << "STOCHINT\n";
  os << "# first 24 chars ignored\n";
  
  for(int i=0,to = topo.num_solute_atoms(); i<to; ++i){
    os.setf(std::ios::fixed, std::ios::floatfield);
    os << std::setw(5)  << solute.atom(i).residue_nr+1 << " "
       << std::setw(5)  << std::left 
       << residue_name[solute.atom(i).residue_nr] << " "
       << std::setw(6)  << std::left << solute.atom(i).name << std::right
       << std::setw(6)  << i+1;
    os.setf(std::ios::scientific, std::ios::floatfield);
    os << std::setw(m_width) << conf.current().stochastic_integral(i)(0)
       << std::setw(m_width) << conf.current().stochastic_integral(i)(1)
       << std::setw(m_width) << conf.current().stochastic_integral(i)(2)
       << "\n";
  }

  int index = topo.num_solute_atoms();
  int res_nr = 1;

  for(unsigned int s=0; s < topo.num_solvents(); ++s){

    for(unsigned int m=0; m < topo.num_solvent_molecules(s); ++m, ++res_nr){
      
      for(unsigned int a=0; a < topo.solvent(s).num_atoms(); ++a, ++index){
        os.setf(std::ios::fixed, std::ios::floatfield);
	os << std::setw(5)  << res_nr
	   << ' ' << std::setw(5)  << std::left
	   << residue_name[topo.solvent(s).atom(a).residue_nr] << " "
	   << std::setw(6)  << std::left 
	   << topo.solvent(s).atom(a).name << std::right
	   << std::setw(6)  << index + 1;
        os.setf(std::ios::scientific, std::ios::floatfield);
	os << std::setw(m_width) << conf.current().stochastic_integral(index)(0)
	   << std::setw(m_width) << conf.current().stochastic_integral(index)(1)
	   << std::setw(m_width) << conf.current().stochastic_integral(index)(2)
	   << "\n";
      }
    }
  }
  os.setf(std::ios::fixed, std::ios::floatfield);
  os << "# seed\n" << std::setw(10) << std::right 
     << conf.current().stochastic_seed << "\n";
  os << "END\n";
}

void io::Out_Configuration::_print_pertdata(topology::Topology const &topo,
					    std::ostream &os)
{
  os.setf(std::ios::fixed, std::ios::floatfield);
  os.precision(7);
  os << "PERTDATA" << std::endl
     << std::setw(15) << topo.lambda() << std::endl
     << "END" << std::endl;
}

void io::Out_Configuration::write_replica_step
(
 simulation::Simulation const & sim,
 util::Replica_Data const & replica_data,
 output_format const form
 )
{
  DEBUG(10, "REPLICA");
  
  // print to all trajectories
  if (form == reduced){

    if(m_every_pos && (sim.steps() % m_every_pos) == 0)
      print_REMD(*m_pos_traj, replica_data, sim.param());
    
    if (m_every_vel && (sim.steps() % m_every_vel) == 0)
      print_REMD(*m_vel_traj, replica_data, sim.param());
    
    if(m_every_force && ((sim.steps()) % m_every_force) == 0){
      if(sim.steps())
	print_REMD(*m_force_traj, replica_data, sim.param());
    }
    
    if(m_every_energy && (sim.steps() % m_every_energy) == 0){
      if(sim.steps())
	print_REMD(*m_energy_traj, replica_data, sim.param());
    }
    
    if(m_every_free_energy && (sim.steps() % m_every_free_energy) == 0){
      if(sim.steps())
	print_REMD(*m_free_energy_traj, replica_data, sim.param());
    }

    if (m_every_blockaverage && (sim.steps() % m_every_blockaverage) == 0){
      if(m_write_blockaverage_energy){
	if(sim.steps())
	  print_REMD(*m_blockaveraged_energy, replica_data, sim.param());
      }

      if(m_write_blockaverage_free_energy){
	if(sim.steps())
	  print_REMD(*m_blockaveraged_free_energy, replica_data, sim.param());
      }
    }
  } // reduced

  else if(form == final && m_final){
    print_REMD(*m_final_conf, replica_data, sim.param());
    
    // forces and energies still go to their trajectories
    if (m_every_force && ((sim.steps()) % m_every_force) == 0)
      if (sim.steps())
        print_REMD(*m_force_traj, replica_data, sim.param());

    if(m_every_energy && (sim.steps() % m_every_energy) == 0)
      print_REMD(*m_energy_traj, replica_data, sim.param());

    if(m_every_free_energy && (sim.steps() % m_every_free_energy) == 0)
      print_REMD(*m_free_energy_traj, replica_data, sim.param());
  } // final

  else{

    // not reduced or final (so: decorated)

    if(m_every_pos && (sim.steps() % m_every_pos) == 0)
      print_REMD(*m_pos_traj, replica_data, sim.param());
    
    if (m_every_vel && (sim.steps() % m_every_vel) == 0)
      print_REMD(*m_vel_traj, replica_data, sim.param());
    
    if(m_every_force && (sim.steps() % m_every_force) == 0){
      if (sim.steps())
	print_REMD(*m_force_traj, replica_data, sim.param());
    }
  } // decorated

}

void io::Out_Configuration::_print_jvalue(simulation::Parameter const & param,
					  configuration::Configuration const &conf,
					  topology::Topology const &topo,
					  std::ostream &os)
{
  DEBUG(10, "JVALUE Averages and LE data");

  if (param.jvalue.mode != simulation::restr_inst) {
    os << "JVALUERESEXPAVE\n";
    os.setf(std::ios::fixed, std::ios::floatfield);
    os.precision(7);
    std::vector<double>::const_iterator av_it = conf.special().jvalue_av.begin(),
            av_to = conf.special().jvalue_av.end();
    for( ;av_it != av_to; ++av_it){
      os << std::setw(15) << *av_it << "\n";
    }
    os << "END\n";
  }

  if (param.jvalue.le){
    os << "JVALUERESEPS\n";
    os.setf(std::ios::scientific, std::ios::floatfield); 
    os.precision(7);
    std::vector<std::vector<double> >::const_iterator
        le_it = conf.special().jvalue_epsilon.begin(),
        le_to = conf.special().jvalue_epsilon.end();
    
    for( ;le_it != le_to; ++le_it){ 
      for(unsigned int i=0; i < le_it->size(); ++i)
	os << std::setw(15) << (*le_it)[i];   
      os << "\n";
    }
    os << "END\n";
  }
}

void io::Out_Configuration::_print_distance_restraint_averages(
					  configuration::Configuration const &conf,
					  topology::Topology const &topo,
					  std::ostream &os)
{
  DEBUG(10, "distance restraint averages");
  
  std::vector<double>::const_iterator it = conf.special().distanceres_av.begin(),
          to = conf.special().distanceres_av.end();
  
  os.setf(std::ios::fixed, std::ios::floatfield);
  os.precision(m_distance_restraint_precision);
  
  os << "DISRESEXPAVE" << std::endl;
  int i;
  for(i = 1 ; it != to; ++it, ++i) {
    os << std::setw(m_width) << *it;
    if (i % 5 == 0) 
      os << std::endl;
  }
  if (i % 5 != 0) 
    os << std::endl;  
  
  os << "END" << std::endl;
}

void io::Out_Configuration::_print_pscale_jrest(configuration::Configuration const &conf,
						topology::Topology const &topo,
						std::ostream &os)
{
  DEBUG(10, "PSCALE JREST data");
  
  std::vector<topology::jvalue_restraint_struct>::const_iterator
    jval_it = topo.jvalue_restraints().begin(),
    jval_to = topo.jvalue_restraints().end();
  
  os << "PSCALEJREST\n";

  for(int i=0; jval_it != jval_to; ++jval_it, ++i){
    
    os << std::setw(15) << jval_it->i+1
       << std::setw(10) << jval_it->j+1
       << std::setw(10) << jval_it->k+1
       << std::setw(10) << jval_it->l+1
       << std::setw(10) << conf.special().pscale.scaling[i]
       << std::setw(15) << conf.special().pscale.t[i]
       << "\n";
  }
  os << "END\n";
}

void io::Out_Configuration::_print_replica_information
(
 std::vector<util::Replica_Data> const & replica_data,
 std::ostream &os
 )
{
  DEBUG(10, "replica information");
  
  std::vector<util::Replica_Data>::const_iterator
    it = replica_data.begin(),
    to = replica_data.end();
  
  os << "REPDATA\n";

  for(int i=0; it != to; ++it, ++i){
    
    os << std::setw(6) << it->ID
       << " "
       << std::setw(6) << it->run
       << std::setw(6) << it->Ti
       << std::setw(6) << it->li
       << std::setw(18) << it->epot_i
       << std::setw(6) << it->Tj
       << std::setw(6) << it->lj
       << std::setw(18) << it->epot_j
       << std::setw(18) << it->probability
       << std::setw(4) << it->switched
       << std::setw(18) << it->time
       << std::setw(4) << it->state
       << "\n";
  }
  os << "END\n";
}

static void _print_energyred_helper(std::ostream & os, configuration::Energy const &e)
{

  const int numenergygroups = unsigned(e.bond_energy.size());
  const int numbaths = unsigned(e.kinetic_energy.size());
  // const int energy_group_size = numenergygroups * (numenergygroups + 1) /2;

  os << "# totals\n";
  
  os << std::setw(18) << e.total << "\n"
     << std::setw(18) << e.kinetic_total << "\n"
     << std::setw(18) << e.potential_total << "\n"
     << std::setw(18) << e.bond_total << "\n"
     << std::setw(18) << e.angle_total << "\n"
     << std::setw(18) << e.improper_total << "\n"
     << std::setw(18) << e.dihedral_total << "\n"
     << std::setw(18) << e.lj_total << "\n"
     << std::setw(18) << e.crf_total << "\n"
     << std::setw(18) << e.ls_realspace_total << "\n"
     << std::setw(18) << e.ls_kspace_total << "\n"
     << std::setw(18) << e.ls_a_term_total << "\n"
     << std::setw(18) << e.ls_self_total << "\n"
     << std::setw(18) << e.ls_surface_total << "\n"
     << std::setw(18) << e.constraints_total << "\n"
     << std::setw(18) << e.posrest_total << "\n"
     << std::setw(18) << e.distanceres_total << "\n"
     << std::setw(18) << e.dihrest_total << "\n"
     << std::setw(18) << e.jvalue_total << "\n"
     << std::setw(18) << e.self_total << "\n" // self energy from polarization
     //<< std::setw(18) << 0.0 << "\n" // local elevation
     << std::setw(18) << e.eds_vr << "\n"; // eds energy of reference state
     //<< std::setw(18) << 0.0 << "\n"; // path integral
  // << std::setw(18) << e.entropy_term << "\n"; // dH/dl * H
  
  // put eds V_R energy here
  
  os << "# baths\n";
  os << numbaths << "\n";

  for(int i=0; i < numbaths; ++i){
    os << std::setw(18) << e.kinetic_energy[i]
       << std::setw(18) << e.com_kinetic_energy[i]
       << std::setw(18) << e.ir_kinetic_energy[i] << "\n";
  }

  os << "# bonded\n";
  os << numenergygroups << "\n";
  for(int i=0; i<numenergygroups; i++){
    os << std::setw(18) << e.bond_energy[i] 
       << std::setw(18) << e.angle_energy[i] 
       << std::setw(18) << e.improper_energy[i] 
       << std::setw(18) << e.dihedral_energy[i] << "\n";
  }

  os << "# nonbonded\n";
  for(int i=0; i<numenergygroups; i++){
    for(int j=i; j<numenergygroups; j++){
      os << std::setw(18) << e.lj_energy[j][i]
	 << std::setw(18) << e.crf_energy[j][i]
         << std::setw(18) << e.ls_real_energy[j][i] 
         << std::setw(18) << e.ls_k_energy[j][i] << "\n";
    }
  }

  os << "# special\n";
  for(int i=0; i<numenergygroups; i++){
    os << std::setw(18) << e.constraints_energy[i] 
       << std::setw(18) << e.posrest_energy[i] 
       << std::setw(18) << e.distanceres_energy[i] // disres
       << std::setw(18) << e.dihrest_energy[i] // dihedral res
       << std::setw(18) << 0.0 // jval
       << std::setw(18) << 0.0 // local elevation
       << std::setw(18) << 0.0 << "\n"; // path integral
  }
  
  // eds energy of end states
  os << "# eds\n";
  os << "# numstates\n";
  const unsigned int numstates = e.eds_vi.size();
  os << numstates << "\n";
  os << std::setw(18) << "# total" 
     << std::setw(18) << "nonbonded"
     << std::setw(18) << "special\n"; 
  for(unsigned i = 0; i < e.eds_vi.size(); i++){
    os << std::setw(18) << e.eds_vi[i] 
       << std::setw(18) << e.eds_vi[i] - e.eds_vi_special[i]
       << std::setw(18) << e.eds_vi_special[i] << "\n";
  }
  
  // write eds energies (vr,{V_i}) here
}

static void _print_volumepressurered_helper(std::ostream &os,
					    double mass,
					    simulation::Multibath const & m,
					    std::vector<double> const & s,
					    configuration::Energy const & e,
					    math::Box const & b,
					    math::boundary_enum t,
					    math::Matrix const & p,
					    math::Matrix const & v,
					    math::Matrix const & k)
{
  const int numbaths = int(m.size());

  os << "# mass\n";
  os << std::setw(18) << mass << "\n";
  
  os << "# temperature\n";
  os << numbaths << "\n";
  
  for(int i=0; i < numbaths; ++i){
    if (m[i].dof)
      os << std::setw(18) << 2 * e.kinetic_energy[i] / math::k_Boltzmann / m[i].dof;
    else
      os << std::setw(18) << 0.0;
    if (m[i].com_dof)
      os << std::setw(18) << 2 * e.com_kinetic_energy[i] / math::k_Boltzmann / m[i].com_dof;
    else
      os << std::setw(18) << 0.0;
    if (m[i].ir_dof)
      os << std::setw(18) << 2 * e.ir_kinetic_energy[i] / math::k_Boltzmann / m[i].ir_dof;
    else
      os << std::setw(18) << 0.0;
    
    if (s.size())
      os << std::setw(18) << s[i] << "\n";
    else
      os << std::setw(18) << m[i].scale << "\n";
    
  }

  os << "# volume\n";
  os << std::setw(18) << math::volume(b, t) << "\n";

  os << std::setw(18) << b(0)(0) << std::setw(18) << b(0)(1) << std::setw(18) << b(0)(2) << "\n"
     << std::setw(18) << b(1)(0) << std::setw(18) << b(1)(1) << std::setw(18) << b(1)(2) << "\n"
     << std::setw(18) << b(2)(0) << std::setw(18) << b(2)(1) << std::setw(18) << b(2)(2) << "\n";
  
  os << "# pressure\n";
  os << std::setw(18) << (p(0,0) + p(1,1) + p(2,2)) / 3.0 << "\n";
  os << std::setw(18) << (v(0,0) + v(1,1) + v(2,2)) / 3.0 << "\n";
  os << std::setw(18) << (k(0,0) + k(1,1) + k(2,2)) / 3.0 << "\n";

  os << std::setw(18) << p(0,0) << std::setw(18) << p(0,1) << std::setw(18) << p(0,2) << "\n"
     << std::setw(18) << p(1,0) << std::setw(18) << p(1,1) << std::setw(18) << p(1,2) << "\n"
     << std::setw(18) << p(2,0) << std::setw(18) << p(2,1) << std::setw(18) << p(2,2) << "\n";

  os << std::setw(18) << v(0,0) << std::setw(18) << v(0,1) << std::setw(18) << v(0,2) << "\n"
     << std::setw(18) << v(1,0) << std::setw(18) << v(1,1) << std::setw(18) << v(1,2) << "\n"
     << std::setw(18) << v(2,0) << std::setw(18) << v(2,1) << std::setw(18) << v(2,2) << "\n";

  os << std::setw(18) << k(0,0) << std::setw(18) << k(0,1) << std::setw(18) << k(0,2) << "\n"
     << std::setw(18) << k(1,0) << std::setw(18) << k(1,1) << std::setw(18) << k(1,2) << "\n"
     << std::setw(18) << k(2,0) << std::setw(18) << k(2,1) << std::setw(18) << k(2,2) << "\n";
  
}

void io::Out_Configuration
::_print_position_restraints(simulation::Simulation const & sim,
		             topology::Topology const &topo, 
		             std::ostream &os)
{
  os.setf(std::ios::fixed, std::ios::floatfield);
  os.precision(m_precision);
  
  os << "REFPOSITION\n";
  
  std::vector<topology::position_restraint_struct>::const_iterator it
          = topo.position_restraints().begin(),
          to = topo.position_restraints().end();
  topology::Solute const &solute = topo.solute();
  std::vector<std::string> const &residue_name = topo.residue_names();
  
  for(; it != to; ++it) {
    if (it->seq < topo.num_solute_atoms()) {
      os << std::setw(5)  << solute.atom(it->seq).residue_nr+1 << " "
         << std::setw(5)  << std::left << residue_name[solute.atom(it->seq).residue_nr] << " "
         << std::setw(6)  << std::left << solute.atom(it->seq).name << std::right
         << std::setw(6)  << it->seq+1
         << std::setw(m_width) << it->pos(0)
         << std::setw(m_width) << it->pos(1)
         << std::setw(m_width) << it->pos(2)
         << "\n";
    } else { // just writing out dummy values for first 17 chars
      os << std::setw(5)  << "0" << " "
	 << std::setw(5)  << std::left 
	 << "SOLV" << " "
	 << std::setw(6)  << std::left << "AT" << std::right
	 << std::setw(6)  << it->seq + 1
	 << std::setw(m_width) << it->pos(0)
	 << std::setw(m_width) << it->pos(1)
	 << std::setw(m_width) << it->pos(2)
	 << "\n";      
    }
  }
  
  os << "END\n";
  
  if (sim.param().posrest.posrest == simulation::posrest_bfactor) {
    it = topo.position_restraints().begin();
    os << "BFACTOR\n";
    
    for(; it != to; ++it) {
      if (it->seq < topo.num_solute_atoms()) {
        os << std::setw(5)  << solute.atom(it->seq).residue_nr+1 << " "
           << std::setw(5)  << std::left << residue_name[solute.atom(it->seq).residue_nr] << " "
           << std::setw(6)  << std::left << solute.atom(it->seq).name << std::right
           << std::setw(6)  << it->seq+1
           << std::setw(m_width) << it->bfactor
           << "\n";
      } else {
        os << std::setw(5)  << "0" << " "
           << std::setw(5)  << std::left
           << "SOLV" << " "
           << std::setw(6)  << std::left << "AT" << std::right
           << std::setw(6)  << it->seq + 1
           << std::setw(m_width) << it->bfactor
           << "\n";
      }
    }
    
    os << "END\n";
  }

}

void io::Out_Configuration::
_print_nose_hoover_chain_variables(const simulation::Multibath & multibath,
        std::ostream & os) {
  os.setf(std::ios::fixed, std::ios::floatfield);
  os.precision(m_precision);
  
  os << "NHCVARIABLES\n";
  for(simulation::Multibath::const_iterator it = multibath.begin(), to = multibath.end();
    it != to; ++it) {
    const std::vector<double> & zeta = it->zeta;
    for(std::vector<double>::const_iterator z_it = zeta.begin(), z_to = zeta.end();
      z_it != z_to; ++z_it) {
      os << std::setw(m_width) << *z_it;
    }
    os << "\n";
  }  
  os << "END\n";
}

void io::Out_Configuration::
_print_rottrans(configuration::Configuration const &conf,
        simulation::Simulation const &sim,
        std::ostream &os) {
  os.setf(std::ios::fixed, std::ios::floatfield);
  os.precision(m_precision);
  
  os << "ROTOTRANSREF\n";
  os << "# translational matrix\n";
  for(unsigned int i = 0; i < 3; ++i) {
    for(unsigned int j = 0; j < 3; ++j) {
      os << std::setw(m_width) << conf.special().rottrans_constr.theta_inv_trans(i,j);
    }
    os << "\n";
  }
  os << "# rotational matrix\n";
  for(unsigned int i = 0; i < 3; ++i) {
    for(unsigned int j = 0; j < 3; ++j) {
      os << std::setw(m_width) << conf.special().rottrans_constr.theta_inv_rot(i,j);
    }
    os << "\n";
  }
  os << "# reference positions\n";
  unsigned int last = sim.param().rottrans.last;
  const math::VArray & refpos = conf.special().rottrans_constr.pos;
  for(unsigned int i = 0; i < last; ++i ) {
    os << std::setw(m_width) << refpos(i)(0) 
       << std::setw(m_width) << refpos(i)(1) 
       << std::setw(m_width) << refpos(i)(2) << "\n";
  }
  os << "END\n";
}
