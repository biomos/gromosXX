/**
 * @file out_configuration.cc
 * definition of the Out_Configuration methods.
 */

#include <stdheader.h>

#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <simulation/simulation.h>
#include <configuration/configuration.h>

#include <math/periodicity.h>
#include <math/volume.h>

#include <io/print_block.h>
#include <io/argument.h>

#include <sstream>

#include "out_configuration.h"

#include <util/replica_data.h>
#include <util/template_split.h>
#include <util/debug.h>

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
  : m_output(os),
    m_final(false),
    m_replica(false),
    m_every_pos(0),
    m_every_vel(0),
    m_every_force(0),
    m_every_energy(0),
    m_every_free_energy(0),
    m_every_blockaverage(0),
    m_write_blockaverage_energy(false),
    m_write_blockaverage_free_energy(false),
    m_precision(9),
    m_force_precision(9),
    m_width(15),
    m_force_width(18),
    m_title(title)
{
  _print_title(m_title, "output file", os);
}

io::Out_Configuration::~Out_Configuration()
{
  // std::cout << "out_configuration destructor" << std::endl;
  
  if (m_every_pos){
    m_pos_traj.flush();
    m_pos_traj.close();
  }
  if (m_final){
    m_final_conf.flush();
    m_final_conf.close();
  }

  if (m_every_vel){
    m_vel_traj.flush();
    m_vel_traj.close();
  }
    
  if (m_every_force){
    m_force_traj.flush();
    m_force_traj.close();
  }
    
  if (m_every_energy){
    m_energy_traj.flush();
    m_energy_traj.close();
  }

  if (m_every_free_energy){
    m_free_energy_traj.flush();
    m_free_energy_traj.close();
  }

  if (m_write_blockaverage_energy){
	  m_blockaveraged_energy.flush();
	  m_blockaveraged_energy.close();
  }
  if (m_write_blockaverage_free_energy){
	  m_blockaveraged_free_energy.flush();
	  m_blockaveraged_free_energy.close();
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
  if (args.count("fin") > 0)
    final_configuration(args["fin"]);
  else io::messages.add("argument fin for final configuration required!",
			"Out_Configuration",
			io::message::error);

  if (args.count("trj") > 0)
    trajectory(args["trj"], param.write.position);
  else if (param.write.position)
    io::messages.add("write trajectory but no trj argument",
		     "Out_Configuration",
		     io::message::error);

  if (args.count("trv") > 0)
    velocity_trajectory(args["trv"], param.write.velocity);
  else if (param.write.velocity)
    io::messages.add("write velocity trajectory but no trv argument",
		     "Out_Configuration",
		     io::message::error);

  if (args.count("trf") > 0)
    force_trajectory(args["trf"], 1);

  if (args.count("tre") > 0)
    energy_trajectory(args["tre"], param.write.energy);
  else if (param.write.energy)
    io::messages.add("write energy trajectory but no tre argument",
		     "Out_Configuration",
		     io::message::error);

  if (args.count("trg") > 0)
    free_energy_trajectory(args["trg"], param.write.free_energy);
  else if (param.write.free_energy)
    io::messages.add("write free energy trajectory but no trg argument",
		     "Out_Configuration",
		     io::message::error);

  if (args.count("bae") > 0)
    block_averaged_energy(args["bae"], param.write.block_average);
  else if (param.write.block_average && param.write.energy)
    io::messages.add("write block averaged energy but no bae argument",
		     "Out_Configuration",
		     io::message::error);

  if (param.perturbation.perturbation){
    if (args.count("bag") > 0)
      block_averaged_free_energy(args["bag"], 
				 param.write.block_average);
    else if (param.write.block_average && param.write.free_energy)
      io::messages.add("write block averaged free energy "
			"but no bag argument",
		       "Out_Configuration",
		       io::message::error);
  }

  if (param.replica.num_T * param.replica.num_l){
    m_replica = true;
  }
  
}

void io::Out_Configuration::init
(
 std::string fin, std::string trj, std::string trv, 
 std::string trf, std::string tre, std::string trg,
 std::string bae, std::string bag,
 simulation::Parameter const & param
 )
{
  if (fin != "")
    final_configuration(fin);
  else io::messages.add("argument fin for final configuration required!",
			"Out_Configuration",
			io::message::error);

  if (trj != "")
    trajectory(trj, param.write.position);
  else if (param.write.position)
    io::messages.add("write trajectory but no trj argument",
		     "Out_Configuration",
		     io::message::error);

  if (trv != "")
    velocity_trajectory(trv, param.write.velocity);
  else if (param.write.velocity)
    io::messages.add("write velocity trajectory but no trv argument",
		     "Out_Configuration",
		     io::message::error);

  // use same step as position writing
  if (trf != "")
    force_trajectory(trf, param.write.position);

  if (tre != "")
    energy_trajectory(tre, param.write.energy);
  else if (param.write.energy)
    io::messages.add("write energy trajectory but no tre argument",
		     "Out_Configuration",
		     io::message::error);

  if (trg != "")
    free_energy_trajectory(trg, param.write.free_energy);
  else if (param.write.free_energy)
    io::messages.add("write free energy trajectory but no trg argument",
		     "Out_Configuration",
		     io::message::error);

  if (bae != "")
    block_averaged_energy(bae, param.write.block_average);
  else if (param.write.block_average && param.write.energy)
    io::messages.add("write block averaged energy but no bae argument",
		     "Out_Configuration",
		     io::message::error);

  if (param.perturbation.perturbation){
    if (bag != "")
      block_averaged_free_energy(bag, param.write.block_average);
    else if (param.write.block_average && param.write.free_energy)
      io::messages.add("write block averaged free energy "
			"but no bag argument",
		       "Out_Configuration",
		       io::message::error);
  }

}

void io::Out_Configuration::write(configuration::Configuration &conf,
				  topology::Topology const &topo,
				  simulation::Simulation const &sim,
				  output_format const form)
{
  // standard trajectories
  if (form == reduced){

    if(m_every_pos && (sim.steps() % m_every_pos) == 0){
      // don't write starting configuration if analyzing a trajectory
      if (sim.steps() || !sim.param().analyze.analyze){
	_print_timestep(sim, m_pos_traj);
	
	if (sim.param().write.solute_only)
	  _print_positionred(conf, topo,  topo.num_solute_atoms(), m_pos_traj);
	else
	  _print_positionred(conf, topo,  topo.num_atoms(), m_pos_traj);
	
	if (conf.boundary_type != math::vacuum)
	  _print_box(conf, m_pos_traj);
      }
    }
    
    if (m_every_vel && (sim.steps() % m_every_vel) == 0){
      _print_timestep(sim, m_vel_traj);
      _print_velocityred(conf, m_vel_traj);
    }
    
    if(m_every_force && ((sim.steps()) % m_every_force) == 0){
      if(sim.steps()){
	_print_old_timestep(sim, m_force_traj);
	_print_forcered(conf, m_force_traj);
      }
    }
    
    if(m_every_energy && (sim.steps() % m_every_energy) == 0){
      if(sim.steps()){
	_print_old_timestep(sim, m_energy_traj);
	_print_energyred(conf, m_energy_traj);
	_print_volumepressurered(topo, conf, sim, m_energy_traj);
      }
    }
    
    if(m_every_free_energy && (sim.steps() % m_every_free_energy) == 0){
      if(sim.steps()){
	_print_old_timestep(sim, m_free_energy_traj);
	_print_free_energyred(conf, topo, m_free_energy_traj);
      }
    }

    if (m_every_blockaverage && (sim.steps() % m_every_blockaverage) == 0){
      
      if(m_write_blockaverage_energy){
	if(sim.steps()){
	  _print_old_timestep(sim, m_blockaveraged_energy);
	  _print_blockaveraged_energyred(conf, m_blockaveraged_energy);
	  _print_blockaveraged_volumepressurered(conf, sim, m_blockaveraged_energy);
	}
      }

      if(m_write_blockaverage_free_energy){
	if(sim.steps()){
	  _print_old_timestep(sim, m_blockaveraged_free_energy);
	  _print_blockaveraged_free_energyred(conf, sim.param().perturbation.dlamt,
					      m_blockaveraged_free_energy);
	}
      }
      conf.current().averages.block().zero();
    }
    
  }
  else if(form == final && m_final){
    _print_timestep(sim, m_final_conf);
    _print_position(conf, topo, m_final_conf);

    if(sim.param().minimise.ntem == 0)
      _print_velocity(conf, topo, m_final_conf);

    _print_box(conf, m_final_conf);

    if(sim.param().constraint.solute.algorithm
       == simulation::constr_flexshake){
      _print_flexv(conf, topo, m_final_conf);
    }

    if(sim.param().jvalue.mode != simulation::restr_off){
      _print_jvalue(conf, topo, m_final_conf);
    }

    if(sim.param().pscale.jrest){
      _print_pscale_jrest(conf, topo, m_final_conf);
    }
    
    // forces and energies still go to their trajectories
    if (m_every_force && ((sim.steps()) % m_every_force) == 0){
      _print_old_timestep(sim, m_force_traj);
      _print_forcered(conf, m_force_traj);
    }

    if(m_every_energy && (sim.steps() % m_every_energy) == 0){
      _print_old_timestep(sim, m_energy_traj);
      _print_energyred(conf, m_energy_traj);
      _print_volumepressurered(topo, conf, sim, m_energy_traj);
    }

    if(m_every_free_energy && (sim.steps() % m_every_free_energy) == 0){
      _print_old_timestep(sim, m_free_energy_traj);
      _print_free_energyred(conf, topo, m_free_energy_traj);
    }

  }
  else{

    // not reduced or final (so: decorated)

    if(m_every_pos && (sim.steps() % m_every_pos) == 0){
      _print_timestep(sim, m_pos_traj);
      _print_position(conf, topo, m_pos_traj);
      if (conf.boundary_type != math::vacuum)
	_print_box(conf, m_pos_traj);
    }
    
    if (m_every_vel && (sim.steps() % m_every_vel) == 0){
      _print_timestep(sim, m_vel_traj);
      _print_velocity(conf, topo, m_vel_traj);
    }
    
    if(m_every_force && (sim.steps() % m_every_force) == 0){
      if (sim.steps()){
	_print_timestep(sim, m_force_traj);	
	_print_force(conf, topo, m_force_traj);
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
      _print_timestep(sim, m_pos_traj);
      _print_replica_information(replica_data, m_pos_traj);
      
      for(unsigned int i=0; i<conf.size(); ++i){
	_print_positionred(conf[i], topo,  topo.num_atoms(), m_pos_traj);
	_print_velocityred(conf[0], m_vel_traj);
	
	if (conf[i].boundary_type != math::vacuum)
	  _print_box(conf[i], m_pos_traj);
      }
    }
  }
  else if(form == final && m_final){
    for(unsigned int i=0; i<conf.size(); ++i){
      
      m_final_conf << "REPLICAFRAME\n"
		   << std::setw(12) << i + 1
		   << "\nEND\n";

      if (i==0){
	_print_timestep(sim, m_final_conf);
	_print_replica_information(replica_data, m_final_conf);
      }
      
      _print_position(conf[i], topo, m_final_conf);
      _print_velocity(conf[i], topo, m_final_conf);
      _print_box(conf[i], m_final_conf);
    }
    
    /*
    if(sim.param().jvalue.mode != simulation::restr_off){
      _print_jvalue(conf[0], topo, m_final_conf);
    }
    if(sim.param().pscale.jrest){
      _print_pscale_jrest(conf[0], topo, m_final_conf);
    }
    */
  }
  // done writing replicas!
}

void io::Out_Configuration
::final_configuration(std::string name)
{
  m_final_conf.open(name.c_str());
  _print_title(m_title, "final configuration", m_final_conf);
  m_final = true;
}

void io::Out_Configuration
::trajectory(std::string name, int every)
{
  m_pos_traj.open(name.c_str());
  m_every_pos = every;
  _print_title(m_title, "position trajectory", m_pos_traj);
}

void io::Out_Configuration
::velocity_trajectory(std::string name, int every)
{
  m_vel_traj.open(name.c_str());
  m_every_vel = every;
  _print_title(m_title, "velocity trajectory", m_vel_traj);
}

void io::Out_Configuration
::force_trajectory(std::string name, int every)
{
  m_force_traj.open(name.c_str());
  m_every_force = every;
  _print_title(m_title, "force trajectory", m_force_traj);
}

void io::Out_Configuration
::energy_trajectory(std::string name, int every)
{
  m_energy_traj.open(name.c_str());
  m_every_energy = every;
  _print_title(m_title, "energy trajectory", m_energy_traj);
}

void io::Out_Configuration
::free_energy_trajectory(std::string name, int every)
{
  m_free_energy_traj.open(name.c_str());
  m_every_free_energy = every;
  _print_title(m_title, "free energy trajectory", m_free_energy_traj);
}

void io::Out_Configuration
::block_averaged_energy(std::string name, int every)
{
  m_blockaveraged_energy.open(name.c_str());
  if (m_every_blockaverage && m_every_blockaverage != every){
    io::messages.add("overwriting how often block averages are written out illegal",
		     "Out_Configuration",
		     io::message::error);
  }
  m_every_blockaverage = every;
  m_write_blockaverage_energy = true;
  _print_title(m_title, "block averaged energies", m_blockaveraged_energy);
}

void io::Out_Configuration
::block_averaged_free_energy(std::string name, int every)
{
  m_blockaveraged_free_energy.open(name.c_str());
  if (m_every_blockaverage && m_every_blockaverage != every){
    io::messages.add("overwriting how often block averages are written out illegal",
		     "Out_Configuration",
		     io::message::error);
  }
  m_every_blockaverage = every;
  m_write_blockaverage_free_energy = true;
  _print_title(m_title, "block averaged free energies", m_blockaveraged_free_energy);
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
::_print_velocityred(configuration::Configuration const &conf,
		     std::ostream &os)
{
  os.setf(std::ios::fixed, std::ios::floatfield);
  os.precision(m_precision);
  
  os << "VELOCITYRED\n";
  
  math::VArray const &vel = conf.current().vel;
  
  for(int i=0,to = vel.size(); i<to; ++i){

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
	       std::ostream &os)
{
  os.setf(std::ios::fixed, std::ios::floatfield);
  os.precision(m_force_precision);
  
  os << "FORCE\n";
  
  math::VArray const &force = conf.current().force;
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

}

void io::Out_Configuration
::_print_forcered(configuration::Configuration const &conf,
		  std::ostream &os)
{
  os.setf(std::ios::fixed, std::ios::floatfield);
  os.precision(m_force_precision);
  
  os << "FORCERED\n";
  
  math::VArray const & force = conf.current().force;
  
  for(int i=0,to = force.size(); i<to; ++i){
    
    os << std::setw(m_force_width) << force(i)(0)
       << std::setw(m_force_width) << force(i)(1)
       << std::setw(m_force_width) << force(i)(2)
       << "\n";

    if((i+1)% 10 == 0) os << '#' << std::setw(10) << i+1 << "\n";
  }
  
  os << "END\n";
  
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

void io::Out_Configuration
::_print_box(configuration::Configuration const &conf,
	     std::ostream &os)
{
  os.setf(std::ios::fixed, std::ios::floatfield);
  os.precision(m_precision);
  
  os << "TRICLINICBOX\n";
  
  math::Box const &box = conf.current().box;
  
  os << std::setw(5) << conf.boundary_type << "\n";
  
  for(int i=0,to = 3; i<to; ++i){
    
    os << std::setw(m_width) << box(0)(i)
       << std::setw(m_width) << box(1)(i)
       << std::setw(m_width) << box(2)(i)
       << "\n";
  }
  
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
  if ((sim.steps() % sim.param().print.stepblock) == 0){

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
      print_REMD(m_pos_traj, replica_data, sim.param());
    
    if (m_every_vel && (sim.steps() % m_every_vel) == 0)
      print_REMD(m_vel_traj, replica_data, sim.param());
    
    if(m_every_force && ((sim.steps()) % m_every_force) == 0){
      if(sim.steps())
	print_REMD(m_force_traj, replica_data, sim.param());
    }
    
    if(m_every_energy && (sim.steps() % m_every_energy) == 0){
      if(sim.steps())
	print_REMD(m_energy_traj, replica_data, sim.param());
    }
    
    if(m_every_free_energy && (sim.steps() % m_every_free_energy) == 0){
      if(sim.steps())
	print_REMD(m_free_energy_traj, replica_data, sim.param());
    }

    if (m_every_blockaverage && (sim.steps() % m_every_blockaverage) == 0){
      if(m_write_blockaverage_energy){
	if(sim.steps())
	  print_REMD(m_blockaveraged_energy, replica_data, sim.param());
      }

      if(m_write_blockaverage_free_energy){
	if(sim.steps())
	  print_REMD(m_blockaveraged_free_energy, replica_data, sim.param());
      }
    }
  } // reduced

  else if(form == final && m_final){
    print_REMD(m_final_conf, replica_data, sim.param());
    
    // forces and energies still go to their trajectories
    if (m_every_force && ((sim.steps()) % m_every_force) == 0)
      print_REMD(m_force_traj, replica_data, sim.param());

    if(m_every_energy && (sim.steps() % m_every_energy) == 0)
      print_REMD(m_energy_traj, replica_data, sim.param());

    if(m_every_free_energy && (sim.steps() % m_every_free_energy) == 0)
      print_REMD(m_free_energy_traj, replica_data, sim.param());
  } // final

  else{

    // not reduced or final (so: decorated)

    if(m_every_pos && (sim.steps() % m_every_pos) == 0)
      print_REMD(m_pos_traj, replica_data, sim.param());
    
    if (m_every_vel && (sim.steps() % m_every_vel) == 0)
      print_REMD(m_vel_traj, replica_data, sim.param());
    
    if(m_every_force && (sim.steps() % m_every_force) == 0){
      if (sim.steps())
	print_REMD(m_force_traj, replica_data, sim.param());
    }
  } // decorated

}

void io::Out_Configuration::_print_jvalue(configuration::Configuration const &conf,
					  topology::Topology const &topo,
					  std::ostream &os)
{
  DEBUG(10, "JVALUE Aversages");
  
  std::vector<double>::const_iterator av_it = conf.special().jvalue_av.begin();
  std::vector<topology::jvalue_restraint_struct>::const_iterator
    jval_it = topo.jvalue_restraints().begin(),
    jval_to = topo.jvalue_restraints().end();
  
  os << "JVALRESEXPAVE03\n";

  for( ; jval_it != jval_to; ++jval_it, ++av_it){
    
    os << std::setw(15) << jval_it->i+1
       << std::setw(10) << jval_it->j+1
       << std::setw(10) << jval_it->k+1
       << std::setw(10) << jval_it->l+1
       << std::setw(20) << *av_it
       << "\n";
  }
  os << "END\n";

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
       << std::setw(6) << it->run
       << std::setw(6) << it->Ti
       << std::setw(6) << it->li
       << std::setw(18) << it->epot_i
       << std::setw(6) << it->Tj
       << std::setw(6) << it->lj
       << std::setw(18) << it->epot_j
       << std::setw(18) << it->probability
       << std::setw(4) << it->switched
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
     << std::setw(18) << e.constraints_total << "\n"
     << std::setw(18) << e.posrest_total << "\n"
     << std::setw(18) << e.distrest_total << "\n"
     << std::setw(18) << 0.0 << "\n" // dihedral res
     << std::setw(18) << e.jvalue_total << "\n"
     << std::setw(18) << 0.0 << "\n" // local elevation
     << std::setw(18) << 0.0 << "\n"; // path integral
  
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
	 << std::setw(18) << e.crf_energy[j][i] << "\n";

    }
  }

  os << "# special\n";
  for(int i=0; i<numenergygroups; i++){
    os << std::setw(18) << e.constraints_energy[i] 
       << std::setw(18) << e.posrest_energy[i] 
       << std::setw(18) << e.distrest_energy[i] // disres
       << std::setw(18) << 0.0 // dihedral res
       << std::setw(18) << 0.0 // jval
       << std::setw(18) << 0.0 // local elevation
       << std::setw(18) << 0.0 << "\n"; // path integral
  }
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
