/**
 * @file in_configuration.cc
 * implements methods of In_Configuration.
 */

#include <stdheader.h>

#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <simulation/multibath.h>
#include <simulation/parameter.h>
#include <simulation/simulation.h>
#include <configuration/configuration.h>

#include <io/blockinput.h>
#include <io/instream.h>
#include <io/configuration/inframe.h>

#include <util/generate_velocities.h>
#include <util/replica_data.h>

#include <math/volume.h>

#include "in_configuration.h"

#undef MODULE
#undef SUBMODULE
#define MODULE io
#define SUBMODULE configuration

static std::set<std::string> block_read;

/**
 * read in a trajectory.
 */
void io::In_Configuration::read(configuration::Configuration &conf, 
				topology::Topology &topo, 
				simulation::Simulation & sim,
				std::ostream & os)
{
  if (!quiet)
    os << "\nCONFIGURATION\n";

  simulation::Parameter const & param = sim.param();
  
  // resize the configuration
  conf.resize(topo.num_atoms());

  DEBUG(8, "reading in a frame");
  read_frame();

  block_read.clear();

  read_time(topo, conf, sim, os);
  read_position(topo, conf, sim, os);
  read_velocity(topo, conf, sim, os);
  read_box(topo, conf, sim, os);
  read_jvalue(topo, conf, sim, os);
  read_pscale(topo, conf, sim, os);
  read_flexv(topo, conf, sim, os);
  
  // and set the boundary type!
  conf.boundary_type = param.boundary.boundary;


  // print some information
  if (!quiet){
    os << "\n\t";
    switch(conf.boundary_type){
      case math::vacuum:
	os << "PBC            = vacuum\n";
	break;
      case math::rectangular:
	os << "PBC            = rectangular\n";
	break;
      case math::triclinic:
	os << "PBC            = triclinic\n";
	break;
      case math::truncoct:
	os << "PBC            = truncoct\n";
	break;
      default:
	os << "wrong periodic boundary conditions!";
	io::messages.add("wrong PBC!", "In_Configuration", io::message::error);
    }
  }

  if (!quiet){
    os << "\ttotal mass     = " << math::sum(topo.mass()) << "\n"
       << "\tvolume         = " << math::volume(conf.current().box, conf.boundary_type);
    
    if (conf.boundary_type != math::vacuum)
      os << "\n\tdensity        = " 
	 << math::sum(topo.mass()) / math::volume(conf.current().box,
							 conf.boundary_type);
    
    os << "\n\n";
  }
  

  // warn for unread input data
  for(std::map<std::string, std::vector<std::string> >::const_iterator
	it = m_block.begin(),
	to = m_block.end();
      it != to;
      ++it){

    if (block_read.count(it->first) == 0 && it->second.size()){
      io::messages.add("block " + it->first + " not read in!",
		       "In_Configuration",
		       io::message::warning);
    }
  }
  
  if (!quiet)
    os << "\n\nEND\n\n";
}

/**
 * read in a trajectory.
 */
void io::In_Configuration::read_replica
(
 std::vector<configuration::Configuration> & conf, 
 topology::Topology &topo, 
 simulation::Simulation & sim,
 std::vector<util::Replica_Data> & replica_data,
 std::ostream & os)
{
  if (!quiet)
    os << "\nREPLICA CONFIGURATION\n";
  
  simulation::Parameter const & param = sim.param();

  const int switch_T = sim.param().replica.num_T;
  const int switch_l = sim.param().replica.num_l;

  const int rep_num = switch_T * switch_l;
  conf.resize(rep_num);
  replica_data.resize(rep_num);

  block_read.clear();
  
  DEBUG(8, "reading in a frame");
  read_frame();

  // is it a multi-configuration (replica) file?
  std::vector<std::string> buffer = m_block["REPLICAFRAME"];
  
  if(!buffer.size()){
    conf[0].resize(topo.num_atoms());

    read_time(topo, conf[0], sim, os);
    read_position(topo, conf[0], sim, os);
    read_velocity(topo, conf[0], sim, os);
    read_box(topo, conf[0], sim, os);
    read_jvalue(topo, conf[0], sim, os);
    read_pscale(topo, conf[0], sim, os);
    read_flexv(topo, conf[0], sim, os);
  
    conf[0].boundary_type = param.boundary.boundary;

    for(unsigned int i=1; i<conf.size(); ++i)
      conf[i] = conf[0];

    // setup replica information
    {
      int i=0;
      for(int l=0; l < switch_l; ++l){
	for(int t=0; t < switch_T; ++t, ++i){
	  
	  replica_data[i].ID = i;
	  replica_data[i].run = 0;
	  
	  replica_data[i].Ti = t;
	  replica_data[i].Tj = t;
	  
	  replica_data[i].li = l;
	  replica_data[i].lj = l;
	  
	  replica_data[i].epot_i = 0.0;
	  replica_data[i].epot_j = 0.0;
	  
	  replica_data[i].state = util::waiting;
	  replica_data[i].probability = 0.0;
	  replica_data[i].switched = false;
	}
      }
    }
  }
  else{
    read_replica_information(replica_data, os);

    for(unsigned int i=0; i<conf.size(); ++i){
      conf[i].resize(topo.num_atoms());
      
      read_time(topo, conf[i], sim, os);
      read_position(topo, conf[i], sim, os);
      read_velocity(topo, conf[i], sim, os);
      read_box(topo, conf[i], sim, os);
      read_jvalue(topo, conf[i], sim, os);
      read_pscale(topo, conf[i], sim, os);
      read_flexv(topo, conf[i], sim, os);
      
      conf[i].boundary_type = param.boundary.boundary;
      
      // warn for unread input data
      for(std::map<std::string, std::vector<std::string> >::const_iterator
	    it = m_block.begin(),
	    to = m_block.end();
	  it != to;
	  ++it){
	
	if (block_read.count(it->first) == 0 && it->second.size()){
	  io::messages.add("block " + it->first + " not read in!",
			   "In_Configuration",
			   io::message::warning);
	}
      }
      
      block_read.clear();
      read_frame();
    }
  }

  // print some information
  if (!quiet){
    os << "\n\t";
    switch(conf[0].boundary_type){
      case math::vacuum:
	os << "PBC            = vacuum\n";
	break;
      case math::rectangular:
	os << "PBC            = rectangular\n";
	break;
      case math::triclinic:
	os << "PBC            = triclinic\n";
	break;
      case math::truncoct:
	os << "PBC            = truncoct\n";
	break;
      default:
	os << "wrong periodic boundary conditions!";
	io::messages.add("wrong PBC!", "In_Configuration", io::message::error);
    }
  }

  if (!quiet){
    os << "\ttotal mass     = " << math::sum(topo.mass()) << "\n";
    
    if (conf[0].boundary_type != math::vacuum){
      const double v = math::volume(conf[0].current().box, 
				    conf[0].boundary_type);
      os << "\tvolume         = " 
	 << std::setw(18) << v
	 << "\n\tdensity        = " 
	 << std::setw(18) << math::sum(topo.mass()) / v;
    }
    
    os << "\n\n";
  }

  // warn for unread input data
  for(std::map<std::string, std::vector<std::string> >::const_iterator
	it = m_block.begin(),
	to = m_block.end();
      it != to;
      ++it){

    if (block_read.count(it->first) == 0 && it->second.size()){
      io::messages.add("block " + it->first + " not read in!",
		       "In_Configuration",
		       io::message::warning);
    }
  }
  
  if (!quiet)
    os << "\n\nEND\n\n";
}

/**
 * read in a trajectory.
 */
bool io::In_Configuration::read_next
(
 topology::Topology &topo,
 configuration::Configuration &conf,
 simulation::Simulation & sim,
 std::ostream & os
 )
{
  if (!quiet)
    os << "\nread next frame\n";

  read_frame();

  block_read.clear();

  // ignore errors reading time step
  read_time_step(topo, conf, sim, os);
  
  if (!(
	read_position(topo, conf, sim, os) &&
	read_box(topo, conf, sim, os)

	// read_velocity(topo, conf, sim, os);
	// read_jvalue(topo, conf, sim, os);
	// read_pscale(topo, conf, sim, os);
	// read_flexv(topo, conf, sim, os);

	)){
    std::cout << "could not read frame!!!" << std::endl;
    return false;
  }

  // print some information
  if (!quiet){
    os << "\n\t";
    os << "time : " << sim.time()
       << "\n\t"
       << "step : " << sim.steps();

    os << "\n\n";
  }

  // warn for unread input data
  for(std::map<std::string, std::vector<std::string> >::const_iterator
	it = m_block.begin(),
	to = m_block.end();
      it != to;
      ++it){

    if (block_read.count(it->first) == 0 && it->second.size()){
      os << "\tblock " + it->first + " not read in!\n";
    }
  }

  return true;
}

bool io::In_Configuration::read_position
(
 topology::Topology &topo, 
 configuration::Configuration &conf, 
 simulation::Simulation & sim,
 std::ostream & os)
{
  // read positions
  std::vector<std::string> buffer;
  buffer = m_block["POSITION"];
  if (buffer.size()){
    
    check_coordinates(topo, conf, sim, buffer.size() - 1, os);
    
    if (!quiet)
      os << "\treading POSITION...\n";
    _read_position(conf.current().pos, buffer, topo.num_atoms());
    block_read.insert("POSITION");
  }
  else{
    buffer = m_block["POSITIONRED"];
    if (buffer.size()){

      check_coordinates(topo, conf, sim, buffer.size() - 1, os);

      if (!quiet)
	os << "\treading POSITIONRED...\n";
      _read_positionred(conf.current().pos, buffer, topo.num_atoms());
      block_read.insert("POSITIONRED");
    }
    else{
      io::messages.add("no POSITION / POSITIONRED block found in input configuration",
		       "in_configuration",
		       io::message::error);
      return false;
    }
  }
  return true;
}

bool io::In_Configuration::read_velocity
(
 topology::Topology &topo, 
 configuration::Configuration &conf, 
 simulation::Simulation & sim,
 std::ostream & os)
{
  // read velocities
  std::vector<std::string> buffer;
  if(!sim.param().start.generate_velocities && !sim.param().minimise.ntem){
    buffer = m_block["VELOCITY"];
    if (buffer.size()){
      if (!quiet)
	os << "\treading VELOCITY...\n";
      _read_velocity(conf.current().vel, buffer, topo.num_atoms());
      block_read.insert("VELOCITY");
    }
    else{
      buffer = m_block["VELOCITYRED"];
      if (buffer.size()){
	if (!quiet)
	  os << "\treading VELOCITYRED...\n";
	_read_velocityred(conf.current().vel, buffer, topo.num_atoms());
	block_read.insert("VELOCITYRED");
      }
      else{
	io::messages.add("no VELOCITY / VELOCITYRED block found in input configuration",
			 "in_configuration",
			 io::message::error);
	conf.current().vel = 0.0;
	return false;
      }
    }
    // store also in old velocities (for initial temperature calculation)
    conf.old().vel = conf.current().vel;
  }
  else{
    // generate initial velocities
    util::generate_velocities(sim.param().start.tempi, 
			      topo.mass(),
			      conf.current().vel,
			      conf.old().vel,
			      sim.param().start.ig,
			      os);
  }
  return true;
}

bool io::In_Configuration::read_box
(
 topology::Topology &topo, 
 configuration::Configuration &conf, 
 simulation::Simulation & sim,
 std::ostream & os)
{
  // read box
  std::vector<std::string> buffer;
  if(sim.param().boundary.boundary != math::vacuum){
    buffer = m_block["TRICLINICBOX"];
    if (buffer.size()){
      if (!quiet)
	os << "\treading TRICLINICBOX...\n";
      _read_box(conf.current().box, buffer, sim.param().boundary.boundary);
      conf.old().box = conf.current().box;
      block_read.insert("TRICLINICBOX");
    }
    else{
      buffer = m_block["BOX"];
      if (buffer.size() && (sim.param().boundary.boundary == math::rectangular ||
			    sim.param().boundary.boundary == math::truncoct)){
	if (!quiet)
	  os << "\treading BOX...\n";
	_read_g96_box(conf.current().box, buffer);
	conf.old().box = conf.current().box;
	block_read.insert("BOX");
      }
      else{
	io::messages.add("no TRICLINICBOX / BOX (for rectangular/truncoct "
			 "boundary conditions)\n"
			 "\tblock found in input configuration",
			 "in_configuration",
			 io::message::error);
	return false;
      }
    }    
  }
  return true;
}

bool io::In_Configuration::read_jvalue
(
 topology::Topology &topo, 
 configuration::Configuration &conf, 
 simulation::Simulation & sim,
 std::ostream & os)
{
  std::vector<std::string> buffer;

  if (sim.param().jvalue.mode != simulation::restr_off){

    if (sim.param().jvalue.mode == simulation::restr_inst){
      if (sim.param().jvalue.read_av)
	io::messages.add("instantaneous J-value restraints, ignoring reading of averages",
			 "in_configuration",
			 io::message::warning);
    }
    else if (!sim.param().jvalue.read_av && sim.param().jvalue.mode != simulation::restr_inst){

      buffer = m_block["JVALRESEXPAVE03"];
      if (buffer.size()){
	block_read.insert("JVALRESEXPAVE03");

	io::messages.add("re-initialising J-restraint averages, non-continuous simulation",
			 "in_configuration",
			 io::message::warning);
      }
      else{
	io::messages.add("initialising J-restraint averages",
			 "in_configuration",
			 io::message::notice);
      }
    }
    else {

      buffer = m_block["JVALRESEXPAVE03"];
      if (buffer.size())
      {
	block_read.insert("JVALRESEXPAVE03");
	_read_jvalue_av(buffer, conf.special().jvalue_av, topo.jvalue_restraints());
      }
      else{
	io::messages.add("reading in of J-restraints averages requested "
			 "but JVALRESEXPAVE03 block not found",
			 "in_configuration",
			 io::message::error);
	return false;
      }
    }
  } // jvalue averages
  return true;
}

bool io::In_Configuration::read_pscale
(
 topology::Topology &topo, 
 configuration::Configuration &conf, 
 simulation::Simulation & sim,
 std::ostream & os)
{
  std::vector<std::string> buffer;
  if (sim.param().pscale.jrest){

    if (!sim.param().pscale.read_data){
      buffer = m_block["PSCALEJREST"];
      if (buffer.size()){
	block_read.insert("PSCALEJREST");
	io::messages.add("re-initialising J-restraints periodic scaling data, non-continuous simulation",
			 "in_configuration",
			 io::message::warning);
      }
      else{
	io::messages.add("initialising J-restraints periodic scaling data",
			 "in_configuration",
			 io::message::notice);
      }
    }
    else {
      buffer = m_block["PSCALEJREST"];
      if (buffer.size())
      {
	block_read.insert("PSCALEJREST");
	_read_pscale_jrest(buffer, conf.special().pscale, topo.jvalue_restraints());
      }
      else{
	io::messages.add("reading in of J-restraints periodic scaling data requested "
			 "but PSCALEJREST block not found",
			 "in_configuration",
			 io::message::error);
	return false;
      }
    }
  } // PSCALE JREST
  return true;
}

bool io::In_Configuration::read_flexv
(
 topology::Topology &topo, 
 configuration::Configuration &conf, 
 simulation::Simulation & sim,
 std::ostream & os)
{
  std::vector<std::string> buffer;
  if (sim.param().constraint.solute.algorithm == simulation::constr_flexshake){
    conf.special().flexible_constraint.flexible_vel.
      resize(topo.solute().distance_constraints().size()+
	     topo.perturbed_solute().distance_constraints().size());
    const unsigned int numb = unsigned(sim.param().multibath.multibath.size());

    conf.special().flexible_constraint.flexible_ekin.resize(numb);

    buffer = m_block["FLEXV"];
    if (buffer.size() && sim.param().constraint.solute.flexshake_readin){
      block_read.insert("FLEXV");
      if (!quiet)
	os << "\treading FLEXV...\n";
      _read_flexv(conf.special().flexible_constraint.flexible_vel, buffer, 
		  topo.solute().distance_constraints(),
		  topo.perturbed_solute().distance_constraints());

    }
    else{
      if (sim.param().constraint.solute.flexshake_readin){
	io::messages.add("no FLEXV block found but reading in of constraint velocities requested",
			 "in_configuration",
			 io::message::error);
	return false;
      }
      else
	io::messages.add("no FLEXV block found, assuming SHAKEd positions (and velocities)",
			 "in_configuration",
			 io::message::notice);
    }
  }
  return true;
}

bool io::In_Configuration::read_time
(
 topology::Topology &topo, 
 configuration::Configuration &conf, 
 simulation::Simulation & sim,
 std::ostream & os)
{
  // current time
  std::vector<std::string> buffer;
  if (sim.param().step.t0 == -1){
    buffer = m_block["TIMESTEP"];
    if (!buffer.size()){
      io::messages.add("Requested time information from TIMESTEP block, "
		       "but block not found", "in_configuration", io::message::error);
      return false;
    }
    else{
      _read_time(buffer, sim.time());
      block_read.insert("TIMESTEP");
    }
  }
  else{
    buffer = m_block["TIMESTEP"];
    if (buffer.size())
      block_read.insert("TIMESTEP");
  }
  return true;
}

bool io::In_Configuration::read_time_step
(
 topology::Topology &topo, 
 configuration::Configuration &conf, 
 simulation::Simulation & sim,
 std::ostream & os)
{
  // current time
  std::vector<std::string> buffer;

  buffer = m_block["TIMESTEP"];
  if (!buffer.size()){
    io::messages.add("TIMESTEP block not found",
		     "in_configuration",
		     io::message::error);
    return false;
  }
  else{
    _read_time_step(buffer, sim);
    block_read.insert("TIMESTEP");
  }
  return true;
}

bool io::In_Configuration::_read_positionred(math::VArray &pos, 
					     std::vector<std::string> &buffer,
					     int const num)
{
  DEBUG(8, "read positionred");
  
  // no title in buffer!
  std::vector<std::string>::const_iterator it = buffer.begin(),
    to = buffer.end()-1;
  
  int i;

  if (pos.size() < unsigned(num)){
    io::messages.add("configuration: too many coordinates for given topology",
		     "in_configuration",
		     io::message::critical);
    std::cout << "position size is : " << pos.size() << " and num is : " << num << std::endl;
    return false;
  }

  for(i=0; it != to; ++i, ++it){
    if (i >= num){
      io::messages.add("configuration file does not match topology: "
		       "too many coordinates in POSITIONRED block",
		       "In_Configuration",
		       io::message::error);
      break;
    }
    
    _lineStream.clear();
    _lineStream.str(*it);
    _lineStream >> pos(i)(0) >> pos(i)(1) >> pos(i)(2);
    
    if(_lineStream.fail()){
      io::messages.add("bad line in POSITIONRED block",
		       "In_Configuration",
		       io::message::error);      
      return false;
    }
  }

  if (i != num){
    io::messages.add("configuration file does not match topology: "
		     "not enough coordinates in POSITIONRED block",
		     "In_Configuration",
		     io::message::error);
    return false;
  }

  return true;
  
}

bool io::In_Configuration::_read_position(math::VArray &pos, std::vector<std::string> &buffer,
					  int const num)
{
  DEBUG(8, "read position");

  // no title in buffer!
  std::vector<std::string>::const_iterator it = buffer.begin(),
    to = buffer.end()-1;
  
  std::string s1, s2;
  int i, n, nr;

  for(i=0; it != to; ++i, ++it){
    if (i >= num){
      io::messages.add("configuration file does not match topology: "
		       "too many coordinates in POSITION block",
		       "In_Configuration",
		       io::message::error);
      break;
    }

    _lineStream.clear();
    _lineStream.str(*it);
    // ignore first 4 fields
    _lineStream >> n >> s1 >> s2 >> nr;
    _lineStream >> pos(i)(0) >> pos(i)(1) >> pos(i)(2);
    
    if(_lineStream.fail()){
      io::messages.add("bad line in POSITION block",
		       "In_Configuration",
		       io::message::critical);
      return false;
    }
  }

  if (i != num){
    io::messages.add("configuration file does not match topology: "
		     "not enough coordinates in POSITION block",
		     "In_Configuration",
		     io::message::error);
    return false;
  }

  return true;
  
}

bool io::In_Configuration::_read_velocityred(math::VArray &vel, 
					     std::vector<std::string> &buffer,
					     int const num)
{
  DEBUG(8, "read velocityred");

  // no title in buffer!
  std::vector<std::string>::const_iterator it = buffer.begin(),
    to = buffer.end()-1;
  
  int i;
  
  for(i=0; it != to; ++i, ++it){
    if (i >= num){
      io::messages.add("configuration file does not match topology: "
		       "too many coordinates in VELOCITYRED block",
		       "In_Configuration",
		       io::message::error);
      break;
    }

    _lineStream.clear();
    _lineStream.str(*it);
    _lineStream >> vel(i)(0) >> vel(i)(1) >> vel(i)(2);
    
    if(_lineStream.fail()){
      io::messages.add("bad line in VELOCITYRED block",
		       "In_Configuration",
		       io::message::error);      
      return false;
    }

  }

  if (i != num){
    io::messages.add("configuration file does not match topology: "
		     "not enough coordinates in VELOCITYRED block",
		     "In_Configuration",
		     io::message::error);
    return false;
  }
  
  return true;
  
}

bool io::In_Configuration::_read_velocity(math::VArray &vel, 
					  std::vector<std::string> &buffer,
					  int const num)
{
  DEBUG(8, "read velocity");

  // no title in buffer!
  std::vector<std::string>::const_iterator it = buffer.begin(),
    to = buffer.end()-1;
  
  std::string s1, s2;
  int i, n, nr;

  for(i=0; it != to; ++i, ++it){
    if (i >= num){
      io::messages.add("configuration file does not match topology: "
		       "too many coordinates in VELOCITY block",
		       "In_Configuration",
		       io::message::error);
      break;
    }
   
    _lineStream.clear();
    _lineStream.str(*it);
    // ignore first 4 fields
    _lineStream >> n >> s1 >> s2 >> nr;
    _lineStream >> vel(i)(0) >> vel(i)(1) >> vel(i)(2);
    
    if(_lineStream.fail()){
      io::messages.add("bad line in VELOCITY block",
		       "In_Configuration",
		       io::message::critical);
      
      return false;
    }
  }

  if (i != num){
    io::messages.add("configuration file does not match topology: "
		     "not enough coordinates in VELOCITY block",
		     "In_Configuration",
		     io::message::error);
    return false;
  }
  
  return true;
  
}

bool io::In_Configuration::_read_box(math::Box &box, std::vector<std::string> &buffer,
				     math::boundary_enum const boundary)
{
  DEBUG(8, "read triclinic box");

  // no title in buffer!
  std::vector<std::string>::const_iterator it = buffer.begin(),
    to = buffer.end()-1;
  
  int bound;
  _lineStream.clear();
  _lineStream.str(*it);
  _lineStream >> bound;

  ++it;
  
  int i;
  
  for(i=0; it != to; ++i, ++it){

    if (i>=3){
      io::messages.add("bad line in TRICLINICBOX block","In_Configuration", io::message::error);
      return false;
    }
    
    _lineStream.clear();
    _lineStream.str(*it);

    _lineStream >> box(0)(i) >> box(1)(i) >> box(2)(i);
    
    if(_lineStream.fail())
      return false;
    if (!_lineStream.eof()) {
      std::string msg = "Warning, end of line not reached, but should have been: \n" + *it +  "\n";
      DEBUG(10, msg);
    }
  }
  
  // and check the boundary condition...
  if (math::boundary_enum(bound) != boundary){
    io::messages.add("Boundary condition from input file and from TRICLINICBOX do not match!"
		     " - using input file",
		     "In_Configuration", io::message::warning);
  }

  return true;
  
}

bool io::In_Configuration::_read_g96_box(math::Box &box, std::vector<std::string> &buffer)
{
  DEBUG(8, "read g96 box");

  // no title in buffer!
  std::vector<std::string>::const_iterator it = buffer.begin();
  
  _lineStream.clear();
  _lineStream.str(*it);
  
  _lineStream >> box(0)(0) >> box(1)(1) >> box(2)(2);
  box(0)(1) = box(0)(2) = box(1)(0) = box(1)(2) = box(2)(0) = box(2)(1) = 0.0;
    
  if(_lineStream.fail()){
    io::messages.add("failed to read BOX block",
		     "In_Configuration",
		     io::message::error);
    return false;
  }
  
  if (!_lineStream.eof()) {
    std::string msg = "Warning, end of line not reached, but should have been: \n" + *it +  "\n";
    DEBUG(10, msg);
  }

  return true;
  
}

bool io::In_Configuration::_read_flexv
(
 std::vector<double> &flexv,
 std::vector<std::string> &buffer,
 std::vector<topology::two_body_term_struct> const & constr,
 std::vector<topology::perturbed_two_body_term_struct> const & pert_constr)
{
  DEBUG(8, "read flexv");

  // no title in buffer!
  std::vector<std::string>::const_iterator it = buffer.begin(),
    to = buffer.end()-1;
  
  std::vector<topology::two_body_term_struct>::const_iterator 
    constr_it = constr.begin(),
    constr_to = constr.end();
  std::vector<topology::perturbed_two_body_term_struct>::const_iterator 
    pert_constr_it = pert_constr.begin(),
    pert_constr_to = pert_constr.end();

  std::vector<double>::iterator flexv_it = flexv.begin();

  int i, j, c, pc;
  double v, l;
  
  for(c=0; (it != to) && (constr_it != constr_to); ++it, ++constr_it, ++flexv_it, ++c){

    _lineStream.clear();
    _lineStream.str(*it);

    _lineStream >> i >> j >> l >> v;
    
    --i;
    --j;
    
    if(_lineStream.fail()){
      io::messages.add("failed to read FLEXV block",
		       "In_Configuration",
		       io::message::error);
      return false;
    }
    
    if (!_lineStream.eof()) {
      std::string msg = "Warning, end of line not reached, but should have been: \n" + *it +  "\n";
      DEBUG(10, msg);
    }

    if (i != int(constr_it->i) || j != int(constr_it->j))
      io::messages.add("wrong order in FLEXV block, constraints do not match",
		       "In_Configuration",
		       io::message::error);

    *flexv_it = v;

  }

  // and the perturbed constraints
  for(pc=0; (it != to) && (pert_constr_it != pert_constr_to);
       ++it, ++pert_constr_it, ++flexv_it, ++pc){

    _lineStream.clear();
    _lineStream.str(*it);

    _lineStream >> i >> j >> v;
    
    --i;
    --j;
    
    if(_lineStream.fail()){
      io::messages.add("Failed to read (perturbed) FLEXV block",
		       "In_Configuration",
		       io::message::error);
      return false;
    }

    if (!_lineStream.eof()) {
      std::string msg = "Warning, end of line not reached, but should have been: \n" + *it +  "\n";
      DEBUG(10, msg);
    }

    if (i != int(pert_constr_it->i) || j != int(pert_constr_it->j))
      io::messages.add("wrong order in FLEXV block, perturbed constraints do not match",
		       "In_Configuration",
		       io::message::error);

    *flexv_it = v;

  }
  
  if (c && !quiet)
    std::cout << "\tvelocities for " << c << " flexible constraints read in\n";
  if (pc && !quiet)
    std::cout << "\tvelocities for " << pc << " perturbed flexible constraints read in\n";
  
  return true;
  
}


bool io::In_Configuration::
_read_jvalue_av(std::vector<std::string> &buffer,
		std::vector<double> & jval_av,
		std::vector<topology::jvalue_restraint_struct> const & jval_res)
{
  DEBUG(8, "read jvalue averages");

  // no title in buffer!
  std::vector<std::string>::const_iterator it = buffer.begin(),
    to = buffer.end()-1;
  
  std::vector<topology::jvalue_restraint_struct>::const_iterator 
    jval_it = jval_res.begin(),
    jval_to = jval_res.end();

  jval_av.clear();
  
  int i, j, k, l;
  double av;
  
  if (buffer.size() - 1 != jval_res.size()){
    io::messages.add("number of J-restraints does not match with number of "
		     "continuation data", "in_configuration",
		     io::message::error);
    return false;
  }
  
  for( ; (it != to) && (jval_it != jval_to); ++it, ++jval_it){

    _lineStream.clear();
    _lineStream.str(*it);

    _lineStream >> i >> j >> k >> l >> av;
    
    if (int(jval_it->i) != i-1 ||
	int(jval_it->j) != j-1 ||
	int(jval_it->k) != k-1 ||
	int(jval_it->l) != l-1){

      io::messages.add("Wrong J-Value in JVALAVERAGE block",
		       "In_Configuration",
		       io::message::error);
      return false;

    }
    
    jval_av.push_back(av);
  }
  
  if (jval_it != jval_to || it != to){
    io::messages.add("Wrong number of J-Values in JVALAVERAGE block",
		     "In_Configuration",
		     io::message::error);
    return false;
  }
  
  return true;
}


bool io::In_Configuration::
_read_pscale_jrest(std::vector<std::string> &buffer,
		   configuration::Configuration::special_struct::pscale_struct &pscale,
		   std::vector<topology::jvalue_restraint_struct> const & jval_res)
{
  DEBUG(8, "read pscale jrest data");

  // no title in buffer!
  std::vector<std::string>::const_iterator it = buffer.begin(),
    to = buffer.end()-1;
  
  std::vector<topology::jvalue_restraint_struct>::const_iterator 
    jval_it = jval_res.begin(),
    jval_to = jval_res.end();

  if (buffer.size() - 1 != jval_res.size()){
    io::messages.add("number of J-restraints does not match with number of "
		     "J-restraints for periodic scaling", "in_configuration",
		     io::message::error);
    return false;
  }
  
  pscale.t.clear();
  pscale.scaling.clear();
  
  int i, j, k, l, s;
  double t;
  
  for( ; (it != to) && (jval_it != jval_to); ++it, ++jval_it){

    _lineStream.clear();
    _lineStream.str(*it);

    _lineStream >> i >> j >> k >> l >> s >> t;
    
    if (int(jval_it->i) != i-1 ||
	int(jval_it->j) != j-1 ||
	int(jval_it->k) != k-1 ||
	int(jval_it->l) != l-1){

      io::messages.add("Wrong J-Restraint in PSCALEJREST block",
		       "In_Configuration",
		       io::message::error);
      DEBUG(8, "wrong J-Restraint in PSCALEJREST block!");
      return false;
    }
    
    pscale.t.push_back(t);
    pscale.scaling.push_back(s);

    DEBUG(10, "\tt = " << t << "\tscaling = " << s);
  }
  
  if (jval_it != jval_to || it != to){
    io::messages.add("Wrong number of J-Restraints in PSCALEJREST block",
		     "In_Configuration",
		     io::message::error);
    return false;
  }
  
  return true;
}

bool io::In_Configuration::
_read_time(std::vector<std::string> &buffer,
	   double & t)
{
  DEBUG(8, "read time");

  // no title in buffer!
  std::vector<std::string>::const_iterator it = buffer.begin(),
    to = buffer.end()-1;
  
  _lineStream.clear();
  _lineStream.str(*it);

  int i;
  _lineStream >> i >> t;
  
  if (_lineStream.fail() || t < 0){
    io::messages.add("Could not read time from configuration file",
		     "In_Configuration", io::message::error);
    return false;
  }

  return true;
}

bool io::In_Configuration::
_read_time_step(std::vector<std::string> &buffer,
		simulation::Simulation & sim)
{
  DEBUG(8, "read time step");

  // no title in buffer!
  std::vector<std::string>::const_iterator it = buffer.begin(),
    to = buffer.end()-1;
  
  _lineStream.clear();
  _lineStream.str(*it);

  int i;
  double t;
  
  _lineStream >> i >> t;
  
  if (_lineStream.fail() || t < 0 || i < 0){
    io::messages.add("Could not read time from configuration file",
		     "In_Configuration", io::message::error);
    return false;
  }

  sim.steps() = i;
  sim.time() = t;
  
  return true;
}

bool io::In_Configuration::read_replica_information
(
 std::vector<util::Replica_Data> & replica_data,
 std::ostream & os
 )
{
  DEBUG(8, "read replica information");
  std::vector<std::string> buffer;
  buffer = m_block["REPDATA"];
  if (!m_block.size()){
    io::messages.add("no REPDATA block in first REPLICAFRAME!",
		     "In_Configuration",
		     io::message::error);
    return false;
  }
  
  block_read.insert("REPDATA");
  
  if (!quiet)
    os << "\treading REPDATA...\n";

  // no title in buffer!
  std::vector<std::string>::const_iterator it = buffer.begin(),
    to = buffer.end()-1;
  
  for(int i=0; it!=to; ++it, ++i){
    _lineStream.clear();
    _lineStream.str(*it);
  
    assert(unsigned(i) < replica_data.size());

    int st;
    
    _lineStream >> replica_data[i].ID 
		>> replica_data[i].run
		>> replica_data[i].Ti
		>> replica_data[i].li
		>> replica_data[i].epot_i
		>> replica_data[i].Tj
		>> replica_data[i].lj
		>> replica_data[i].epot_j
		>> replica_data[i].probability
		>> replica_data[i].switched
		>> st;

    replica_data[i].state = (util::state_enum)st;

    if (_lineStream.fail()){
      io::messages.add("Could not read replica information (REPDATA)",
		       "In_Configuration", io::message::error);
      return false;
    }
  }

  return true;
}

bool io::In_Configuration::check_coordinates
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation & sim,
 int num_coords,
 std::ostream & os)
{
  // do a quick estimate of the number of solvents
  if (topo.num_solvents() == 1){
    const unsigned int coords = num_coords - topo.num_solute_atoms();
    const unsigned int mols = coords / topo.solvent(0).num_atoms();
    
    if (mols * topo.solvent(0).num_atoms() != coords){
      io::messages.add("wrong number of coordinates!",
		       "in_configuration",
		       io::message::warning);
      os << "\twrong number of coordinates: " << mols * topo.solvent(0).num_atoms()
	 << " != " << coords << "\n";
    }
    if (topo.num_solvent_atoms() != coords){
      std::ostringstream os;
      os << "[Frame " << sim.steps() << "] resolvating: " 
	 << topo.num_solvent_atoms() / topo.solvent(0).num_atoms()
	 << " -> " << mols << " solvents";
      
      io::messages.add(os.str(),
		       "in_configuration",
		       io::message::notice);
      
      os << "\tresolvating! expected solvent atoms = " << topo.num_solvent_atoms()
	 << " got " << coords << std::endl;
      os << "\t(total coords: " << num_coords
	 << "\n\t solute coords: " << topo.num_solute_atoms()
	 << "\n\t solvent mols: " << mols
	 << "\n\t)" << std::endl;
      
      topo.resolvate(0, mols);
      conf.resize(num_coords);
      
      if (sim.multibath().bath_index()[sim.multibath().bath_index().size()-1].last_atom
	  != topo.num_atoms() - 1){
	sim.multibath().bath_index()[sim.multibath().bath_index().size()-1].last_atom
	  = topo.num_atoms() - 1;
	os << "\tadjusting temperature bath index\n";
      }
    }
  }
  
  return true;
}

