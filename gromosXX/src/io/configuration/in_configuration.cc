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
#include <math/transformation.h>
#include <util/umbrella_weight.h>

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
  DEBUG(7, "read configuration");
  
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
  read_cos_position(topo, conf, sim, os);
  read_velocity(topo, conf, sim, os);
  read_lattice_shifts(topo, conf, sim, os);
  read_box(topo, conf, sim, os);
  read_jvalue(topo, conf, sim, os);
  read_xray(topo, conf, sim, os);
  read_pscale(topo, conf, sim, os);
  read_flexv(topo, conf, sim, os);
  read_stochastic_integral(topo, conf, sim, os);
  read_perturbation(topo, sim, os);
  read_distance_restraint_averages(topo, conf, sim, os);
  read_nose_hoover_chains(topo, conf, sim, os);
  read_rottrans(topo, conf, sim, os);
  read_position_restraints(topo, conf, sim, os);
  read_leusbias(topo, conf, sim, os);
  
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
  
  DEBUG(8, "configuration read");
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
  
  DEBUG(8, "read replica");
  simulation::Parameter const & param = sim.param();

  const int switch_T = sim.param().replica.num_T;
  const int switch_l = sim.param().replica.num_l;
	
  DEBUG(9, "switch_T = " << switch_T << " switch_l = " << switch_l);

  const int rep_num = switch_T * switch_l;
  DEBUG(9, "=> " << rep_num << " replicas");

  if (rep_num < 1){
    io::messages.add("replica exchange with < 1 replica!",
		     "in_configuration",
		     io::message::error);
    return;
  }

  DEBUG(9, "resizing conf: " << rep_num);
  conf.resize(rep_num);
  DEBUG(9, "resizing replica data: " << rep_num);
  replica_data.resize(rep_num);

  block_read.clear();
  
  DEBUG(8, "reading in a frame");
  read_frame();

  // is it a multi-configuration (replica) file?
  std::vector<std::string> buffer = m_block["REPLICAFRAME"];
  
  if(!buffer.size()){
	DEBUG(8, "no REPLICAFRAME block. reading all into first configuration");
    assert(conf.size() > 0);
    conf[0].resize(topo.num_atoms());

    read_time(topo, conf[0], sim, os);
    read_position(topo, conf[0], sim, os);
    read_cos_position(topo, conf[0], sim, os);
    read_velocity(topo, conf[0], sim, os);
    read_lattice_shifts(topo, conf[0], sim, os);
    read_box(topo, conf[0], sim, os);
    read_jvalue(topo, conf[0], sim, os);
    read_pscale(topo, conf[0], sim, os);
    read_flexv(topo, conf[0], sim, os);
    read_stochastic_integral(topo, conf[0], sim, os);
    read_perturbation(topo, sim, os);
    read_distance_restraint_averages(topo, conf[0], sim, os);
    read_nose_hoover_chains(topo, conf[0], sim, os);
    read_rottrans(topo, conf[0], sim, os);
  
	DEBUG(10, "setting boundary type");
    conf[0].boundary_type = param.boundary.boundary;

	DEBUG(10, "copying configurations");
    for(unsigned int i=1; i<conf.size(); ++i)
      conf[i] = conf[0];

	DEBUG(10, "setting up replica information");
    // setup replica information
    {
      int i=0;
      for(int l=0; l < switch_l; ++l){
	for(int t=0; t < switch_T; ++t, ++i){
	  
	  assert(replica_data.size() > unsigned(i));
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
	  
	  replica_data[i].time = param.step.t0;
	}
      }
    }
  }
  else{
	DEBUG(8, "got REPLICAFRAME, it's a replica exchange configuration file!");
    read_replica_information(replica_data, os);

	DEBUG(8, "REPLICAFRAME block. reading " << conf.size() << " configurations");
    for(unsigned int i=0; i<conf.size(); ++i){
	  DEBUG(9, "topo contains " << topo.num_atoms() << " atoms");
      conf[i].resize(topo.num_atoms());
      
      read_time(topo, conf[i], sim, os);
      read_position(topo, conf[i], sim, os);
      read_cos_position(topo, conf[i], sim, os);
      read_velocity(topo, conf[i], sim, os);
      read_lattice_shifts(topo, conf[i], sim, os);
      read_box(topo, conf[i], sim, os);
      read_jvalue(topo, conf[i], sim, os);
      read_pscale(topo, conf[i], sim, os);
      read_flexv(topo, conf[i], sim, os);
      read_stochastic_integral(topo, conf[i], sim, os);
      read_perturbation(topo, sim, os);
      read_distance_restraint_averages(topo, conf[i], sim, os);
      read_nose_hoover_chains(topo, conf[i], sim, os);
      read_rottrans(topo, conf[i], sim, os);
      
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
	// read_stochastic_integral(topo, conf, sim, os);

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

bool io::In_Configuration::read_cos_position
(
 topology::Topology &topo, 
 configuration::Configuration &conf, 
 simulation::Simulation & sim,
 std::ostream & os)
{
    // read virtual sites for polarisation
  std::vector<std::string> buffer;

  if (sim.param().polarise.cos) {
    buffer = m_block["COSDISPLACEMENTS"];
    if (buffer.size()){
      
      check_coordinates(topo, conf, sim, buffer.size() - 1, os);
      
      if (!quiet)
        os << "\treading COSDISPLACEMENTS...\n";
      _read_cos_position(conf.current().posV, buffer, topo.num_atoms());
      
      conf.old().posV = conf.current().posV;
      
      block_read.insert("COSDISPLACEMENTS");
    }
    
    else{
      io::messages.add("no COSDISPLACEMENTS block found in input configuration."
              " Setting COS position to zero.",
              "in_configuration",
              io::message::notice);
      
      conf.current().posV = 0.0;
      conf.old().posV = 0.0;
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
    util::generate_velocities(sim.param(), sim.param().start.tempi, 
			      topo.mass(),
			      conf.current().vel,
			      conf.old().vel,
			      sim.param().start.ig,
			      os);
  }
  return true;
}

bool io::In_Configuration::read_lattice_shifts
(
 topology::Topology &topo, 
 configuration::Configuration &conf, 
 simulation::Simulation & sim,
 std::ostream & os)
{
  // read lattice_shifts
  std::vector<std::string> buffer;
  buffer = m_block["LATTICESHIFTS"];
  if (sim.param().start.read_lattice_shifts) {
    if (buffer.size()) {
      if (!quiet)
        os << "\treading LATTICESHIFTS...\n";
      _read_lattice_shifts(conf.special().lattice_shifts, buffer, topo.num_atoms());
      block_read.insert("LATTICESHIFTS");
    } else {
      io::messages.add("no LATTICESHIFTS block found in input configuration",
              "in_configuration",
              io::message::error);
      conf.special().lattice_shifts = 0.0;
      return false;
    }
  } else {
    if (buffer.size()) {
      io::messages.add("LATTICESHIFTS block provied but shifts reset to zero.",
              "in_configuration", io::message::warning);
      conf.special().lattice_shifts = 0.0;
    }
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
  conf.current().phi = conf.current().theta = conf.current().psi = 0.0;
  if(sim.param().boundary.boundary != math::vacuum){
    //read in GENBOX
    buffer = m_block["GENBOX"];
    if(buffer.size()){
       if (!quiet)
           os <<"\treading GENBOX...\n";
       _read_genbox(conf.current().box,conf.current().phi, conf.current().theta,
              conf.current().psi, buffer, sim.param().boundary.boundary);
       conf.old().box = conf.current().box;
       conf.old().phi = conf.current().phi;
       conf.old().theta = conf.current().theta;
       conf.old().psi = conf.current().psi;
       block_read.insert("GENBOX");
    } else {
    buffer = m_block["TRICLINICBOX"];
    if (buffer.size()){
      if (!quiet)
	os << "\treading TRICLINICBOX...\n";
      _read_box(conf.current().box,conf.current().phi, conf.current().theta,
              conf.current().psi,buffer, sim.param().boundary.boundary);
      conf.old().box = conf.current().box;
      conf.old().phi = conf.current().phi;
      conf.old().theta = conf.current().theta;
      conf.old().psi = conf.current().psi;
      block_read.insert("TRICLINICBOX");
      io::messages.add("TRICLINICBOX given"
		     " - output will be GENBOX",
		     "In_Configuration", io::message::notice);
    }
    else{
      buffer = m_block["BOX"];
      if (buffer.size() && (sim.param().boundary.boundary == math::rectangular ||
			    sim.param().boundary.boundary == math::truncoct)){
	if (!quiet)
	  os << "\treading BOX...\n";
          _read_g96_box(conf.current().box, buffer);
          conf.old().box = conf.current().box;
          conf.old().phi = conf.current().phi;
          conf.old().theta = conf.current().theta;
          conf.old().psi = conf.current().psi;
          io::messages.add("BOX given"
		     " - output will be GENBOX",
		     "In_Configuration", io::message::notice);
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
  }
  /* rotate the
   *postion
   *velocities
   *cos position
   *into the frame of the box
   */
 // math::Matrixl Rmat(math::transpose(math::rmat(conf.current().phi,
 //         conf.current().theta, conf.current().psi)));
  math::Matrixl Rmat((math::rmat(conf.current().phi,
          conf.current().theta, conf.current().psi)));
  DEBUG(10, "box \n " << math::m2s(math::Matrix(conf.current().box)));
          
  DEBUG(10, "Transformation Matrix \n" << math::m2s(Rmat))
  for (int i = 0, to = topo.num_atoms(); i < to; ++i) {
    DEBUG(10, "Position Cartesian: " << math::v2s(conf.current().pos(i)));
    conf.current().pos(i) = math::Vec(math::product(math::transpose(Rmat), conf.current().pos(i)));
    DEBUG(10, "Position Rotated  : " << math::v2s(conf.current().pos(i)));
    conf.current().posV(i) = math::Vec(math::product(math::transpose(Rmat), conf.current().posV(i)));
    conf.current().vel(i) = math::Vec(math::product(math::transpose(Rmat), conf.current().vel(i)));
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

  if (sim.param().jvalue.mode != simulation::jvalue_restr_off){

    if (sim.param().jvalue.mode == simulation::jvalue_restr_inst ||
        sim.param().jvalue.mode == simulation::jvalue_restr_inst_weighted){
      if (sim.param().jvalue.read_av)
	io::messages.add("instantaneous J-value restraints, ignoring reading of averages",
			 "in_configuration",
			 io::message::warning);
    }
    else if (!sim.param().jvalue.read_av){

      buffer = m_block["JVALUERESEXPAVE"];
      if (buffer.size()){
	block_read.insert("JVALUERESEXPAVE");

	io::messages.add("re-initialising J-value averages, non-continuous simulation",
			 "in_configuration",
			 io::message::warning);
      }
      else{
	io::messages.add("initialising J-value averages",
			 "in_configuration",
			 io::message::notice);
      }
    }
    else {

      buffer = m_block["JVALUERESEXPAVE"];
      if (buffer.size())
      {
	block_read.insert("JVALUERESEXPAVE");
	_read_jvalue_av(buffer, conf.special().jvalue_av, topo.jvalue_restraints());
      }
      else{
	io::messages.add("reading in of J-value averages requested "
			 "but JVALUERESEXPAVE block not found",
			 "in_configuration",
			 io::message::error);
	return false;
      }
    }
  } // jvalue averages
  
  if (sim.param().jvalue.le) {
    buffer = m_block["JVALUERESEPS"];
    if (!sim.param().jvalue.read_av){
        if (buffer.size()){
	block_read.insert("JVALUERESEPS");
	io::messages.add("re-initialising J-value local elevation epsilons, non-continuous simulation",
			 "in_configuration",
			 io::message::warning);
      } else{
	io::messages.add("initializing J-value local elevation epsilons",
			 "in_configuration",
			 io::message::notice);
      }
    } else {
      if (buffer.size()) {
	block_read.insert("JVALUERESEPS");
	_read_jvalue_le(buffer, conf.special().jvalue_epsilon, 
                        topo.jvalue_restraints(), 
                        sim.param().jvalue.ngrid);
      } else{
	io::messages.add("reading in of J-value local elevation epsilons "
                         "requested but JVALUERESEPS block not found",
                         "in_configuration",
                         io::message::error);
	return false;
      }
    }
  } // jvalue local elevation
  return true;
}
bool io::In_Configuration::read_xray
(
        topology::Topology &topo,
        configuration::Configuration &conf,
        simulation::Simulation & sim,
        std::ostream & os) {
  std::vector<std::string> buffer;

  if (sim.param().xrayrest.xrayrest != simulation::xrayrest_off) {

    if (sim.param().xrayrest.xrayrest == simulation::xrayrest_inst) {
      if (!sim.param().xrayrest.readavg)
        io::messages.add("instantaneous Xray restraints, ignoring reading of averages",
              "in_configuration",
              io::message::warning);
    } else if (!sim.param().xrayrest.readavg) {

      io::messages.add("initialising Xray-restraint averages",
              "in_configuration",
              io::message::notice);

    } else {

      buffer = m_block["XRAYRESEXPAVE"];
      if (buffer.size()) {
        block_read.insert("XRAYRESEXPAVE");
        io::messages.add("reading Xray-Restraint-Averages from last configuration file",
                "in_configuration",
                io::message::warning);
        _read_xray_av(buffer, conf.special().xray_rest, topo.xray_restraints());
      } else {
        io::messages.add("reading in of Xray-restraints averages requested "
                "but XRAYRESEXPAVE block not found",
                "in_configuration",
                io::message::error);
        return false;
      }
    }
  } // xray averages

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
      buffer = m_block["JVALUEPERSCALE"];
      if (buffer.size()){
	block_read.insert("JVALUEPERSCALE");
	io::messages.add("re-initialising J-value restraints periodic scaling data, non-continuous simulation",
			 "in_configuration",
			 io::message::warning);
      }
      else{
	io::messages.add("initialising J-value restraints periodic scaling data",
			 "in_configuration",
			 io::message::notice);
      }
    }
    else {
      buffer = m_block["JVALUEPERSCALE"];
      if (buffer.size())
      {
	block_read.insert("JVALUEPERSCALE");
	_read_pscale_jrest(buffer, conf.special().pscale, topo.jvalue_restraints());
      }
      else{
	io::messages.add("reading in of J-value restraints periodic scaling data requested "
			 "but JVALUEPERSCALE block not found",
			 "in_configuration",
			 io::message::error);
	return false;
      }
    }
  } // JVALUEPERSCALE
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

bool io::In_Configuration::read_stochastic_integral
(
 topology::Topology &topo, 
 configuration::Configuration &conf, 
 simulation::Simulation & sim,
 std::ostream & os)
{
  std::vector<std::string> buffer;
  if (sim.param().stochastic.sd){

    conf.old().stochastic_integral.resize(topo.num_atoms());
    conf.current().stochastic_integral.resize(topo.num_atoms());
    
    buffer = m_block["STOCHINT"];
    if (sim.param().stochastic.generate_integral == false) {
      if (buffer.size()){
        block_read.insert("STOCHINT");
        if (!quiet)
          os << "\treading STOCHINT...\n";
        
        _read_stochastic_integral(conf.current().stochastic_integral, buffer,
                topo.num_atoms(), conf.current().stochastic_seed);
      }
      else {
        io::messages.add("could not read stochastic integrals from configuration",
                "In_Configuration", io::message::error);
      }
    }
  }
  return true;
}

bool io::In_Configuration::read_perturbation
(
 topology::Topology &topo, 
 simulation::Simulation & sim,
 std::ostream & os)
{
  std::vector<std::string> buffer;
  if (sim.param().perturbation.perturbation){    
    buffer = m_block["PERTDATA"];
    if (sim.param().perturbation.read_initial) {
      if (buffer.size()){
        block_read.insert("PERTDATA");
        if (!quiet)
          os << "\treading PERTDATA...\n";
        
        _read_pertdata(topo, buffer);
      } else {
        io::messages.add("could not read perturbation data (PERTDATA) from configuration",
                "In_Configuration", io::message::error);
      }
    }
  }
  return true;
}

bool io::In_Configuration::read_distance_restraint_averages
(
 topology::Topology &topo, 
 configuration::Configuration &conf, 
 simulation::Simulation & sim,
 std::ostream & os)
{
  std::vector<std::string> buffer;
  if (sim.param().distanceres.distanceres < 0){
    
    buffer = m_block["DISRESEXPAVE"];
    if (buffer.size()){
      block_read.insert("DISRESEXPAVE");
      if (!quiet)
	os << "\treading DISRESEXPAVE...\n";

      if (sim.param().distanceres.read)
        _read_distance_restraint_averages(buffer, topo.distance_restraints(),
                                          conf.special().distanceres_av);
      else
        io::messages.add("distance restraint averages found but not read.",
                         "in_configuration",
                         io::message::warning);
    }
    else{
      if (sim.param().distanceres.read)
        io::messages.add("no DISRESEXPAVE block in configuration.",
                         "in_configuration",
                         io::message::error);
    }
  }
  return true;
}

bool io::In_Configuration::read_nose_hoover_chains
(
 topology::Topology &topo, 
 configuration::Configuration &conf, 
 simulation::Simulation & sim,
 std::ostream & os)
{
  std::vector<std::string> buffer;
  if (sim.param().multibath.nosehoover > 1){
    
    buffer = m_block["NHCVARIABLES"];
    if (buffer.size()){
      block_read.insert("NHCVARIABLES");
      if (!quiet)
	os << "\treading NHCVARIABLES...\n";

      if (sim.param().start.read_nosehoover_chains)
        _read_nose_hoover_chain_variables(buffer, sim.multibath());
      else
        io::messages.add("Nose-Hoover-Chains variables found but not read.",
                         "in_configuration",
                         io::message::warning);
    }
    else{
      if (sim.param().start.read_nosehoover_chains)
        io::messages.add("no NHCVARIABLES block in configuration.",
                         "in_configuration",
                         io::message::error);
    }
  }
  return true;
}

bool io::In_Configuration::read_rottrans
(
 topology::Topology &topo, 
 configuration::Configuration &conf, 
 simulation::Simulation & sim,
 std::ostream & os)
{
  std::vector<std::string> buffer;
  if (sim.param().rottrans.rottrans){
    
    buffer = m_block["ROTTRANSREFPOS"];
    if (buffer.size()){
      block_read.insert("ROTTRANSREFPOS");
      if (!quiet)
	os << "\treading ROTTRANSREFPOS...\n";

      if (sim.param().start.read_rottrans)
        _read_rottrans(buffer, sim.param().rottrans.last, conf.special().rottrans_constr);
      else
        io::messages.add("Initial settings for roto-translational constraints found but not read.",
                         "in_configuration",
                         io::message::warning);
    }
    else{
      if (sim.param().start.read_rottrans)
        io::messages.add("no ROTTRANSREFPOS block in configuration.",
                         "in_configuration",
                         io::message::error);
    }
  }
  return true;
}

bool io::In_Configuration::read_position_restraints
(
 topology::Topology &topo,
 configuration::Configuration &conf, 
 simulation::Simulation & sim,
 std::ostream & os)
{
  if (sim.param().posrest.posrest == simulation::posrest_off)
    return true;
  
  bool result = true;
  
  std::vector<std::string> buffer;
  if (!sim.param().posrest.read) { // read means from spec file!
    buffer = m_block["REFPOSITION"];
    if (buffer.size()) {
      block_read.insert("REFPOSITION");
      if (!quiet)
	os << "\treading REFPOSITION...\n";

      conf.special().reference_positions.resize(topo.num_atoms());
      result = result && 
              _read_position(conf.special().reference_positions, buffer,
              topo.num_atoms(), "REFPOSITION");
    } else {
      io::messages.add("no REFPOSITION block in configuration.",
                       "in_configuration", io::message::error);
      return false;
    }
    
    if (sim.param().posrest.posrest == simulation::posrest_bfactor) {
      buffer = m_block["BFACTOR"];
      if (buffer.size()) {
        block_read.insert("BFACTOR");
        if (!quiet)
          os << "\treading BFACTOR...\n";

        conf.special().bfactors.resize(topo.num_atoms());
        result = result && 
                _read_bfactor(conf.special().bfactors, buffer, false);
      } else {
        io::messages.add("no BFACTOR block in configuration.",
                "in_configuration", io::message::error);
        return false;
      }
    }
  }
  
  return result;
}

bool io::In_Configuration::read_leusbias
(
 topology::Topology &topo,
 configuration::Configuration &conf,
 simulation::Simulation & sim,
 std::ostream & os)
{
  if (sim.param().localelev.localelev == simulation::localelev_off)
    return true;

  std::vector<std::string> buffer;
  if (!sim.param().localelev.read) { // read means from spec file!
    buffer = m_block["LEUSBIAS"];
    if (buffer.size()) {
      block_read.insert("LEUSBIAS");
      if (!quiet)
	os << "\treading LEUSBIAS...\n";

      bool result = _read_leusbias(conf.special().umbrellas, buffer, false);
      if (!quiet) {
        os << "\t\t" << conf.special().umbrellas.size() <<
                " umbrellas found.\n";
      }
      return result;
    } else {
      io::messages.add("no LEUSBIAS block in configuration.",
                       "in_configuration", io::message::error);
      return false;
    }
  } else {
    buffer = m_block["LEUSBIAS"];
    if (buffer.size()) {
      block_read.insert("LEUSBIAS");
      io::messages.add("LEUSBIAS block in configuration is ignored",
                       "in_configuration", io::message::warning);
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

bool io::In_Configuration::_read_cos_position(math::VArray &pos, 
					     std::vector<std::string> &buffer,
					     int const num)
{
  DEBUG(8, "read COSDISPLACEMENTS");
  
  // no title in buffer!
  std::vector<std::string>::const_iterator it = buffer.begin(),
    to = buffer.end()-1;
  
  int i;

  if (pos.size() < unsigned(num)){
    io::messages.add("configuration: too many cos coordinates for given topology",
		     "in_configuration",
		     io::message::critical);
    std::cout << "cos position size is : " << pos.size() << " and num is : " << 
              num << std::endl;
    return false;
  }

  for(i=0; it != to; ++i, ++it){
    if (i >= num){
      io::messages.add("configuration file does not match topology: "
		       "too many coordinates in COSDISPLACEMENTS block",
		       "In_Configuration",
		       io::message::error);
      break;
    }
    
    _lineStream.clear();
    _lineStream.str(*it);
    _lineStream >> pos(i)(0) >> pos(i)(1) >> pos(i)(2);
    
    if(_lineStream.fail()){
      io::messages.add("bad line in COSDISPLACEMENTS block",
		       "In_Configuration",
		       io::message::error);      
      return false;
    }
  }

  if (i != num){
    io::messages.add("configuration file does not match topology: "
		     "not enough coordinates in COSDISPLACEMENTS block",
		     "In_Configuration",
		     io::message::error);
    return false;
  }

  return true;
}

bool io::In_Configuration::_read_position(math::VArray &pos, std::vector<std::string> &buffer,
					  int const num, std::string blockname)
{
  DEBUG(8, "read position");

  // no title in buffer!
  std::vector<std::string>::const_iterator it = buffer.begin(),
    to = buffer.end()-1;
  
  std::string s1, s2;
  int i, n, nr;

  std::istringstream _lineStream;

  for(i=0; it != to; ++i, ++it){
    if (i >= num){
      io::messages.add("configuration file does not match topology: "
		       "too many coordinates in "+blockname+" block",
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
      io::messages.add("bad line in "+blockname+" block",
		       "In_Configuration",
		       io::message::critical);
      return false;
    }
  }

  if (i != num){
    io::messages.add("configuration file does not match topology: "
		     "not enough coordinates in "+blockname+" block",
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

bool io::In_Configuration::_read_lattice_shifts(math::VArray &shift, 
					  std::vector<std::string> &buffer,
					  int const num)
{
  DEBUG(8, "read lattice shifts");

  // no title in buffer!
  std::vector<std::string>::const_iterator it = buffer.begin(),
    to = buffer.end()-1;
  
  std::string s1, s2;
  int i;

  for(i=0; it != to; ++i, ++it){
    if (i >= num){
      io::messages.add("configuration file does not match topology: "
		       "too many coordinates in LATTICESHIFTS block",
		       "In_Configuration", io::message::error);
      break;
    }
   
    _lineStream.clear();
    _lineStream.str(*it);
    _lineStream >> shift(i)(0) >> shift(i)(1) >> shift(i)(2);
    
    if(_lineStream.fail()){
      io::messages.add("bad line in LATTICESHIFTS block",
		       "In_Configuration",
		       io::message::critical);
      return false;
    }
  }

  if (i != num){
    io::messages.add("configuration file does not match topology: "
		     "not enough coordinates in LATTICESHIFTS block",
		     "In_Configuration",
		     io::message::error);
    return false;
  }
  return true;
}

bool io::In_Configuration::_read_genbox(math::Box &box, double &phi, 
        double &theta, double &psi,
        std::vector<std::string> &buffer,
				     math::boundary_enum const boundary)
{
  DEBUG(8, "read genbox");

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

    if (i>=4){
      io::messages.add("bad line in GENBOX block","In_Configuration", io::message::error);
      return false;
    }
    
    _lineStream.clear();
    _lineStream.str(*it);
    //point of origin is ignored so far
    if(i<3)
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
    io::messages.add("Boundary condition from input file and from GENBOX do not match!"
		     " - using input file",
		     "In_Configuration", io::message::warning);
  }
  //change from genbox to triclinicbox...
  long double a, b, c, alpha, beta, gamma;// phi, theta, psi;
  a = box(0)(0);
  b = box(1)(0);
  c = box(2)(0);
  alpha = math::Pi*box(0)(1)/180;
  beta  = math::Pi*box(1)(1)/180;
  gamma = math::Pi*box(2)(1)/180;
  if ( boundary == math::rectangular  &&(box(0)(1)!=90.0 || box(1)(1)!=90.0 ||box(2)(1)!=90.0) ){
    io::messages.add("Rectangular box, but box angles != 90",
		     "In_Configuration", io::message::error);
  }
  long double cosdelta=(cosl(alpha)-cosl(beta)*cosl(gamma))/(sinl(beta)*sinl(gamma));
  long double sindelta=sqrtl(1-cosdelta*cosdelta);  
   
  math::Vecl SBx(a, 0.0, 0.0);
  math::Vecl SBy(b*cosl(gamma), 
          b*sinl(gamma), 
          0.0);
  math::Vecl SBz(c*cosl(beta), 
          c*cosdelta*sinl(beta), 
          c*sindelta*sinl(beta));
  
  phi   = math::Pi*box(0)(2)/180;
  theta = math::Pi*box(1)(2)/180;
  psi   = math::Pi*box(2)(2)/180;

  box(0)=math::Vec(SBx);
  box(1)=math::Vec(SBy);
  box(2)=math::Vec(SBz);
  
  /* stay in the frame of the box -> don't rotate
  phi   = math::Pi*box(0)(2)/180;
  theta = math::Pi*box(1)(2)/180;
  psi   = math::Pi*box(2)(2)/180;
  
  const math::Vecl Rx(cosl(theta)*cosl(phi),
          cosl(theta)*sinl(phi),
          -sinl(theta));
  const math::Vecl Ry(sinl(psi)*sinl(theta)*cosl(phi)-cosl(psi)*sinl(phi),
          sinl(psi)*sinl(theta)*sinl(phi)+cosl(psi)*cosl(phi),
          sinl(psi)*cosl(theta));
  const math::Vecl Rz(cosl(psi)*sinl(theta)*cosl(phi)+sinl(psi)*sinl(phi),
          cosl(psi)*sinl(theta)*sinl(phi)+(-sinl(psi)*cosl(phi)),
          cosl(psi)*cosl(theta));
  
  math::Matrixl Rmat(Rx,Ry,Rz);
  
  // we have to convert Vecl to Vec by hand.
  box(0)=math::Vec(product(Rmat,SBx));
  box(1)=math::Vec(product(Rmat,SBy));
  box(2)=math::Vec(product(Rmat,SBz));
 
  */

  
  // and check the boundary condition...
  if (math::boundary_enum(bound) != boundary){
    io::messages.add("Boundary condition from input file and from TRICLINICBOX do not match!"
		     " - using input file",
		     "In_Configuration", io::message::warning);
  }
  for(int i=0; i<3; i++){
      for(int j=0; j<3; j++){
          if(fabs(box(i)(j))<math::epsilon)
              box(i)(j)=0.0;                    
      }
  }
  
  return true;
}


bool io::In_Configuration::_read_box(math::Box &box, double &phi, double &theta,
        double &psi,std::vector<std::string> &buffer,
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
  if ( boundary == math::rectangular  &&
          ( math::dot(box(0),box(1))!=0.0 
          || math::dot(box(0),box(2))!=0.0  
          || math::dot(box(1),box(2))!=0.0 ) ){
    io::messages.add("Rectangular box, but box angles != 90",
		     "In_Configuration", io::message::error);
  }
  
    //find phi, theta and psi
    math::Matrixl Rmat = (math::rmat(box));
    long double R11R21 = sqrtl(Rmat(0, 0) * Rmat(0, 0) + Rmat(0, 1) * Rmat(0, 1));
    if (R11R21 == 0.0) {
      theta = -math::sign(Rmat(0, 2)) * M_PI / 2;
      psi = 0.0;
      phi = -math::sign(Rmat(1, 0)) * acosl(math::costest(Rmat(1, 1)));
    } else {
      theta = -math::sign(Rmat(0, 2)) * acosl(math::costest(R11R21));
      long double costheta = cosl(theta);
      psi = math::sign(Rmat(1, 2) / costheta) * acosl(math::costest(Rmat(2, 2) / costheta));
      phi = math::sign(Rmat(0, 1) / costheta) * acosl(math::costest(Rmat(0, 0) / costheta));

  }

  //rotate box to frame of box
  DEBUG(10, "Rmat        : \n" << math::m2s(Rmat));
  DEBUG(10, "Rmat,transpo: \n" << math::m2s(math::transpose(Rmat)));
  DEBUG(10, "original box: \n" << math::m2s(math::Matrix(box)));
  //box=math::product(math::transpose(Rmat),box);
 // box=math::product((Rmat),box);
  
  math::Box m(0.0);
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      for (int k = 0; k < 3; ++k)
        m(j)(i) += Rmat(i, k) * box(j)(k);
  box = m;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
         if (fabs(box(i)(j)) <= math::epsilon)
              box(i)(j) = 0.0;
  DEBUG(10, "box in frame of box: \n" << math::m2s(math::Matrix(box)));
  DEBUG(10, "box(0): "<< math::v2s(box(0)));
  DEBUG(10, "box(0,1): " << box(0)(1));
  //box from scratch...
  math::Matrixl Smat = (math::smat(box, boundary));
  DEBUG(10, "Smat \n" << math::m2s(Smat));
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

bool io::In_Configuration::_read_stochastic_integral
(
 math::VArray & stochastic_integral,
 std::vector<std::string> &buffer,
 int num_atoms, std::string &seed)
{
  DEBUG(8, "read stochastic integral");

  // no title in buffer!
  std::vector<std::string>::const_iterator it = buffer.begin(),
    to = buffer.end()-2;
  
  int i, j, c;
  std::string a, r;

  for(c=0; it != to; ++it, ++c){

    if (c >= num_atoms){
      std::cout << "too much data in STOCHINT block! got: " << c
		<< " expected: " << num_atoms << std::endl;
      
      io::messages.add("too much data in STOCHINT block!",
		       "in_configuration",
		       io::message::error);
      return false;
    }
    
    _lineStream.clear();
    _lineStream.str(*it);

    _lineStream >> i >> a >> r >> j
		>> stochastic_integral(c)(0)
		>> stochastic_integral(c)(1)
		>> stochastic_integral(c)(2);

    if(_lineStream.fail()){
      std::cout << "failed to read stochastic integral (STOCHINT block) line: " << c
		<< "\n" << *it << std::endl;
      
      io::messages.add("failed to read stochastic integral (STOCHINT block)",
		       "In_Configuration",
		       io::message::error);
      return false;
    }
    
    if (!_lineStream.eof()) {
      std::string msg = "Warning, end of line not reached, but should have been: \n" + *it +  "\n";
      DEBUG(10, msg);
    }
  }

  if (c && !quiet)
    std::cout << "\tstochastic integrals for " << c << " atoms read\n";

  if (c != num_atoms){
    std::cout << "not enough stochastic integrals in STOCHINT block! : "
	      << " got: " << c << " expected: " << num_atoms << std::endl;
    
    io::messages.add("not enough stochastic integrals or seed missing in "
                     "STOCHINT block!", "In_Configuration", io::message::error);
    return false;
  }
  
  // get the seed and trimm of whitespace
  seed = *it;
  std::string::size_type pos = seed.find_last_not_of(' ');
  if(pos != std::string::npos) {
    seed.erase(pos + 1);
    pos = seed.find_first_not_of(' ');
    if(pos != std::string::npos) seed.erase(0, pos);
  }
  else seed.erase(seed.begin(), seed.end());
  
  DEBUG(12, "trimmed seed: '" << seed << "'");
  return true;
}

bool io::In_Configuration::_read_pertdata(topology::Topology & topo,
        std::vector<std::string> & buffer)
{
  std::string s;
  _lineStream.clear();
  _lineStream.str(concatenate(buffer.begin(), buffer.end(), s));
  
  double lambda;
  _lineStream >> lambda;
  if (_lineStream.fail()) {
    io::messages.add("Bad line in PERTDATA block.", "In_Configuration",
            io::message::error);
    return false;
  }
  
  // we have to set lambda twice to make sure old_lambda() is also
  // set to lambda.
  topo.lambda(lambda);
  topo.lambda(lambda);
  topo.update_for_lambda();
  return true;
}

bool io::In_Configuration::_read_distance_restraint_averages
(
 std::vector<std::string> &buffer,
 const std::vector<topology::distance_restraint_struct> &distanceress,
 std::vector<double> &distanceres_av
 )
{
  DEBUG(8, "read distance restaint averages");
  
  std::vector<topology::distance_restraint_struct>::const_iterator 
    distanceress_it = distanceress.begin(),
    distanceress_to = distanceress.end();
  
  std::string s;
  
  _lineStream.clear();
  _lineStream.str(concatenate(buffer.begin(), buffer.end(), s));
  
  for( ;distanceress_it != distanceress_to; ++distanceress_it){
    double ave;
    _lineStream >> ave;
    
    if (_lineStream.fail() || ave < 0.0) {
      io::messages.add("Could not read averages from DISRESEXPAVE block",
		       "In_Configuration",
		       io::message::error);
      return false;
    }
    
    distanceres_av.push_back(ave);  
  }
  
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
  
  double av;
  
  if (buffer.size() - 1 != jval_res.size()){
    std::cout << "JVALUERESEXPAVE: " << buffer.size() - 1
	      << " but restraints: " << jval_res.size()
	      << std::endl;

    io::messages.add("number of J-restraints does not match with number of "
		     "continuation data", "in_configuration",
		     io::message::error);
    return false;
  }
  
  for( ; (it != to) && (jval_it != jval_to); ++it, ++jval_it){

    _lineStream.clear();
    _lineStream.str(*it);

    _lineStream >> av;
    if (_lineStream.fail()){
      io::messages.add("Bad value in JVALUERESEXPAVE block",
		       "In_Configuration", io::message::error);
      return false;
    }
    jval_av.push_back(av);
  }
  
  if (jval_it != jval_to || it != to){
    io::messages.add("Wrong number of J-Values in JVALUERESEXPAVE block",
		     "In_Configuration",
		     io::message::error);
    return false;
  }
  
  return true;
}

bool io::In_Configuration::
_read_jvalue_le(std::vector<std::string> &buffer,
		std::vector<std::vector<double> > & jval_epsilon,
		std::vector<topology::jvalue_restraint_struct> const & jval_res,
                unsigned int const & grid_size)
{
  DEBUG(8, "read jvalue local elevation epsilon");

  // no title in buffer!
  std::vector<std::string>::const_iterator it = buffer.begin(),
    to = buffer.end()-1;
  
  std::vector<topology::jvalue_restraint_struct>::const_iterator 
    jval_it = jval_res.begin(),
    jval_to = jval_res.end();

  jval_epsilon.clear();
  
  if (buffer.size() - 1 != jval_res.size()){
    std::cout << "JVALRESEPSILON size: " << buffer.size() - 1
	      << " but restraints size: " << jval_res.size()
	      << std::endl;

    io::messages.add("number of J-restraints does not match with number of "
		     "LE continuation data", "in_configuration",
		     io::message::error);
    return false;
  }
  
  for( ; (it != to) && (jval_it != jval_to); ++it, ++jval_it){
    std::vector<double> eps(grid_size, 0.0);
    _lineStream.clear();
    _lineStream.str(*it);
    
    for(unsigned int i = 0; i < grid_size; ++i) 
      _lineStream >> eps[i];
    
    if (_lineStream.fail()){
      io::messages.add("Bad value in JVALRESEPSILON block",
		       "In_Configuration", io::message::error);
      return false;
    }
    jval_epsilon.push_back(eps);
  }
  
  if (jval_it != jval_to || it != to){
    io::messages.add("Wrong number of J-Values in JVALRESEPSILON block",
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
		     "periodic scaling data", "in_configuration",
		     io::message::error);
    return false;
  }
  
  pscale.t.clear();
  pscale.scaling.clear();

  for( ; (it != to) && (jval_it != jval_to); ++it, ++jval_it){
    _lineStream.clear();
    _lineStream.str(*it);
    int s;
    double t;
    _lineStream >> s >> t;

    if (_lineStream.fail()) {
      io::messages.add("Bad line in JVALUEPERSCALE block."
		     "periodic scaling", "in_configuration",
		     io::message::error);
    return false;
    }
    pscale.t.push_back(t);
    pscale.scaling.push_back(s);

    DEBUG(10, "\tt = " << t << "\tscaling = " << s);
  }
  
  if (jval_it != jval_to || it != to){
    io::messages.add("Wrong number of J-Restraints in JVALUEPERSCALE block",
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
		>> replica_data[i].time
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
  if (topo.num_solute_atoms() != 0) {
    const unsigned int coords = num_coords;
    if (topo.num_solute_atoms() > coords) {
      // Here, we have to abort: seg faults in resolvate.
      return false;
    }
  }
  
  // do a quick estimate of the number of solvents
  if (sim.param().system.nsm > 0 &&
      topo.num_solvents() == 1){
    const unsigned int coords = num_coords - topo.num_solute_atoms();
    
    if (topo.num_solvent_atoms() != coords) {
      // resolvating is very error prone. We disable it here
      return false;
    }
      /*
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
       */
  }
  return true;
}

bool io::In_Configuration::_read_bfactor(
math::SArray & b,
std::vector<std::string> &buffer, bool hasTitle)
{
  DEBUG(8, "read bfactor");
  // no title in buffer?
  std::vector<std::string>::const_iterator it = buffer.begin() + (hasTitle ? 1 : 0),
    to = buffer.end()-1;
  std::string s1, s2;
  int i, n, nr;
  int num = b.size();

  std::istringstream _lineStream;

  for(i=0; it != to; ++i, ++it){
    if (i >= num){
      io::messages.add("configuration file does not match topology: "
		       "too many B-factors in BFACTOR block",
		       "In_Configuration",
		       io::message::error);
      break;
    }

    _lineStream.clear();
    _lineStream.str(*it);
    // ignore first 4 fields
    _lineStream >> n >> s1 >> s2 >> nr;
    _lineStream >> b(i);

    if(_lineStream.fail()){
      io::messages.add("bad line in BFACTOR block",
		       "In_Configuration",
		       io::message::critical);
      return false;
    }
  }

  if (i != num){
    io::messages.add("configuration file does not match topology: "
		     "not enough B-factors in BFACTOR block",
		     "In_Configuration",
		     io::message::error);
    return false;
  }
  return true;
  
}

bool io::In_Configuration::_read_leusbias(
            std::vector<util::Umbrella> & umbrellas,
            std::vector<std::string> &buffer, bool reset) {
  DEBUG(8, "read LEUSBIAS");
  std::istringstream _lineStream;
  std::string s;
  _lineStream.clear();
  _lineStream.str(concatenate(buffer.begin(), buffer.end()-1, s));

  unsigned int num_umb;
  _lineStream >> num_umb;
  if (_lineStream.fail()) {
    io::messages.add("bad line in LEUSBIAS block: NUMUMB",
            "In_Configuration",
            io::message::error);
    return false;
  }
  for(unsigned int x = 0; x < num_umb; ++x) {
    int id, dim;
    _lineStream >> id >> dim;

    if(_lineStream.fail()){
      io::messages.add("bad line in LEUSBIAS block: NDIM",
		       "In_Configuration",
		       io::message::error);
      return false;
    }
    // create the umbrella
    util::Umbrella u(id, dim);

    _lineStream >> u.force_constant;
    if(_lineStream.fail()){
      io::messages.add("bad line in LEUSBIAS block: CLES",
		       "In_Configuration",
		       io::message::error);
      return false;
    }

    // loop over dimensions
    for (unsigned int i = 0; i < u.dim(); ++i) {
      int type, form;
      _lineStream >> type >> form >> u.width[i] >> u.cutoff[i]
              >> u.num_grid_points[i] >> u.grid_min[i] >> u.grid_max[i];
      if (_lineStream.fail()) {
        io::messages.add("LEUSBIAS block: Could not read umbrella definition",
                "In_Configuration",
                io::message::error);
        return false;
      }
      u.variable_type[i] = util::Umbrella::variable_type_enum(type);
      DEBUG(10, "variable type: " << u.variable_type[i]);
      u.functional_form[i] = util::Umbrella::functional_form_enum(form);
      DEBUG(10, "functional form: " << u.functional_form[i]);
    } // for dimensions (grid properties)

    if (reset) {
      umbrellas.push_back(u);
      continue;
    }

    // the rest of this block is stored and read later. We have to save the
    // stream position and the block.
    // the rest ist read in util::Umbrella::read_configuration. This
    // method must not be called before the weight factory is set!
    u.configuration_block_pos = _lineStream.tellg();
    u.configuration_block = _lineStream.str();
    umbrellas.push_back(u);
  } // for umbrellas
  return true;
}

bool io::In_Configuration::_read_nose_hoover_chain_variables(
std::vector<std::string> &buffer,
simulation::Multibath & multibath) {
  DEBUG(8, "read nose hoover chain variables");

  // no title in buffer!
  std::vector<std::string>::const_iterator it = buffer.begin(),
    to = buffer.end()-1;
  simulation::Multibath::iterator b_it = multibath.begin(),
          b_to = multibath.end();

  const unsigned int num_zeta = b_it->zeta.size();
  
  for (; it != to && b_it != b_to; ++it, ++b_it) {
    _lineStream.clear();
    _lineStream.str(*it);

    for(unsigned int i = 0; i < num_zeta; ++i)
      _lineStream >> b_it->zeta[i];

    if (_lineStream.fail()) {
      io::messages.add("Could not read Nose-Hoover-Chains from configuration file",
              "In_Configuration", io::message::error);
      return false;
    }
  }

  if (b_it != b_to) {
    io::messages.add("Could not read Nose-Hoover-Chains: Not enough lines (baths)",
            "In_Configuration", io::message::error);
    return false;
  }
  
  if (it != to) {
    io::messages.add("Could not read Nose-Hoover-Chains: Too many lines (baths)",
            "In_Configuration", io::message::error);
    return false;
  }

  return true;  
}

bool io::In_Configuration::_read_rottrans(
std::vector<std::string> &buffer, unsigned int last,
configuration::Configuration::special_struct::rottrans_constr_struct & rottrans) {
  DEBUG(8, "read configuration for roto-translational constraints");

  // no title in buffer!
  std::vector<std::string>::const_iterator it = buffer.begin(),
    to = buffer.end()-1;
  
  unsigned int i;
  for (i = 0; it != to && i < 3; ++it, ++i) {
    _lineStream.clear();
    _lineStream.str(*it);

    for(unsigned int j = 0; j < 3; ++j)
      _lineStream >> rottrans.theta_inv_trans(i,j);

    if (_lineStream.fail()) {
      io::messages.add("ROTTRANSREFPOS block: Could not read translation matrix.",
              "In_Configuration", io::message::error);
      return false;
    }
  }
  if (i != 3) {
    io::messages.add("ROTTRANSREFPOS block: Could not read translation matrix.",
            "In_Configuration", io::message::error);
    return false;
  }
  
  for (i = 0; it != to && i < 3; ++it, ++i) {
    _lineStream.clear();
    _lineStream.str(*it);

    for(unsigned int j = 0; j < 3; ++j)
      _lineStream >> rottrans.theta_inv_rot(i,j);

    if (_lineStream.fail()) {
      io::messages.add("ROTTRANSREFPOS block: Could not read rotation matrix.",
              "In_Configuration", io::message::error);
      return false;
    }
  }
  if (i != 3) {
    io::messages.add("ROTTRANSREFPOS block: Could not read rotation matrix.",
            "In_Configuration", io::message::error);
    return false;
  }
  
  // make sure there is enough space for the positions
  rottrans.pos.resize(last);
  
  for (i = 0; it != to && i < last; ++it, ++i) {
    _lineStream.clear();
    _lineStream.str(*it);

    for(unsigned int j = 0; j < 3; ++j)
      _lineStream >> rottrans.pos(i)(j);

    if (_lineStream.fail()) {
      io::messages.add("ROTTRANSREFPOS block: Could not read reference positions.",
              "In_Configuration", io::message::error);
      return false;
    }
  }  
  
  if (it != to) {
    io::messages.add("ROTTRANSREFPOS block: Too many lines (reference positions)",
            "In_Configuration", io::message::error);
    return false;
  }
  
  if (i != last) {
    io::messages.add("ROTTRANSREFPOS block: Not enough lines (reference positions)",
            "In_Configuration", io::message::error);
    return false;
  }  

  return true;    
}

bool io::In_Configuration::
_read_xray_av(std::vector<std::string> &buffer,
        std::vector<configuration::Configuration::special_struct::xray_struct> & xray_av,
        std::vector<topology::xray_restraint_struct> const & xray_res) {
  DEBUG(8, "read xray averages");

  // no title in buffer!
  std::vector<std::string>::const_iterator it = buffer.begin(),
          to = buffer.end() - 1;

  std::vector<topology::xray_restraint_struct>::const_iterator
  xray_it = xray_res.begin(),
          xray_to = xray_res.end();

  xray_av.clear();

  double av, phase_av;

  if (buffer.size() - 1 != xray_res.size()) {
    std::cout << "XRAYRESEXPAVE: " << buffer.size() - 1
            << " but restraints: " << xray_res.size()
            << std::endl;

    io::messages.add("number of Xray-restraints does not match with number of "
            "continuation data", "in_configuration",
            io::message::error);
    return false;
  }

  for (; (it != to) && (xray_it != xray_to); ++it, ++xray_it) {

    _lineStream.clear();
    _lineStream.str(*it);

    _lineStream >> av >> phase_av;
    if (_lineStream.fail()) {
      io::messages.add("Bad value in XRAYRESEXPAVE block",
              "In_Configuration", io::message::error);
      return false;
    }
    configuration::Configuration::special_struct::xray_struct tempstruct = {av, 0.0, phase_av, 0.0};
    xray_av.push_back(tempstruct);
  }

  if (xray_it != xray_to || it != to) {
    io::messages.add("Wrong number of Xray-Av's in XRAYRESEXPAVE block",
            "In_Configuration",
            io::message::error);
    return false;
  }

  return true;
}

