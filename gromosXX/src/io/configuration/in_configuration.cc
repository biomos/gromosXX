/**
 * @file in_configuration.cc
 * implements methods of In_Configuration.
 */

#include <util/stdheader.h>

#include <topology/core/core.h>

#include <topology/solute.h>
#include <topology/solvent.h>
#include <topology/perturbed_atom.h>
#include <topology/perturbed_solute.h>

#include <topology/topology.h>
#include <simulation/multibath.h>
#include <simulation/parameter.h>
#include <configuration/energy.h>
#include <configuration/energy_average.h>
#include <configuration/configuration.h>

#include <io/blockinput.h>
#include <io/instream.h>
#include <io/configuration/inframe.h>

#include <util/generate_velocities.h>
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
				topology::Topology const &topo, 
				simulation::Parameter const &param)
{
  if (!quiet)
    std::cout << "\nCONFIGURATION\n";
  
  // resize the configuration
  conf.resize(topo.num_atoms());

  // possibly resize the dihedral angle monitoring array
  // initialize or set to such a value that it is recalculated in
  // the first step, initialization would require an additional function
  // which can be done only after reading of the coordinates. Would be a
  // nicer solution, but also requires the parameters...
  if(param.print.monitor_dihedrals)
    conf.special().dihedral_angle_minimum.resize
      (topo.solute().dihedrals().size(), 4*math::Pi);
    
  DEBUG(8, "reading in a frame");
  read_frame();
  DEBUG(8, "frame read");

  block_read.clear();
  
  std::vector<std::string> buffer;

  // read positions
  buffer = m_block["POSITION"];
  if (buffer.size()){
    if (!quiet)
      std::cout << "\treading POSITION...\n";
    _read_position(conf.current().pos, buffer, topo.num_atoms());
    block_read.insert("POSITION");
  }
  else{
    buffer = m_block["POSITIONRED"];
    if (buffer.size()){
      if (!quiet)
	std::cout << "\treading POSITIONRED...\n";
      _read_positionred(conf.current().pos, buffer, topo.num_atoms());
      block_read.insert("POSITIONRED");
    }
    else{
      io::messages.add("no POSITION / POSITIONRED block found in input configuration",
		       "in_configuration",
		       io::message::error);
    }
  }
  
  // read velocities
  if(!param.start.generate_velocities){
    buffer = m_block["VELOCITY"];
    if (buffer.size()){
      if (!quiet)
	std::cout << "\treading VELOCITY...\n";
      _read_velocity(conf.current().vel, buffer, topo.num_atoms());
      block_read.insert("VELOCITY");
    }
    else{
      buffer = m_block["VELOCITYRED"];
      if (buffer.size()){
	if (!quiet)
	  std::cout << "\treading VELOCITYRED...\n";
	_read_velocityred(conf.current().vel, buffer, topo.num_atoms());
	block_read.insert("VELOCITYRED");
      }
      else{
	io::messages.add("no VELOCITY / VELOCITYRED block found in input configuration",
			 "in_configuration",
			 io::message::error);
	conf.current().vel = 0.0;
      }
    }
    // store also in old velocities (for initial temperature calculation)
    conf.old().vel = conf.current().vel;
  }
  else{
    // generate initial velocities
    util::generate_velocities(param.start.tempi, 
			      topo.mass(),
			      conf.current().vel,
			      conf.old().vel,
			      param.start.ig);
  }
  
  // read box
  if(param.boundary.boundary != math::vacuum){
    buffer = m_block["TRICLINICBOX"];
    if (buffer.size()){
      if (!quiet)
	std::cout << "\treading TRICLINICBOX...\n";
      _read_box(conf.current().box, buffer, param.boundary.boundary);
      conf.old().box = conf.current().box;
      block_read.insert("TRICLINICBOX");
    }
    else{
      buffer = m_block["BOX"];
      if (buffer.size() && param.boundary.boundary == math::rectangular){
	if (!quiet)
	  std::cout << "\treading BOX...\n";
	_read_g96_box(conf.current().box, buffer);
	conf.old().box = conf.current().box;
	block_read.insert("BOX");
      }
      else{
	io::messages.add("no TRICLINICBOX / BOX (for rectangular boundary conditions)\n"
			 "\tblock found in input configuration",
			 "in_configuration",
			 io::message::error);
      }
    }    
  }
  // and set the boundary type!
  conf.boundary_type = param.boundary.boundary;

  // resize the energy arrays
  const size_t num = topo.energy_groups().size();
  const size_t numb = param.multibath.multibath.size();

  DEBUG(5, "number of energy groups: " << num 
	<< "\nnumber of baths: " << numb);

  conf.current().energies.resize(num, numb);
  conf.current().energy_averages.resize(num, numb);
  conf.current().perturbed_energy_derivatives.resize(num, numb);
  conf.current().perturbed_energy_derivative_averages.resize(num, numb);

  conf.old().energies.resize(num, numb);
  conf.old().energy_averages.resize(num, numb);
  conf.old().perturbed_energy_derivatives.resize(num, numb);
  conf.old().perturbed_energy_derivative_averages.resize(num, numb);

  // resize some special data
  conf.special().rel_mol_com_pos.resize(topo.num_atoms());

  // print some information
  if (!quiet){
    std::cout << "\n\t";
    switch(conf.boundary_type){
      case math::vacuum:
	std::cout << "PBC            = vacuum\n";
	break;
      case math::rectangular:
	std::cout << "PBC            = rectangular\n";
	break;
      case math::triclinic:
	std::cout << "PBC            = triclinic\n";
	break;
      default:
	std::cout << "wrong periodic boundary conditions!";
	io::messages.add("wrong PBC!", "In_Configuration", io::message::error);
    }
  }

  if (!quiet){
    std::cout << "\ttotal mass     = " << math::sum(topo.mass()) << "\n"
	      << "\tvolume         = " << math::volume(conf.current().box, conf.boundary_type);
    
    if (conf.boundary_type != math::vacuum)
      std::cout << "\n\tdensity        = " 
		<< math::sum(topo.mass()) / math::volume(conf.current().box,
							 conf.boundary_type);
    
    std::cout << "\n\n";
  }
  
  if (param.constraint.solute.algorithm == simulation::constr_flexshake){
    conf.special().flexible_vel.resize(topo.solute().distance_constraints().size()+
				       topo.perturbed_solute().distance_constraints().size());
    conf.special().flexible_ekin.resize(numb);

    buffer = m_block["FLEXV"];
    if (buffer.size() && param.constraint.solute.flexshake_readin){
      block_read.insert("FLEXV");
      if (!quiet)
	std::cout << "\treading FLEXV...\n";
      _read_flexv(conf.special().flexible_vel, buffer, 
		  topo.solute().distance_constraints(),
		  topo.perturbed_solute().distance_constraints());

    }
    else{
      if (param.constraint.solute.flexshake_readin)
	io::messages.add("no FLEXV block found but reading in of constraint velocities requested",
			 "in_configuration",
			 io::message::error);
      else
	io::messages.add("no FLEXV block found, assuming SHAKE'n positions (and velocities)",
			 "in_configuration",
			 io::message::notice);
    }
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
    std::cout << "END\n\n";

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
      throw std::runtime_error("bad line in POSITIONRED block");
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
      throw std::runtime_error("bad line in POSITION block");
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
      throw std::runtime_error("bad line in VELOCITYRED block");
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
      
      throw std::runtime_error("bad line in VELOCITY block");
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
      io::messages.add("bad line in BOX block","In_Configuration", io::message::error);
      throw std::runtime_error("bad line in TRICLINICBOX block");
    }
    
    _lineStream.clear();
    _lineStream.str(*it);

    _lineStream >> box(0)(i) >> box(1)(i) >> box(2)(i);
    
    if(_lineStream.fail())
      throw std::runtime_error("bad line in TRICLINICBOX block");
    if (!_lineStream.eof()) {
      std::string msg = "Warning, end of line not reached, but should have been: \n" + *it +  "\n";
      DEBUG(10, msg);
    }
  }
  
  // and check the boundary condition...
  if (bound == 0 && boundary != math::vacuum){
    io::messages.add("Boundary condition from input file and from TRICLINICBOX do not match!",
		     "In_Configuration", io::message::warning);
  }
  else if (bound > 0 && boundary != math::rectangular){
    io::messages.add("Boundary condition from input file and from TRICLINICBOX do not match!",
		     "In_Configuration", io::message::warning);
  }
  else if (bound < 0 && boundary != math::triclinic){
    io::messages.add("Boundary condition from input file and from TRICLINICBOX do not match!",
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
    
  if(_lineStream.fail())
    throw std::runtime_error("bad line in TRICLINICBOX block");
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
  double v;
  
  for(c=0; (it != to) && (constr_it != constr_to); ++it, ++constr_it, ++flexv_it, ++c){

    _lineStream.clear();
    _lineStream.str(*it);

    _lineStream >> i >> j >> v;
    
    --i;
    --j;
    
    if(_lineStream.fail())
      throw std::runtime_error("bad line in FLEXV block");
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
    
    if(_lineStream.fail())
      throw std::runtime_error("bad line in FLEXV block");

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
