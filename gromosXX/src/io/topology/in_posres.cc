/**
 * @file in_posres.cc
 * implements methods of In_Topology.
 */


#include <stdheader.h>

#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <simulation/simulation.h>
#include <interaction/interaction_types.h>
#include <configuration/configuration.h>

#include <io/instream.h>
#include <io/blockinput.h>
#include <io/configuration/in_configuration.h>
#include <io/configuration/out_configuration.h>

#include "in_posres.h"

#undef MODULE
#undef SUBMODULE
#define MODULE io
#define SUBMODULE topology


void 
io::In_Posres::read(topology::Topology& topo,
		    simulation::Simulation & sim,
		    std::ostream & os){
  
  DEBUG(7, "reading in a position restraints file");

  if (!quiet)
    os << "POSITION RESTRAINTS\n";
  
  std::vector<std::string> buffer;
  std::vector<std::string>::const_iterator it;

  { // POSRESSPEC

    buffer = m_block["POSRESSPEC"];
    DEBUG(10, "POSRESSPEC block : " << buffer.size());
    
    if (!buffer.size()){
      io::messages.add("no POSRESSPEC block in position restraints file",
		       "in_posres", io::message::error);
      return;
    }

    std::vector<std::string>::const_iterator it = buffer.begin()+1,
      to = buffer.end()-1;

    DEBUG(10, "reading in POSRESSPEC data");
    
    unsigned int num_posres;
    _lineStream.clear();
    _lineStream.str(*it);    
    _lineStream >> num_posres;
    if (_lineStream.fail()) {
       io::messages.add("Could not read number of position restraints in "
                        "POSRESSPEC block.", "in_posres", io::message::error);
      return;     
    }

    // jump to next line
    ++it;

    unsigned int i, nr;
    for(i=0; it != to; ++i, ++it){
      
      DEBUG(11, "\tnr " << i);
      
      std::string line(*it);
      if (line.length() < 17) {
        io::messages.add("line too short in POSRESSPEC block", "In_Posres",
                io::message::error);
      }
      
      // the first 17 chars are ignored
      line.erase(line.begin(), line.begin() + 17);
      
      _lineStream.clear();
      _lineStream.str(line);

      _lineStream >> nr;
      
      DEBUG(11, "\t" << nr);
      
      if(_lineStream.fail()){
        io::messages.add("bad line in POSRESSPEC block",
                "In_Posres", io::message::error);
        return;
      }
      
      if (nr < 1 || nr > topo.num_atoms()){
        io::messages.add("illegal atom in POSRESSPEC block",
                "In_Posres", io::message::error);
        return;
      }
      
      topo.position_restraints().push_back
      (topology::position_restraint_struct(nr-1, math::Vec(0.0)));
      
    }

    if (i != num_posres) {
      io::messages.add("number of position restraints and lines in POSRESSPEC block "
      "do not correspond.", "In_Posres", io::message::error);
      return;
    }
    
  } // POSRESSPEC  
  
  // read the reference positions from this file
  { // REFPOSITION
    buffer = m_block["REFPOSITION"];
    if (sim.param().posrest.read) {
      if (!buffer.size()) {
        io::messages.add("no REFPOSITION block in position restraints "
                         "specification file", "In_Posres", io::message::error);
        return;
      } else {
        io::In_Configuration::_read_refposition(topo.position_restraints(), 
                                                buffer, true);
      }
    } else {
      if (buffer.size()) {
        io::messages.add("REFPOSITION block ignored in position restraints "
                         "specification file", "In_Posres", io::message::warning);
      }      
    }
  } // REFPOSITION
  
  { // BFACTOR
    buffer = m_block["BFACTOR"];
    if (sim.param().posrest.read && 
        sim.param().posrest.posrest == simulation::posrest_bfactor) {
      if (!buffer.size()) {
        io::messages.add("no BFACTOR block in position restraints "
                         "specification file", "In_Posres", io::message::error);
        return;
      } else {
        io::In_Configuration::_read_bfactor(topo.position_restraints(),
                                            buffer, true);
      }
    } else {
      if (buffer.size()) {
        io::messages.add("BFACTOR block ignored in position restraints "
                         "specification file", "In_Posres", io::message::warning);
      }      
    }  
  } // BFACTOR

  std::vector<topology::position_restraint_struct>::const_iterator
    posres_it = topo.position_restraints().begin(),
    posres_to = topo.position_restraints().end();
  
  if (!quiet){
    switch(sim.param().posrest.posrest){
      case simulation::posrest_off :
	os << "\tPosition restraints OFF\n";
	// how did you get here?
	break;
      case simulation::posrest_on:
      case simulation::posrest_bfactor:
	os << "\tPosition restraints ON\n"
	   << "\t\trestraining to ";
        if (sim.param().posrest.read) 
          os <<  "the following positions:\n";
        else
	  os << "their positions in the startup file.\n";
	break;
      case simulation::posrest_const:
	os << "\tPosition constraints ON\n"
	   << "\t\tconstraining following atoms to ";
        if (sim.param().posrest.read) 
          os <<  "the following positions:\n";
        else
	  os << "their positions in the startup file.\n";
	break;
    }
    
    
    if (sim.param().posrest.read) {
      os << std::setw(10) << "SEQ"
         << std::setw(20) << "POS(X)"
         << std::setw(20) << "POS(Y)"
         << std::setw(20) << "POS(Z)"
         << std::setw(15) << "BFACTOR"
         << "\n";
      
      for(;posres_it != posres_to; ++posres_it)
        os << std::setw(10) << posres_it->seq
           << std::setw(20) << posres_it->pos(0)
           << std::setw(20) << posres_it->pos(1)
           << std::setw(20) << posres_it->pos(2)
           << std::setw(15) << posres_it->bfactor
           << "\n";
    }
    
    if (sim.param().posrest.scale_reference_positions)
      os << "\n\tReference positions are scaled upon pressure scaling.\n";
    
    os << "END\n";
  }
  
}
