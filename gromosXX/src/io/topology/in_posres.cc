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

#include "in_posres.h"

#undef MODULE
#undef SUBMODULE
#define MODULE io
#define SUBMODULE topology


void 
io::In_Posres::read(topology::Topology& topo,
		    configuration::Configuration & conf,
		    simulation::Simulation & sim,
		    std::ostream & os){
  
  DEBUG(7, "reading in a position restraints file");

  os << "POSITION RESTRAINTS\n";
  
  std::vector<std::string> buffer;
  std::vector<std::string>::const_iterator it;
  
  { // POSRES

    buffer = m_block["POSRES"];
    DEBUG(10, "POSRES block : " << buffer.size());
    
    if (!buffer.size()){
      io::messages.add("no POSRES block in position restraints file",
		       "in_posres", io::message::error);
      return;
    }

    std::vector<std::string>::const_iterator it = buffer.begin()+1,
      to = buffer.end()-1;
  
    std::string s1, s2;
    unsigned int i, n, nr;
    math::Vec pos;

    DEBUG(10, "reading in POSRES data");

    for(i=0; it != to; ++i, ++it){
      
      DEBUG(11, "\tnr " << i);
      
      _lineStream.clear();
      _lineStream.str(*it);

      _lineStream >> n >> s1 >> s2 >> nr;
      _lineStream >> pos(0) >> pos(1) >> pos(2);

      DEBUG(11, "\t" << n << "\t" << s1 << "\t" << s2 << "\t" << nr);
    
      if(_lineStream.fail()){
	io::messages.add("bad line in POSRES block",
			 "In_Posres",
			 io::message::error);
      }

      if (nr < 1 || nr > topo.num_atoms()){
	io::messages.add("illegal atom in POSRES block",
			 "In_Posres",
			 io::message::error);
      }

      topo.position_restraints().push_back
	(topology::position_restraint_struct(nr-1, pos));

    }
    
  } // POSRES
    
  
  { // BFACTOR
    buffer = m_block["BFACTOR"];

    if (buffer.size()){
      DEBUG(10, "BFACTOR block");
      os << "\tBFACTOR";

      io::messages.add("reading in atomic bfactors for position restraints",
		       "in_posres", io::message::notice);
      
      it = buffer.begin() + 1;
      _lineStream.clear();
      _lineStream.str(*it);
      unsigned int n, num;
      
      _lineStream >> num;
      ++it;

      if (num != topo.position_restraints().size()){
	io::messages.add("illegal number of bfactors in posres file",
			 "in_posres", io::message::error);
	return;
      }
      
      std::vector<topology::position_restraint_struct>::iterator
	posres_it = topo.position_restraints().begin();
      
      for(n=0; it != buffer.end() - 1; ++it, ++n, ++posres_it){
	unsigned int seq;
	double bfactor;
	
	_lineStream.clear();
	_lineStream.str(*it);
	_lineStream >> seq >> bfactor;
	
	if (_lineStream.fail() || ! _lineStream.eof()){
	  io::messages.add("Bad line in BFACTOR block",
			   "In_Posres", io::message::error);
	}
      
	if (seq > topo.num_atoms() || seq < 1 || seq != posres_it->seq){
	  io::messages.add("Atom number out of range / wrong in BFACTOR block",
			   "In_Posres", io::message::error);
	}
      
	posres_it->bfactor = bfactor;
	
      }
    
      if(n != num){
	io::messages.add("Wrong number of bfactors in BFACTOR block",
			 "In_Posres", io::message::error);
      }
    }
  
  } // BFACTOR

  std::vector<topology::position_restraint_struct>::const_iterator
    posres_it = topo.position_restraints().begin(),
    posres_to = topo.position_restraints().end();
  
  switch(sim.param().posrest.posrest){
    case 0:
      os << "\tPosition restraints OFF\n";
      // how did you get here?
      break;
    case 1:
      os << "\tPosition restraints ON\n"
	 << "\t\trestraining to following positions:\n";
      break;
    case 2:
      os << "\tPosition restraints ON\n"
	 << "\t\trestraining to following positions\n"
	 << "\t\t(using atomic B-factors):\n";
      break;
    case 3:
      os << "\tPosition constraints ON\n"
	 << "\t\tconstraining following atoms"
	 << " to their initial positions:\n";
      break;
  }

  os << std::setw(10) << "SEQ"
     << std::setw(20) << "POS(X)"
     << std::setw(20) << "POS(Y)"
     << std::setw(20) << "POS(Z)"
     << std::setw(15) << "BFACTOR"
     << "\n";
  
  for( ; posres_it != posres_to; ++posres_it)
    os << std::setw(10) << posres_it->seq
       << std::setw(20) << posres_it->pos(0)
       << std::setw(20) << posres_it->pos(1)
       << std::setw(20) << posres_it->pos(2)
       << std::setw(15) << posres_it->bfactor
       << "\n";
  os << "END\n";
  
}
