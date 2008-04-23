/**
 * @file in_friction.cc
 * implements methods of In_Friction.
 */


#include <stdheader.h>

#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <simulation/simulation.h>
#include <interaction/interaction_types.h>
#include <configuration/configuration.h>

#include <io/instream.h>
#include <io/blockinput.h>

#include "in_friction.h"

#undef MODULE
#undef SUBMODULE
#define MODULE io
#define SUBMODULE topology


void 
io::In_Friction::read(topology::Topology& topo,
		    simulation::Simulation & sim,
		    std::ostream & os){
  
  DEBUG(7, "reading in a friciton specification file");
  
  std::vector<std::string> buffer;
  std::vector<std::string>::const_iterator it;
  
  { // FRICTIONSPEC

    buffer = m_block["FRICTIONSPEC"];
    DEBUG(10, "FRICTIONSPEC block : " << buffer.size());
    
    if (!buffer.size()){
      io::messages.add("no FRICTIONSPEC block in friction specification file",
		       "In_Friction", io::message::error);
      return;
    }

    std::vector<std::string>::const_iterator it = buffer.begin()+1,
      to = buffer.end()-1;

    DEBUG(10, "reading in FRICTIONSPEC data");
    unsigned int i;
    for(i=0; it != to; ++i, ++it){
      _lineStream.clear();
      std::string line = *it;
      if (line.length() > 24) {
        line.erase(line.begin(), line.begin() + 24); // erase the first 24 chars
        _lineStream.str(line);
      } else {
	io::messages.add("bad line in FRICTIONSPEC block",
			 "In_Friction",
			 io::message::error);     
        return;
      }

      double gam0;
      _lineStream >> gam0;
    
      if(_lineStream.fail()){
	io::messages.add("bad line in FRICTIONSPEC block",
			 "In_Friction",
			 io::message::error);
        return;
      }

      if (i >= topo.num_atoms()) {
        io::messages.add("too many lines in FRICTIONSPEC block", "In_Friction",
                         io::message::error);
        return;
      }
      topo.stochastic().gamma(i) = gam0;
    }
    
    if (i < topo.num_atoms()) {
      io::messages.add("not enough lines in FRICTIONSPEC block", "In_Friction",
                       io::message::error);
        return;
    }
    
  } // FRICTIONSPEC
  
  if (!quiet){
    os << "FRICTION SPECIFICATION" << std::endl
       << std::setw(10) << "SEQ"
       << std::setw(20) << "GAMMA0" << std::endl;
    for(unsigned int i = 0; i < topo.num_atoms(); ++i) {
      os << std::setw(10) << i+1
         << std::setw(20) << topo.stochastic().gamma(i) << std::endl;
    }
    os << "END" << std::endl;
  }
  
}
