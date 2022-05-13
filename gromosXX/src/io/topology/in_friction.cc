/**
 * @file in_friction.cc
 * implements methods of In_Friction.
 */


#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../interaction/interaction_types.h"
#include "../../configuration/configuration.h"

#include "../../io/instream.h"
#include "../../io/blockinput.h"

#include "in_friction.h"

#undef MODULE
#undef SUBMODULE
#define MODULE io
#define SUBMODULE topology

/**
 * @section frictionspec FRICTIONSPEC block
 * The FRICTIONSPEC block is read from the friction specification file.
 *
 * For every atom, one \c GAM0 friction coefficient is read.
 *
 * @verbatim
FRICTIONSPEC
# the first 24 characters are ignored
# RESIDUE   ATOM           GAM0
    1 HEXA  CH31       1   20.0
    1 HEXA  CH22       2   15.0
    1 HEXA  CH23       3   25.0
    1 HEXA  CH24       4   10.0
    2 HEXA  CH23       9   28.0
    2 HEXA  CH24      10   24.0
    2 HEXA  CH25      11   28.0
    2 HEXA  CH36      12   30.0
END
@endverbatim
 */
void 
io::In_Friction::read(topology::Topology& topo,
		    simulation::Simulation & sim,
		    std::ostream & os){
  
  DEBUG(7, "reading in a friciton specification file");
  
  std::vector<std::string> buffer;
  
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
    unsigned int i = 0;
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

      double gam0 = 0.0;
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
