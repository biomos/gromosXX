/**
 * @file in_jvalue.cc
 * implements methods of In_Jvalue
 */


#include <stdheader.h>

#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <simulation/simulation.h>
#include <interaction/interaction_types.h>
#include <configuration/configuration.h>

#include <io/instream.h>
#include <io/blockinput.h>

#include "in_jvalue.h"

#undef MODULE
#undef SUBMODULE
#define MODULE io
#define SUBMODULE topology

void 
io::In_Jvalue::read(topology::Topology& topo,
		    configuration::Configuration & conf,
		    simulation::Simulation & sim,
		    std::ostream & os){
  
  DEBUG(7, "reading in a jvalue restraints specification file");

  if (!quiet) os << "JVALUE RESTRAINTS\n";
  
  std::vector<std::string> buffer;
  std::vector<std::string>::const_iterator it;
  
  { // JVALUERES
    DEBUG(10, "JVALRESSPEC block");
    buffer = m_block["JVALRESSPEC"];
    
    if (!buffer.size()){
      io::messages.add("no JVALRESSPEC block in jvalue restraints file",
		       "in_jvalue", io::message::error);
      return;
    }

    std::vector<std::string>::const_iterator it = buffer.begin()+1,
      to = buffer.end()-1;
  
    DEBUG(10, "reading in JVALRESSPEC data");
    DEBUG(10, "jvalue_av size = " << conf.special().jvalue_av.size());
    
    int i, j, k, l;
    double K, J;
    double a, b, c, delta;
    int H;
    
    os.precision(2);
    os.setf(std::ios_base::fixed, std::ios_base::floatfield);

    if (!quiet){
      os << std::setw(6) << "i"
	 << std::setw(6) << "j"
	 << std::setw(6) << "k"
	 << std::setw(6) << "l"
	 << std::setw(8) << "K"
	 << std::setw(8) << "J"
	 << std::setw(8) << "a"
	 << std::setw(8) << "b"
	 << std::setw(8) << "c"
	 << std::setw(8) << "delta"
	 << std::setw(7) << "H"
	 << std::setw(8) << "Jav"
	 << "\n";
    }
    
    for(int n=0; it != to; ++i, ++it, ++n){
      
      _lineStream.clear();
      _lineStream.str(*it);

      _lineStream >> i >> j >> k >> l
		  >> K >> J
		  >> delta
		  >> a >> b >> c
		  >> H;

      if(_lineStream.fail()){
	io::messages.add("bad line in JVALRESSPEC block",
			 "In_Jvalue",
			 io::message::error);
      }

      // restr_func_enum:
      // 1 : harmonic
      // 2 : attractive
      // 3 : repulsive
      if (H < 1 || H > 3){
	io::messages.add("bad value for H in JVALRESSPEC block "
			 "(harmonic: 1, attractive: 2, repulsive: 3)",
			 "In_Jvalue",
			 io::message::error);
      }

      topo.jvalue_restraints().push_back
	(topology::jvalue_restraint_struct(i-1, j-1, k-1, l-1,
					   K, J,
					   a, b, c, math::Pi * delta / 180.0,
					   topology::functional_form(H)));

      if (conf.special().jvalue_av.size() < unsigned(n+1))
	conf.special().jvalue_av.push_back(J);

      if (!quiet){
	os << std::setw(6) << i
	   << std::setw(6) << j
	   << std::setw(6) << k
	   << std::setw(6) << l
	   << std::setw(8) << K
	   << std::setw(8) << J
	   << std::setw(8) << a
	   << std::setw(8) << b
	   << std::setw(8) << c
	   << std::setw(8) << delta
	   << std::setw(7) << H
	   << std::setw(8) << conf.special().jvalue_av[n]
	   << "\n";
      }
    }

    assert(conf.special().jvalue_av.size() == topo.jvalue_restraints().size());
    conf.special().jvalue_curr.resize(conf.special().jvalue_av.size());

    if (!quiet) os << "END\n";

  } // JVALRESSPEC
    

  
}
