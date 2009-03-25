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

/**
 * @section jvalresspec JVALRESSPEC block
 * The JVALRESSPEC block is read from the J-value restraint specification file.
 *
 * - Variables \c IPJR, \c JPJR, \c KPJR and \c LPJR are atom sequence numbers 
 *   of the real atoms defining the dihedral angle @f$\phi@f$ that is related 
 *   to the restrained J-value.
 * - Variable \c WJVR is an individual J-value restraint weight factor by which
 *   the J-Value restraining term may be multiplied
 * - Variable \c PJR0 or @f$J_0 \geq 0@f$. In case of a full-harmonic J-value
 *   restraint (\c H = 0), \c PJR0 is the minimum-energy J-value; in case of an 
 *   attractive or repulsive half-harmonic J-value restraint (\c H = +- 1),
 *   \c PJR0 is the upper or lower bound, respectively, beyond which the 
 *   restraining force becomes non-zero.
 * - Variable \c PSJR is the phase shift or difference between the dihedral
 *   angle @f$\theta@f$ formed by the possibly non-existing H-atoms defining the
 *   J-coupling and the dihedral angle i-j-k-l @f$\phi@f$ formed by the real 
 *   atoms that is related to the J-coupling (in degrees). Phase shift
 *   @f$\delta = \theta - \phi@f$.
 * - Variables \c A, \c B and \c C are Karplus parameters @f$a@f$, @f$b@f$ and
 *   @f$c@f$ for the J-coupling constant expressed as a function of
 *   @f$\theta@f$.
 * - Variable \c H is the type of the J-Value restraint:
 *  - -1: half-harmonic repulsive
 *  - 0: full harmonic
 *  - 1: half_harmonic attractive
 *
 * @verbatim
JVALRESSPEC
# For each J-coupling constant restraint the following is to be specified:
# IPJR, JPJR, KPJR, LPJR            atom sequence numbers
# WJVR                              weight factor
# PJR0                              >=0. zero value
# PSJR                              phase shift
# A,B,C                             Karplus parameters
# H                        -1,0,1   type of J-value restraint.
#
#IPJR JPJR KPJR LPJR       WJVR   PJR0 PSJR  A      B     C      H
# 1VAL
  109    1    3    7       10.0   7.3  -60   9.4    -1.1  0.4    0
# 4ALA
   21   23   25   27       10.0   8.6  -60   9.4    -1.1  0.4    0
# 5PHE
   27   29   31   44       10.0   6.8  -60   9.4    -1.1  0.4    0
# 6PHE
   44   46   48   61       10.0   6.6  -60   9.4    -1.1  0.4    0
# 9PHE
   75   77   79   92       10.0   8.3  -60   9.4    -1.1  0.4    0
# 10PHE
   92   94   96  109       10.0   6.7  -60   9.4    -1.1  0.4    0
END
@endverbatim
 */
void 
io::In_Jvalue::read(topology::Topology& topo,
		    configuration::Configuration & conf,
		    simulation::Simulation & sim,
		    std::ostream & os){
  
  DEBUG(7, "reading in a jvalue restraints specification file");

  if (!quiet) os << "JVALUE RESTRAINTS\n";
  
  std::vector<std::string> buffer;
  
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
      /*
       * full harmonic (attractive and repulsive)
        harmonic = 0,
       * half harmonic, attractive
        attractive = 1,
       * half harmonic, repulsive
       repulsive = -1
       */
      if (H < -1 || H > 1){
	io::messages.add("bad value for H in JVALRESSPEC block "
			 "(harmonic: 0, attractive: 1, repulsive: -1)",
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
    if (sim.param().jvalue.le)
      conf.special().jvalue_epsilon.resize(conf.special().jvalue_av.size());

    if (!quiet) os << "END\n";

  } // JVALRESSPEC
    

  
}
