/**
 * @file in_jvalue.cc
 * implements methods of In_Jvalue
 */

#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../interaction/interaction_types.h"
#include "../../configuration/configuration.h"

#include "../../io/instream.h"
#include "../../io/blockinput.h"

#include "in_jvalue.h"

#undef MODULE
#undef SUBMODULE
#define MODULE io
#define SUBMODULE topology

/**
 * @section jvalresspec JVALRESSPEC block
 * The JVALRESSPEC block is read from the @f$^3J@f$-value restraint specification file.
 *
 * - \c IPJV, \c JPJV, \c KPJV and \c LPJV are atom sequence numbers 
 *   of the real atoms present in the simulation that define the dihedral angle 
 *   @f$\phi@f$ related to the restrained @f$^3J@f$-value @f$\phi@f$.
 * - \c WJVR is an individual @f$^3J@f$-value restraint weight factor by which
 *   the restraining term for each @f$^3J@f$-value may be multiplied
 * - \c PJR0 is the experimental or reference @f$^3J@f$-value, @f$J_0@f$. In case of a 
 *   full-harmonic @f$^3J@f$-value restraint (\c NHJV = 0), it is the minimum-energy @f$^3J@f$-value; 
 *   in the case of an attractive or repulsive half-harmonic @f$^3J@f$-value restraint
 *   (\c NHJV = @f$\pm@f$ 1), it is the upper or lower bound, respectively, beyond which the 
 *   restraining force becomes non-zero.
 * - \c PSJR is the phase shift or difference @f$\delta = \theta - \phi@f$
 *   between the dihedral angle @f$\theta@f$ formed by the possibly non-existant
 *   atoms defining the experimental @f$^3J@f$-coupling and the dihedral angle @f$\phi(i-j-k-l)@f$ 
 *   formed by the real atoms present in the simulation (in degrees).
 * - \c AJV, \c BJV and \c CJV are the Karplus parameters @f$a@f$, @f$b@f$ and
 *   @f$c@f$ for the @f$^3J@f$-coupling constant expressed as a function of the dihedral angle
 *   @f$\theta@f$.
 * - \c NHJV is the type of the @f$^3J@f$-value restraint:
 *  - -1: half-harmonic repulsive
 *  -  0: full harmonic [recommended]
 *  -  1: half-harmonic attractive
 *
 *   Note that the half-harmonic forms of the potential are only implemented in analogy to
 *   distance restraining and make little sense for restraining @f$^3J@f$-values, which depend on
 *   a periodic structural parameter. 
 *
 * @verbatim
JVALRESSPEC
# For each J-coupling constant restraint the following is to be specified:
# IPJV, JPJV, KPJV, LPJV            atom sequence numbers
# WJVR                              individual weight factor
# PJR0                              (>=0) reference J-value
# PSJR                              phase shift
# AJV,BJV,CJV                       Karplus parameters
# NHJV                     -1,0,1   type of J-value restraint.
#
#IPJV JPJV KPJV LPJV       WJVR   PJR0 PSJR  AJV    BJV   CJV    NHJV
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
    
    int i = 0, j = 0, k = 0, l = 0;
    double K = 0.0, J = 0.0;
    double a = 0.0, b = 0.0, c = 0.0, delta = 0.0;
    int H = 0;
    
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
	io::messages.add("bad value for NHJV in JVALRESSPEC block "
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
