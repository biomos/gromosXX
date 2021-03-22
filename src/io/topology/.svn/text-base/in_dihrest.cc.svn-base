/**
 * @file in_dihrest.cc
 * implements methods of In_Dihrest
 */


#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../interaction/interaction_types.h"
#include "../../configuration/configuration.h"

#include "../../io/instream.h"
#include "../../io/blockinput.h"

#include "in_dihrest.h"

#undef MODULE
#undef SUBMODULE
#define MODULE io
#define SUBMODULE topology

static std::set<std::string> block_read;

/**
 * @section dihedralresspec DIHEDRALRESSPEC block
 * This block is for unperturbed dihedrals restraints/constraints. It is read 
 * from the dihedral restraint specifcation file.
 *
 * @verbatim
DIHEDRALRESSPEC
# IPLR, JPLR, KPLR, LPLR   atom sequence numbers of the atoms defining the
#                          restrained dihedral i-j-k-l
# WDLR                     individual dihedral restraint weight factor by which
#                          the harmonic dihedral restraining term may be multiplied.
# PDLR                     dihedral angle value (in degrees) at minimum energy of
#                          of the harmonic dihedral restraining term.
# DELTA                    deviation of the zero-energy dihedral angle value after
#                          which the restraint switches to the next periodic
#                          value. The dihedral angle is put in the interval
#                          [ PDLR + DELTA - 360 , PDLR + DELTA ]
#  IPLR  JPLR  KPLR  LPLR  WDLR  PDLR  DELTA
    1     2     3     4    1.0   120.0 180.0
END
@endverbatim
 *
 * @section pertdihresspec PERTDIHRESTSPEC block
 * This block is for perturbed dihedral restraints/constraints. It is read from
 * the dihedral restraints specification file.
 *
 * @verbatim
PERTDIHRESSPEC
# IPLR, JPLR, KPLR, LPLR   atom sequence numbers of the atoms defining the
#                          restrained dihedral i-j-k-l
# APDLR    dihedral angle value (in degrees) at minimum energy of the harmonic
#          dihedral restraining term in state A.
# AWDLR    individual dihedral restraint weight factor by which the harmonic
#          dihedral restraining term may be multiplied in state A.
# BPDLR    as APDLR, but for state B.
# BWDLR    as AWDLR, but for state B.
# M        hidden restraint parameter m and
# N        hidden restraint parameter n of
#          hidden restraint prefactor l^n*(1-l)^m.
# DELTA    deviation of the zero-energy dihedral angle value after which the 
#          restraint switches to the next periodic value. The dihedral angle 
#          is put in the interval
#          [ (1-RLAM)*APDLR + RLAM*BPDLR + DELTA - 360 , 
#                                       (1-RLAM)*APDLR + RLAM*BPDLR + DELTA ]
#
# IPLR  JPLR  KPLR  LPLR  M   N  DELTA  APDLR  AWDLR  BPDLR  BWDLR
  1      2     3     4    2      180.0  120.0    1.0  160.0    1.0
END
@endverbatim
 */
void io::In_Dihrest::read(topology::Topology& topo,
			  simulation::Simulation & sim,
			  std::ostream & os){
  
  DEBUG(7, "reading in a dihedral restraints file");
  
  if (!quiet)
    os << "DIHEDRAL RESTRAINTS\n";
  
  std::vector<std::string> buffer;
  
  { // DIHEDRALRESSPEC
    DEBUG(10, "DIHEDRALRESSPEC block");
    buffer = m_block["DIHEDRALRESSPEC"];
    
    if (!buffer.size()){
      io::messages.add("no DIHEDRALRESSPEC block in dihedral restraints file",
		       "in_dihrest", io::message::error);
      return;
    }

    block_read.insert("DIHEDRALRESSPEC");

    std::vector<std::string>::const_iterator it = buffer.begin()+1,
      to = buffer.end()-1;

    int i, j, k, l;
    double delta, phi, w0;

    DEBUG(10, "reading in DIHREST data");

    if (!quiet){
      switch(sim.param().dihrest.dihrest){
	case simulation::dihedral_restr_off:
	  os << "\tDihedral restraints OFF\n";
	  // how did you get here?
	  break;
	case simulation::dihedral_restr_inst:
	  os << "\tDihedral restraints ON\n"
	     << "\t\t(uniform force constant K)\n";
	  break;
	case simulation::dihedral_restr_inst_weighted:
	  os << "\tDihedral restraints ON\n"
	     << "\t\t(force constant K*w0)\n";
	  break;
	case simulation::dihedral_constr:
	  os << "\tDihedral constraints ON\n";
	  break;
	default:
	  os << "\tDihedral restraints: ERROR\n";
	  io::messages.add("wrong value for method in dihedral restraints block",
			   "in_dihedral", io::message::error);
	  return;
      }
    }
    
    if (!quiet){

      os << std::setw(10) << "i"
	 << std::setw(8) << "j"
	 << std::setw(8) << "k"
	 << std::setw(8) << "l"
	 << std::setw(8) << "delta"
	 << std::setw(8) << "phi"
	 << std::setw(8) << "w0"
	 << "\n";

      os.precision(2);
      os.setf(std::ios::fixed, std::ios::floatfield);

    }
    
    for(int c=0; it != to; ++c, ++it){
      
      DEBUG(11, "\tnr " << c);
      
      _lineStream.clear();
      _lineStream.str(*it);

      _lineStream >> i >> j >> k >> l >> w0 >> phi >> delta;
      if(_lineStream.fail()){
	io::messages.add("bad line in DIHREST block",
			 "In_Dihrest", io::message::error);
      }
      
      topo.dihedral_restraints().push_back
	(topology::dihedral_restraint_struct(i-1, j-1, k-1, l-1,
					     delta * 2 * math::Pi / 360, phi * 2 * math::Pi / 360, w0));
      
      if (!quiet){
	os << std::setw(10) << i
	   << std::setw(8) << j
	   << std::setw(8) << k
	   << std::setw(8) << l
	   << std::setw(8) << delta
	   << std::setw(8) << phi
	   << std::setw(8) <<  w0
	   << "\n";
      }
    }
  } // DIHREST

  { // PERTDIHESPEC
    DEBUG(10, "PERTDIHRESSPEC block");
    buffer = m_block["PERTDIHRESSPEC"];
    
    block_read.insert("PERTDIHRESSPEC");
    if (!buffer.size()){
      return;
    }
    
    sim.param().perturbation.perturbed_par=true;
    std::vector<std::string>::const_iterator it = buffer.begin()+1,
      to = buffer.end()-1;

    int i, j, k, l, m, n;
    double delta, A_phi, A_w0, B_phi, B_w0;
    
    DEBUG(10, "reading in perturbed DIHREST (PERTDIHRESSPEC data");

    if (!quiet){
      switch(sim.param().dihrest.dihrest){
	case 0:
	  os << "\tPerturbed Dihedral restraints OFF\n";
	  // how did you get here?
	  break;
	case 1:
	  os << "\tPerturbed Dihedral restraints ON\n"
	     << "\t\t(using uniform force constant K\n";
	  break;
	case 2:
	  os << "\tPerturbed Dihedral restraints ON\n"
	     << "\t\t(using force constant K*w0)\n";
	  break;
	case 3:
	  os << "\tPerturbed Dihedral constraints ON\n";
	  break;
	default:
	  os << "\tPerturbed Dihedral restraints ERROR\n";
	  io::messages.add("wrong method for dihedral restraints",
			   "in_dihrest", io::message::error);
	  return;
      }
    }

    if (!quiet){
      os << std::setw(10) << "i"
	 << std::setw(8) << "j"
	 << std::setw(8) << "k"
	 << std::setw(8) << "l"
	 << std::setw(8) << "m"
	 << std::setw(8) << "n"
	 << std::setw(8) << "delta"
	 << std::setw(8) << "A_phi"
	 << std::setw(8) << "A_w0"
	 << std::setw(8) << "B_phi"
	 << std::setw(8) << "B_w0"
	 << "\n";

      os.precision(2);
      os.setf(std::ios::fixed, std::ios::floatfield);
    }
    
    for(int c=0; it != to; ++c, ++it){
      
      DEBUG(11, "\tnr " << c);
      
      _lineStream.clear();
      _lineStream.str(*it);

      _lineStream >> i >> j >> k >> l >> m >> n;
      _lineStream >> delta >> A_phi >> A_w0 >> B_phi >> B_w0;
    
      if(_lineStream.fail()){
	io::messages.add("bad line in PERTDIHREST block",
			 "In_Dihrest", io::message::error);
      }

      topo.perturbed_dihedral_restraints().push_back
	(topology::perturbed_dihedral_restraint_struct(i-1, j-1, k-1, l-1, m, n, delta * 2 * math::Pi / 360,
						       A_phi * 2 * math::Pi / 360, A_w0,
						       B_phi * 2 * math::Pi / 360, B_w0 ));

      if (!quiet){
	os << std::setw(10) << i
	   << std::setw(8) << j
	   << std::setw(8) << k
	   << std::setw(8) << l
	   << std::setw(8) << m
	   << std::setw(8) << n
	   << std::setw(8) << delta
	   << std::setw(8) << A_phi
	   << std::setw(8) << A_w0
	   << std::setw(8) << B_phi
	   << std::setw(8) << B_w0    
	   << "\n";
      }
      
    } // PERTDIHRESPEC
    
    if (!quiet) os << "END\n";
  }
  
  for(std::map<std::string, std::vector<std::string> >::const_iterator
	it = m_block.begin(),
	to = m_block.end();
      it != to;
      ++it){
    
    if (block_read.count(it->first) == 0 && it->second.size()){
      io::messages.add("block " + it->first + " not supported!",
		       "In_Dihrest",
		       io::message::warning);
    }
  }
}
