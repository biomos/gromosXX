/**
 * @file in_distanceres.cc
 * implements methods of In_Topology.
 */


#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../interaction/interaction_types.h"
#include "../../configuration/configuration.h"

#include "../../io/instream.h"
#include "../../io/blockinput.h"

#include "in_distanceres.h"

#undef MODULE
#undef SUBMODULE
#define MODULE io
#define SUBMODULE topology

static std::set<std::string> block_read;

/**
 * @section distanceresspec DISTANCERESSPEC block
 * The DISTANCERESSPEC block is read from the distance restraints specification
 * file.
 *
 * \c DISH is the carbon-hydrogen, \c DISC the carbon-carbon distance.
 * See @ref util::virtual_type for valid virtual atom types.
 * \c r0 is the restraint distance, \c w0 a weight factor (multiplied by the force 
 * constant specified in the input file, \c CDIR) and rah specifies the type of 
 * restraint. Possible values for \c rah
 * - -1: half harmonic repulsive 
 * - 0: full harmonic
 * - +1: half harmonic attractive
 *
 * @verbatim
DISTANCERESSPEC
# DISH  DISC
  0.1   0.153
# i  j  k  l  type    i  j  k  l  type    r0    w0    rah
  1  0  0  0  0       10 12 11 13 3       0.2   1.0   0
END
@endverbatim
 * @sa util::virtual_type util::Virtual_Atom
 *
 * @section pertdisresspec PERTDISRESSPEC block
 * The PERTDISRESSPEC block is read from the distance restraints specification
 * file and used for perturbed distance restraints and hidden 
 * restraints. 
 *
 * The format is very similar to the @ref distanceresspec with the difference that
 * one may give values for the A and the B state. The two variables \c n and
 * \c m are the parameters for the hidden restraints. 
 *
 * See:
 * - M. Christen, A.-P.E. Kunz, W.F. van Gunsteren, Sampling of rare events
 *   using hidden restraints, J. Phys. Chem. B 110 (2006) 8488-8498
 *
 * @verbatim
PERTDISRESSPEC
# DISH  DISC
  0.1   0.153
# i  j  k  l  type    i  j  k  l  type n m   A_r0  A_w0  B_r0   B_w0  rah
  1  0  0  0  0       10 12 11 13 3    1 1    0.2   1.0   0.5    2.0   0
END
@endverbatim
 *
 * @section mdisresspec MDISRESSPEC block
 * The MDISRESSPEC block is read from the distance restraints specification
 * file and used for EDS restraints.
 *
 * The format is very similar to the @ref distanceresspec with the difference that
 * one may give values for multipde states.
 * @verbatim
MDISRESSPEC
# DISH  DISC
  0.1   0.153
# N: number of eds states (3 in this example)
# i  j  k  l  type    i  j  k  l  type    r0[1 ... N]    w0[1 ... N]    rah
  1  0  0  0  0       10 12 11 13 3       0.2  0.2  0.2  1.0  0.0 0.0    0
  5  0  0  0  0       10 12 11 13 3       0.2  0.2  0.2  0.0  1.0 0.0    0
  8  0  0  0  0       10 12 11 13 3       0.2  0.2  0.2  0.0  0.0 1.0    0
END
@endverbatim
 *
 * @section dfresspec DFRESSPEC block
 * The DFRESSPEC block is read from the distance restraints specification
 * file and used for distancefield restraints.
 *
 * See:
 * - A. de Ruiter and C. Oostenbrink, Protein-ligand binding from distancefield
 *   distances and Hamiltonian replica exchange simulations, J. Chem. Theory 
 *   Comp. 9 (2013) 883 - 892, doi: 10.1021/ct300967a
 *
 *@verbatim
DFRESSPEC
#   DISH H-C bond length for virtual atoms
#   DISC C-C bond length for virtual atoms
#   PROTEINATOMS > 0 last atom of the host
#   K >= 0.0 Force constant
#   r0 >=0 zero energy distance 
#   TYPE_I Virtual atom type for interaction site I
#   NUM_I  Number of atoms defining interaction site I
#   ATOM_I[0..NUM_I] Index numbers of atoms defining interaction site I
#   TYPE_J Virtual atom type for interaction site J
#   NUM_J  Number of atoms defining interaction site J
#   ATOM_J[0..NUM_J] Index numbers of atoms defining interaction site J
# DISH  DISC
  0.1   0.153
# PROTEINATOMS  K    r0 
  1190          500  0.0 
# TYPE_I  NUM_I  ATOM_I[0] .. ATOM_I[NUM_I]
  -1      7        16  190  249  312  486  632 1208
# TYPE_J  NUM_J  ATOM_J[0] .. ATOM_J[NUM_J]
  -1      2      1194 1203
END
@endverbatim
 *
 * @section pertdfresspec PERTDFRESSPEC block
 * The PERTDFRESSPEC block is read from the distance restraints specification
 * file and used for perturbed distancefield restraints.
 *
 * See:
 * - A. de Ruiter and C. Oostenbrink, Protein-ligand binding from distancefield
 *   distances and Hamiltonian replica exchange simulations, J. Chem. Theory 
 *   Comp. 9 (2013) 883 - 892, doi: 10.1021/ct300967a
 *
 *@verbatim
PERTDFRESSPEC
#   DISH H-C bond length for virtual atoms
#   DISC C-C bond length for virtual atoms
#   PROTEINATOMS > 0 last atom of the host
#   A_r0 >=0 reference distance for state A
#   B_r0 >=0 reference distance for state B
#   K_A >= 0 force constant state A
#   K_B >= 0 force constant state B
#   n >= 0 hidden restraint parameter n 
#   m >= 0 hidden restraint parameter m 
#   TYPE_I Virtual atom type for interaction site I
#   NUM_I  Number of atoms defining interaction site I
#   ATOM_I[0..NUM_I] Index numbers of atoms defining interaction site I
#   TYPE_J Virtual atom type for interaction site J
#   NUM_J  Number of atoms defining interaction site J
#   ATOM_J[0..NUM_J] Index numbers of atoms defining interaction site J
# DISH  DISC
  0.1   0.153
# PROTEINATOMS  A_r0  K_A  B_r0  K_B  n  m
  1190          4.5   500  0.0   500  0  0
# TYPE_I  NUM_I  ATOM_I[0] .. ATOM_I[NUM_I]
  -1      7        16  190  249  312  486  632 1208
# TYPE_J  NUM_J  ATOM_J[0] .. ATOM_J[NUM_J]
  -1      2      1194 1203
END
@endverbatim
 */
void 
io::In_Distanceres::read(topology::Topology& topo,
		      simulation::Simulation & sim,
		      std::ostream & os){
  
  DEBUG(7, "reading in a distance restraints file");

  if (!quiet)
    os << "DISTANCE RESTRAINTS\n";
  
  std::vector<std::string> buffer;

  { // DISTANCERES
    DEBUG(10, "DISTANCERESSPEC block");
    buffer = m_block["DISTANCERESSPEC"];
    block_read.insert("DISTANCERESSPEC");
    if (buffer.size()<=2){
      io::messages.add("no or empty DISTANCERESSPEC block in distance restraints file",
		       "in_distanceres", io::message::warning);
    }
    else{      
      std::vector<std::string>::const_iterator it = buffer.begin()+1,
	to = buffer.end()-1;
      
      double dish,disc;
      bool nr_atoms=true;      
      
      DEBUG(10, "reading in DISTANCERES data");
      
      if (!quiet) {
        switch (sim.param().distanceres.distanceres) {
          case 0:
            os << "\tDistance restraints OFF\n";
            // how did you get here?
            break;
          case 1:
            os << "\tDistance restraints ON\n";
            break;
          case -1:
            os << "\tDistance restraints ON\n"
                    << "\ttime averaging ON\n";
            break;
          case 2:
            os << "\tDistance restraints ON\n"
                    << "\t\t(using force constant K*w0)\n";
            break;
          case -2:
            os << "\tDistance restraints ON\n"
                    << "\ttime averaging ON\n"
                    << "\t\t(using force constant K*w0)\n";
            break;
          default:
            os << "\tDistance restraints ERROR\n";
        }
      }

      
      _lineStream.clear();
      _lineStream.str(*it);
      
             
      _lineStream >> dish >> disc;
      
      if(_lineStream.fail()){
	std::ostringstream msg;
	msg << "bad line in DISTANCERESSPEC block: failed to read in DISH and DISC"  << std::endl;
	io::messages.add(msg.str(),
			 "In_Distanceres",
			 io::message::error);
      }
      
    
      ++it;
      
      if (!quiet){
	os << std::setw(10) << "DISH"
	   << std::setw(10) << "DISC"
	   << "\n" 
	   <<  std::setw(10)<< dish
	   <<  std::setw(10)<< disc
	   << "\n";
	
	os << std::setw(10) << "i"
	   << std::setw(8) << "j"
	   << std::setw(8) << "k"
	   << std::setw(8) << "l"
	   << std::setw(5) << "type"
	   << std::setw(10) << "i"
	   << std::setw(8) << "j"
	   << std::setw(8) << "k"
	   << std::setw(8) << "l"
	   << std::setw(5) << "type"
	   << std::setw(8) << "r0"
	   << std::setw(8) << "w0"
	   << std::setw(4) << "rah"
	   << "\n";
      }
    
      for(unsigned int line_number=2; it != to; ++line_number, ++it){
	
	DEBUG(11, "\tnr " << line_number - 2);
	
	int type1, type2;
	std::vector<int> atom1, atom2;
	double r0,w0;
	int rah;
	
	_lineStream.clear();
	_lineStream.str(*it);
	
	for(unsigned int i = 0; i < io::In_Distanceres::MAX_ATOMS; i++) {
	  int atom;
	  _lineStream >> atom;
	  // -1 because we directly convert to array indices
	  if (atom > 0) atom1.push_back(atom - 1);
          else if (atom <0){
            std::ostringstream msg;
            msg << "COM and COG type not possible for more than " 
                    << io::In_Distanceres::MAX_ATOMS << " atoms"  << std::endl;
            io::messages.add(msg.str(),
			 "In_Distanceres",
			 io::message::error);              
            nr_atoms=false;
          }
	}      
	_lineStream >> type1;
	
	for(unsigned int i = 0; i < io::In_Distanceres::MAX_ATOMS; i++) {
	  int atom;
	  _lineStream >> atom;
	  if (atom > 0) atom2.push_back(atom - 1);
          else if (atom <0){
            std::ostringstream msg;
            msg << "COM and COG type not possible for more than " 
                    << io::In_Distanceres::MAX_ATOMS << " atoms"  << std::endl;
            io::messages.add(msg.str(),
			 "In_Distanceres",
			 io::message::error);              
            nr_atoms=false;
          }
          
	}      
	_lineStream >> type2;
	
	_lineStream >> r0>> w0 >> rah;
	
	
	if(_lineStream.fail()){
	  std::ostringstream msg;
	  msg << "bad line in DISTANCERESSPEC block: " << line_number << std::endl
	      << "          " << *it;
	  io::messages.add(msg.str(),
			   "In_Distanceres",
			   io::message::error);
	}
	
	// g++ 3.2 fix
        if( nr_atoms){
	  util::virtual_type t1 = util::virtual_type(type1);
	  util::virtual_type t2 = util::virtual_type(type2);

	  util::Virtual_Atom v1(t1, atom1, dish, disc);
	  util::Virtual_Atom v2(t2, atom2, dish, disc);
        
	
	  topo.distance_restraints().push_back
	    (topology::distance_restraint_struct(v1,v2,r0,w0,rah));
	
	  if (!quiet){
	    for(unsigned int i = 0; i < io::In_Distanceres::MAX_ATOMS; i++) {
	      // the first element has the width 10, if i is bigger then the number of atoms
	      // specified, just print 0.
	      os << std::setw(i == 0 ? 10 : 8) << (i < atom1.size() ? atom1[i]+1 : 0);
	    }
	    os << std::setw(5) << type1;
	  
	    for(unsigned int i = 0; i < io::In_Distanceres::MAX_ATOMS; i++) {
	      os << std::setw(i == 0 ? 10 : 8) << (i < atom2.size() ? atom2[i]+1 : 0);
	    }
	    os << std::setw(5) << type2
	       << std::setw(8) << r0
	       << std::setw(8) << w0
	       << std::setw(4) << rah
	       << "\n";
          }
	}  
      }
    }
    
  } // DISTANCERES
  
  { // PERTDISRESPEC DISTANCERES
    DEBUG(10, "PERTDISRESSPEC block");
    buffer = m_block["PERTDISRESSPEC"];
    
    block_read.insert("PERTDISRESSPEC");

    // we don't have a PERTDISRESSPEC block
    if (!buffer.size()){
    }
    // check whether there is s.th. in the block
    else if (buffer.size()<=2){
      io::messages.add("empty PERTDISRESSPEC block in distance restraints file",
		       "in_distanceres", io::message::warning);
    }
    else{
      sim.param().perturbation.perturbed_par=true;
      std::vector<std::string>::const_iterator it = buffer.begin()+1,
	to = buffer.end()-1;
      
      double dish, disc;
      bool nr_atoms=true;
      
      DEBUG(10, "reading in DISTANCERES (PERTDISRESSPEC) data");
      
      if (!quiet){
	switch(sim.param().distanceres.distanceres){
	  case 0:
	    os << "\tPerturbed Distance restraints OFF\n";
	    // how did you get here?
	    break;
	  case 1:
	    os << "\tPerturbed Distance restraints ON\n";
	    break;
          case -1:
	    os << "\tPerturbed Distance restraints ON\n"
               << "\t\ttime-averaging ON\n";
	    break;
	  case 2:
	    os << "\tPerturbed Distance restraints ON\n"
	       << "\t\t(using force constant K*w0)\n";
	    break;
          case -2:
	    os << "\tPerturbed Distance restraints ON\n"
               << "\t\ttime-averaging ON\n"
	       << "\t\t(using force constant K*w0)\n";
	    break;
	  default:
	    os << "\tPerturbed Distance restraints ERROR\n";
	}
      }
      
      _lineStream.clear();
      _lineStream.str(*it);
      
      _lineStream >> dish >> disc ;
      
      if(_lineStream.fail()){
	std::ostringstream msg;
	msg << "bad line in PERTDISRESSPEC block: failed to read in DISH and DISC"  << std::endl;
	io::messages.add(msg.str(),
			 "In_Distanceres",
			 io::message::error);
      }
  
      ++it;
      if (!quiet){
	os << std::setw(10) << "DISH"
	   << std::setw(10) << "DISC"
	   << "\n" 
	   <<  std::setw(10)<< dish 
	   <<  std::setw(10)<< disc
	   << "\n";
	
	os << std::setw(10) << "i"
	   << std::setw(8) << "j"
	   << std::setw(8) << "k"
	   << std::setw(8) << "l"
	   << std::setw(5) << "type"
	   << std::setw(10) << "i"
	   << std::setw(8) << "j"
	   << std::setw(8) << "k"
	   << std::setw(8) << "l"
	   << std::setw(5) << "type"
	   << std::setw(5) << "n"
	   << std::setw(5) << "m"
	   << std::setw(8) << "A_r0"
	   << std::setw(8) << "A_w0"
	   << std::setw(8) << "B_r0"
	   << std::setw(8) << "B_w0"
	   << std::setw(4) << "rah"
	   << "\n";
      }
      
      for(unsigned int line_number=0; it != to; ++line_number, ++it){
	
	DEBUG(11, "\tnr " << line_number-2);
	
	int type1, type2;
	int n,m;
	std::vector<int> atom1, atom2;
	double A_r0, B_r0,A_w0, B_w0;
	int rah;
	
	_lineStream.clear();
	_lineStream.str(*it);
	
	for(unsigned int i = 0; i < io::In_Distanceres::MAX_ATOMS; i++) {
	  int atom;
	  _lineStream >> atom;
	  if (atom > 0) atom1.push_back(atom - 1);
          else if (atom <0){
            std::ostringstream msg;
            msg << "COM and COG type not possible for more than " 
                    << io::In_Distanceres::MAX_ATOMS << " atoms"  << std::endl;
            io::messages.add(msg.str(),
			 "In_Distanceres",
			 io::message::error);              
            nr_atoms=false;              
          
          }
	}      
	_lineStream >> type1;
	
	for(unsigned int i = 0; i < io::In_Distanceres::MAX_ATOMS; i++) {
	  int atom;
	  _lineStream >> atom;
	  if (atom > 0) atom2.push_back(atom - 1);
          else if (atom  <0){
            std::ostringstream msg;
            msg << "COM and COG type not possible for more than  " 
                    << io::In_Distanceres::MAX_ATOMS << " atoms"  << std::endl;
            io::messages.add(msg.str(),
			 "In_Distanceres",
			 io::message::error);              
            nr_atoms=false;              
          
          }
	}      
	
	_lineStream >> type2;
	_lineStream >> n >> m >> A_r0 >> A_w0 >> B_r0 >> B_w0 >> rah;
	
	if(_lineStream.fail()){
	  std::ostringstream msg;
	  msg << "bad line in PERTDISRESSPEC block: " << line_number << std::endl
	      << "          " << *it;
	  io::messages.add(msg.str(),
			   "In_Distanceres",
			   io::message::error);
	}
	if( nr_atoms){
	  util::virtual_type t1 = util::virtual_type(type1);
	  util::virtual_type t2 = util::virtual_type(type2);
	
	  util::Virtual_Atom v1(t1, atom1, dish, disc);
	  util::Virtual_Atom v2(t2, atom2, dish, disc);
	
	  topo.perturbed_distance_restraints().push_back
	    (topology::perturbed_distance_restraint_struct(v1,v2,n,m,A_r0,B_r0,A_w0,B_w0, rah));
	
	  if (!quiet){
	    for(unsigned int i = 0; i < io::In_Distanceres::MAX_ATOMS; i++) {
	      os << std::setw(i == 0 ? 10 : 8) << (i < atom1.size() ? atom1[i]+1 : 0);
	    }
	    os << std::setw(5) << type1;
	  
	    for(unsigned int i = 0; i < io::In_Distanceres::MAX_ATOMS; i++) {
	      os << std::setw(i == 0 ? 10 : 8) << (i < atom2.size() ? atom2[i]+1 : 0);
	    }
	    os << std::setw(5) << type2
	       << std::setw(5) << n
	       << std::setw(5) << m
	       << std::setw(8) << A_r0	    
	       << std::setw(8) << A_w0
	       << std::setw(8) << B_r0
	       << std::setw(8) << B_w0
	       << std::setw(8) << rah
	       << "\n";
	  }
        }
      }
      
    }//PERTDISRESPEC DISTANCERES
    
    //if (!quiet) os << "END\n";
    
  }

  { // DISTANCEFIELD RES
    DEBUG(10, "DFRESSPEC block");

    std::vector<std::string> buffer;
    std::string s;

    buffer = m_block["DFRESSPEC"];
    block_read.insert("DFRESSPEC");
    if(!buffer.size()){
      if(sim.param().distancefield.distancefield &&
	 !sim.param().perturbation.perturbation)
	// only need to warn if it was expected
	io::messages.add("no DFRESSPEC block in distance restraints file",
			 "in_distanceres", io::message::warning);
    }
    else if (buffer.size()<=2){
      io::messages.add("empty DFRESSPEC block in distance restraints file",
		       "in_distanceres", io::message::warning);
    }
    else{      
      
      DEBUG(10, "reading in DFRES data");
      
      if (!quiet) {
	switch (sim.param().distancefield.distancefield) {
	  case 0:
	    os << "\tDistancefield restraints OFF\n";
	    // how did you get here?
	    break;
	  case 1:
	    os << "\tDistancefield restraints ON\n";
	    break;
	  default:
	    os << "\tDistancefield restraints ERROR\n";
	}
      }
      
      _lineStream.clear();
      _lineStream.str(concatenate(buffer.begin() + 1, buffer.end() - 1, s));
      
      double dish,disc;
      _lineStream >> dish >> disc;
      
      if(_lineStream.fail()){
	std::ostringstream msg;
	msg << "bad line in DFRESSPEC block: failed to read in DISH and DISC"  << std::endl;
	io::messages.add(msg.str(),
			 "In_Distanceres",
			 io::message::error);
      }
      
      /*      if (!quiet){
	      os << std::setw(10) << "DISH"
	      << std::setw(10) << "DISC"
	      << "\n" 
	      <<  std::setw(10)<< dish
	      <<  std::setw(10)<< disc
	      << "\n";
	      }
      */
      int vtype_i, vtype_j, num, atom;
      std::vector<int> atomi, atomj;
      
      topo.disfield_restraints().on = true;
      
      _lineStream >> topo.disfield_restraints().proteinatoms; 
      topo.disfield_restraints().proteinatoms--;
      
      _lineStream >> topo.disfield_restraints().K
		  >> topo.disfield_restraints().r0;
      
      if (topo.disfield_restraints().proteinatoms < 0) {
	io::messages.add("DFRESSPEC block: PROTEINATOMS must be >= 0.",
			 "In_Distanceres", io::message::error);
      }
      if (topo.disfield_restraints().K < 0.0) {
	io::messages.add("DFRESSPEC block: K must be >= 0.0.",
			 "In_Distanceres", io::message::error);
      }
      if (topo.disfield_restraints().r0 < 0.0) {
	io::messages.add("DFRESSPEC block: r0 must be >= 0.0.",
			 "In_Distanceres", io::message::error);
      }
      
      _lineStream >> vtype_i
		  >> num;
      
      for(int i=0; i< num; i++){
	_lineStream >> atom;
	if (atom <  0 ) {
	  io::messages.add("DFRESSPEC block: ATOM_I must be >= 0.",
			   "In_Distanceres", io::message::error);
	}
	atomi.push_back(atom-1);
      }
      _lineStream >> vtype_j
		  >> num;
      
      for(int i=0; i< num; i++){
	_lineStream >> atom;
	if (atom <  0 ) {
	  io::messages.add("DFRESSPEC block: ATOM_J must be >= 0.",
			   "In_Distancres", io::message::error);
	}
	atomj.push_back(atom-1);
      }
      
      if(_lineStream.fail()){
	std::ostringstream msg;
	msg << "bad line in DFRESSPEC block: " << std::endl;
	io::messages.add(msg.str(),
			 "In_Distanceres",
			 io::message::error);
      }
      
      util::virtual_type t1 = util::virtual_type(vtype_i);
      util::virtual_type t2 = util::virtual_type(vtype_j);
      
      util::Virtual_Atom v1(t1, atomi, dish, disc);
      util::Virtual_Atom v2(t2, atomj, dish, disc);
      
      topo.disfield_restraints().v1 = v1;    
      topo.disfield_restraints().v2 = v2;     
      
    }
    
  } // DISTANCEFIELD
    
  { // PERTDFRESPEC DISTANCERES
    DEBUG(10, "PERTDFRESSPEC block");

    std::vector<std::string> buffer;
    std::string s;

    buffer = m_block["PERTDFRESSPEC"];
    block_read.insert("PERTDFRESSPEC");

    // check whether there is s.th. in the block
    if(!buffer.size()){
     // this is only bad if no DFRESSPEC block was specified either
      // you may be doing a perturbation but not a perturbed disfield
      if(sim.param().distancefield.distancefield && 
	 sim.param().perturbation.perturbation &&
	 !topo.disfield_restraints().on){
	io::messages.add("no (PERT)DFRESSPEC block in distance restraints file",
			 "in_distanceres", io::message::warning);
      }
    }
    else if (buffer.size()<=2){
      io::messages.add("empty PERTDFRESSPEC block in distance restraints file",
		       "in_distanceres", io::message::warning);
    }
    else{
      
      DEBUG(10, "reading in DISTANCERES (PERTDFRESSPEC) data");
      sim.param().perturbation.perturbed_par=true;
      
      if (!quiet){
	switch(sim.param().distancefield.distancefield * sim.param().perturbation.perturbation){
	  case 0:
	    os << "\tPerturbed Distancefield restraints OFF\n";
	    // how did you get here?
	    break;
	  case 1:
	    os << "\tPerturbed Distancefield restraints ON\n";
	    break;
	  default:
	    os << "\tPerturbed Distancefield restraints ERROR\n";
	}
      }
      
      _lineStream.clear();
      _lineStream.str(concatenate(buffer.begin() + 1, buffer.end() - 1, s));
      
      double dish,disc;
      _lineStream >> dish >> disc;
      
      if(_lineStream.fail()){
	std::ostringstream msg;
	msg << "bad line in PERTDFRESSPEC block: failed to read in DISH and DISC"  << std::endl;
	io::messages.add(msg.str(),
			 "In_Distanceres",
			 io::message::error);
      }
  
/*      if (!quiet){
	os << std::setw(10) << "DISH"
	   << std::setw(10) << "DISC"
	   << "\n" 
	   <<  std::setw(10)<< dish 
	   <<  std::setw(10)<< disc
	   << "\n";
      }
*/      
      int vtype_i, vtype_j, num, atom;
      std::vector<int> atomi, atomj;

      topo.perturbed_disfield_restraints().on = true;
      
      _lineStream >> topo.perturbed_disfield_restraints().proteinatoms;
      topo.perturbed_disfield_restraints().proteinatoms--;

      _lineStream >> topo.perturbed_disfield_restraints().A_r0 
                  >> topo.perturbed_disfield_restraints().K_A 
                  >> topo.perturbed_disfield_restraints().B_r0 
                  >> topo.perturbed_disfield_restraints().K_B 
                  >> topo.perturbed_disfield_restraints().n 
                  >> topo.perturbed_disfield_restraints().m;

      if (topo.perturbed_disfield_restraints().proteinatoms < 0) {
        io::messages.add("PERTDFRESSPEC block: PROTEINATOMS must be >= 0.",
                     "In_Distanceres", io::message::error);
      }
      if (topo.perturbed_disfield_restraints().K_A < 0.0 || 
          topo.perturbed_disfield_restraints().K_B < 0.0) {
        io::messages.add("PERTDFRESSPEC block: K must be >= 0.0.",
                     "In_Distanceres", io::message::error);
      }
      if (topo.perturbed_disfield_restraints().A_r0 < 0.0 || 
          topo.perturbed_disfield_restraints().B_r0 < 0.0) {
        io::messages.add("PERTDFRESSPEC block: r0 must be >= 0.0.",
                     "In_Distanceres", io::message::error);
      } 


      _lineStream >> vtype_i
                  >> num;

      for(int i=0; i< num; i++){
        _lineStream >> atom;
        if (atom <  0 ) {
        io::messages.add("PERTDFRESSPEC block: ATOM_I must be >= 0.",
                         "In_Distanceres", io::message::error);
        }
        atomi.push_back(atom-1);
      }
      _lineStream >> vtype_j
                  >> num;

      for(int i=0; i< num; i++){
        _lineStream >> atom;
        if (atom <  0 ) {
          io::messages.add("PERTDFRESSPEC block: ATOM_J must be >= 0.",
                           "In_Distancres", io::message::error);
        }
        atomj.push_back(atom-1);
      }
	
      if(_lineStream.fail()){
        std::ostringstream msg;
        msg << "bad line in PERTDFRESSPEC block: " << std::endl;
        io::messages.add(msg.str(), "In_Distanceres",
	 		io::message::error);
      }

      util::virtual_type t1 = util::virtual_type(vtype_i);
      util::virtual_type t2 = util::virtual_type(vtype_j);
	
      util::Virtual_Atom v1(t1, atomi, dish, disc);
      util::Virtual_Atom v2(t2, atomj, dish, disc);
	
      topo.perturbed_disfield_restraints().v1 = v1;
      topo.perturbed_disfield_restraints().v2 = v2;
	
    }
      
  }//PERTDFRESPEC DISTANCERES
    
  
  { // EDSDISTANCERES
    DEBUG(10, "MDISRESSPEC block");
    buffer = m_block["MDISRESSPEC"];
    block_read.insert("MDISRESSPEC");
    if (buffer.size()<=2){
      io::messages.add("no or empty MDISRESSPEC block in distance restraints file",
		       "in_distanceres", io::message::warning);
    }
    else{
      if (!sim.param().eds.eds){
        io::messages.add("MDISRESSPEC block given but EDS not turned on!",
                         "in_distanceres", io::message::error);
      }
      if(sim.param().distanceres.distanceres < 0){
        io::messages.add("eds perturbed distance restraints not compatible with time averaging!",
                         "in_distanceres", io::message::error);
      }
      std::vector<std::string>::const_iterator it = buffer.begin()+1,
	to = buffer.end()-1;
      
      double dish,disc;
      bool nr_atoms=true;
      
      sim.param().perturbation.perturbed_par=true;
      
      DEBUG(10, "reading in MDISRESSPEC data");
      
      if (!quiet){
	
	switch(sim.param().distanceres.distanceres){
	  case 0:
	    os << "\tEDS distance restraints OFF\n";
	    // how did you get here?
	    break;
	  case 1:
	    os << "\tEDS distance restraints ON\n";
	    
	    break;
	  case 2:
	    os << "\tEDS distance restraints ON\n"
	       << "\t\t(using force constant K*w0)\n";
	    break;
	  default:
	    os << "\tEDS distance restraints ERROR\n";
	}
      }
      
      
      _lineStream.clear();
      _lineStream.str(*it);
      
      _lineStream >> dish >> disc;
      
      if(_lineStream.fail()){
	std::ostringstream msg;
	msg << "bad line in MDISRESSPEC block: failed to read in DISH and  DISC"  << std::endl;
	io::messages.add(msg.str(),
			 "In_Distanceres",
			 io::message::error);
      }     
      ++it;
      
      const unsigned int numstates = sim.param().eds.numstates;
      
      if (!quiet){
	os << std::setw(10) << "DISH"
	   << std::setw(10) << "DISC"
	   << "\n" 
	   <<  std::setw(10)<< dish 
	   <<  std::setw(10)<< disc
	   << "\n";
	
	os << std::setw(10) << "i"
	   << std::setw(8) << "j"
	   << std::setw(8) << "k"
	   << std::setw(8) << "l"
	   << std::setw(5) << "type"
	   << std::setw(10) << "i"
	   << std::setw(8) << "j"
	   << std::setw(8) << "k"
           << std::setw(8) << "l"
           << std::setw(5) << "type";
        
        os << "\t" << std::setw(12) << "r0[";
        for(unsigned int i = 0; i < numstates; i++){
          os << std::setw(8) << i+1;
        }
        os << std::setw(12) << "] w0[";
        for(unsigned int i = 0; i < numstates; i++){
          os << std::setw(10) << i+1;
        }
        os << "] ";
	  
	os << std::setw(4) << "rah"
	   << "\n";
      }
    
      for(unsigned int line_number=2; it != to; ++line_number, ++it){
	
	DEBUG(11, "\tnr " << line_number - 2);
	
	int type1, type2;
	std::vector<int> atom1, atom2;
	std::vector<double> r0(numstates),w0(numstates);
	int rah;
	
	_lineStream.clear();
	_lineStream.str(*it);
	
	for(unsigned int i = 0; i < io::In_Distanceres::MAX_ATOMS; i++) {
	  int atom;
	  _lineStream >> atom;
	  // -1 because we directly convert to array indices
	  if (atom > 0) atom1.push_back(atom - 1); 
          else if (atom <0){
            std::ostringstream msg;
            msg << "COM and COG type not possible for more than " 
                    << io::In_Distanceres::MAX_ATOMS << " atoms"  << std::endl;
            io::messages.add(msg.str(),
			 "In_Distanceres",
			 io::message::error);              
            nr_atoms=false;
          
          }
	}      
	_lineStream >> type1;
	
	for(unsigned int i = 0; i < io::In_Distanceres::MAX_ATOMS; i++) {
	  int atom;
	  _lineStream >> atom;
	  if (atom > 0) atom2.push_back(atom - 1);
          else if (atom <0){
            std::ostringstream msg;
            msg << "COM and COG type not possible for more than " 
                    << io::In_Distanceres::MAX_ATOMS << " atoms"  << std::endl;
            io::messages.add(msg.str(),
			 "In_Distanceres",
			 io::message::error);              
            nr_atoms=false;              
          
          }          
	}      
	_lineStream >> type2;
        
        for(unsigned int i = 0; i < numstates; i++){
          _lineStream >> r0[i];
        }
        for(unsigned int i = 0; i < numstates; i++){
          _lineStream >> w0[i];
        }
        _lineStream >> rah;
	
	
	if(_lineStream.fail()){
	  std::ostringstream msg;
	  msg << "bad line in MDISRESSPEC block: " << line_number << std::endl
	      << "          " << *it;
	  io::messages.add(msg.str(),
			   "In_Distanceres",
			   io::message::error);
	}
	
	// g++ 3.2 fix
        if( nr_atoms){
	  util::virtual_type t1 = util::virtual_type(type1);
	  util::virtual_type t2 = util::virtual_type(type2);
	
	  util::Virtual_Atom v1(t1, atom1, dish, disc);
	  util::Virtual_Atom v2(t2, atom2, dish, disc);
	
	  topo.eds_distance_restraints().push_back
	    (topology::eds_distance_restraint_struct(v1,v2,r0,w0,rah));
	
	  if (!quiet){
	    for(unsigned int i = 0; i < io::In_Distanceres::MAX_ATOMS; i++) {
	      // the first element has the width 10, if i is bigger then the number of atoms
	      // specified, just print 0.
	      os << std::setw(i == 0 ? 10 : 8) << (i < atom1.size() ? atom1[i]+1 : 0);
	    }
	    os << std::setw(5) << type1;
	
            for(unsigned int i = 0; i < io::In_Distanceres::MAX_ATOMS; i++) {
              os << std::setw(i == 0 ? 10 : 8) << (i < atom2.size() ? atom2[i]+1 : 0);
            }
            os << std::setw(5) << type2;
            os << "\t" << std::setw(12) << "   ";
            for(unsigned int i = 0; i < numstates; i++){
              os << std::setw(8) << r0[i];
            }
            os << std::setw(12) << "     ";
            for(unsigned int i = 0; i < numstates; i++){
              os << std::setw(10) << w0[i];
            }
            os << "  ";
          
            os << std::setw(4) << rah << "\n";
          
	  }
        }
      }
      
    }
    if (!quiet) os << "END\n";
  } // EDSDISTANCERES

    
  for(std::map<std::string, std::vector<std::string> >::const_iterator
	it = m_block.begin(),
	to = m_block.end();
      it != to;
      ++it){
    
    if (block_read.count(it->first) == 0 && it->second.size()){
      io::messages.add("block " + it->first + " not supported!",
		       "In_Distanceres",
		       io::message::warning);
    }
  }
  
}
