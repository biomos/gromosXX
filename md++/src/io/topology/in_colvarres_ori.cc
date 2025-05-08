/**
 * @file in_colvarres.cc
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

#include "in_colvarres.h"

#undef MODULE
#undef SUBMODULE
#define MODULE io
#define SUBMODULE topology

static std::set<std::string> block_read;

/**
 * @section contactnumresspec CONTACTNUMRESSPEC block
 * The CONTACTNUMRESSPEC block is read from the colvar restraints specification
 * file.
 *
 * See @ref util::virtual_type for valid virtual atom types.
 *
 * @verbatim
CONTACTNUMRESSPEC
# 
# CONT0: target value of collective variable
# W0: scaling factor for the force constant specified in the COLVARRES block
# RCUT: cutoff (nm)
# N: numerator exponent
# M: denominator exponent
# type: virtual atom type, if not =0, it will try to create one single virtual atom from all given atoms
# G1: >0: number of atoms in group 1, 
#     <0: number of pairs of atoms, each specifying a range of atoms that will be added to group 1
# G2: >0: number of atoms in group 2, 
#     <0: number of pairs of atoms, each specifying a range of atoms that will be added to group 2
# MASS: use only atoms heavier than this
# CONT0  W0   RCUT    N      M       MASS
    90  0.0001   0.7       2    24   1.008
# type  G1   ATOMN(1..G1)
      0 -1      1 15
# type  G2   ATOMN(1..G2)
      0 -1   16  30
END
@endverbatim
 * @verbatim
PERTCONTACTNUMRESSPEC
# 
# CONT0_A, _B: target value in state A and B
# W0_A, _B: scaling factor for the force constant specified in the COLVARRES block for state A and B
# RCUT: cutoff (nm)
# N: numerator exponent
# M: denominator exponent
# type: virtual atom type, if not =0, it will try to create one single virtual atom from all given atoms
# G1: >0: number of atoms in group 1, 
#     <0: number of pairs of atoms, each specifying a range of atoms that will be added to group 1
# G2: >0: number of atoms in group 2, 
#     <0: number of pairs of atoms, each specifying a range of atoms that will be added to group 2
# MASS: use only atoms heavier than this
# CONT0_A  W0_A CONT0_B W0_B   RCUT    N      M       MASS
    90  0.0001     40  0.0001  0.7       2    24   1.008
# type  G1   ATOMN(1..G1)
      0 -1      1 15
# type  G2   ATOMN(1..G2)
      0 -1   16  30
END
@endverbatim
 */
void 
io::In_Colvarres::read(topology::Topology& topo,
		      simulation::Simulation & sim,
		      std::ostream & os){
  
  DEBUG(7, "reading in a colvar restraints file");

  if (!quiet)
    os << "COLVAR RESTRAINTS\n";
  
  std::vector<std::string> buffer;

  { // CONTACTNUMRES
    DEBUG(10, "CONTACTNUMRESSPEC block");
    buffer = m_block["CONTACTNUMRESSPEC"];
    block_read.insert("CONTACTNUMRESSPEC");
    if (buffer.size()<=2){
      io::messages.add("no or empty CONTACTNUMRESSPEC block in colvar restraints file",
		       "in_colvarres", io::message::warning);
    }
    else {            
      double rcut, cont0, w0, masscut;
      int nn, mm;          
      
      DEBUG(10, "reading in CONTACTNUMRES data");
      
      if (!quiet) {
        switch (sim.param().colvarres.colvarres) {
          case 0:
            os << "\tColvar restraints OFF\n";
            // how did you get here?
            break;
          case 1:
            os << "\tColvar restraints ON\n";
            break;
          default:
            os << "\tColvar restraints ERROR\n";
        }
      }
      
      _lineStream.clear();
      std::string s; 
      _lineStream.str(concatenate(buffer.begin() + 1, buffer.end() - 1, s));
      
             
      _lineStream >> cont0>> w0 >> rcut >> nn >> mm >> masscut;
      
      if (_lineStream.fail()) {
	    std::ostringstream msg;
	    msg << "bad line in CONTACTNUMRESSPEC block: failed to read in CONT0, N and M"  << std::endl;
	    io::messages.add(msg.str(),
			 "In_Colvarres",
			 io::message::error);
      }
            
      if (!quiet){
	    os << std::setw(10) << "CONT0"
	       << std::setw(10) << "W0"
	       << std::setw(10) << "RCUT"
	       << std::setw(10) << "N"
	       << std::setw(10) << "M"
	       << std::setw(10) << "MASS"
	       << "\n" 
	       <<  std::setw(10)<< cont0
	       <<  std::setw(10)<< w0
	       <<  std::setw(10)<< rcut
	       <<  std::setw(10)<< nn
	       <<  std::setw(10)<< mm
	       <<  std::setw(10)<< masscut
	       << "\n";
      }
      	
	  // read atom numbers for two groups of atoms
	  std::vector<std::vector<int > > atomgroup(2);
	  std::vector< int> vatypes(2);
	  bool notinsolute=false;
      for (unsigned int ag=0; ag<=1; ag++) {	
	    int natoms;
        _lineStream >> vatypes[ag] >> natoms;
        if (!quiet) {
      	  os << std::setw(8) << " type G"<< ag+1 <<" ATOMN(1..G"<< ag+1<<")\n "
          << std::setw(4) << vatypes[ag] << std::setw(8) << natoms;
        }
	
	    for (unsigned int i = 0; i < abs(natoms); ++i) {  
  	      // single atoms were specified
  	      if (natoms>0){
  	        unsigned int x;
            _lineStream >> x;
	        if(_lineStream.fail()) break;
            if (x-1 <= topo.num_solute_atoms()) {
              if (topo.mass(x-1) > masscut) {
	            // -1 because we directly convert to array indices
                DEBUG(10, "\tsingle atoms specified | adding atom to contact group "<< ag+1 << " " << x - 1);
                atomgroup[ag].push_back(x - 1);
                if (!quiet)	os << "  " << x; 
              } 
            } else {
	          notinsolute=true;        
            }          
          } 
  	      // ranges of atoms were specified
          else if (natoms<0) {
  	        unsigned int start, end;
            _lineStream >> start >> end;
            if (end-1 <= topo.num_solute_atoms()) {
              if (start > end) {
                 io::messages.add("CONTACTNUMRESSPEC block: atom ranges not in order", 
                                       "In_Colvarres", io::message::error);              
              } else {
                for (unsigned int x=start; x<=end; x++) {
                  if (topo.mass(x-1) > masscut) {
	                // -1 because we directly convert to array indices
                    DEBUG(10, "\t ranges specified | adding atom to contact group "<< ag+1 << " " << x - 1);
                    atomgroup[ag].push_back(x - 1);  
                  }        
                } 
                if (!quiet)	os << "  " << start << "  " << end;
              }  
            } else {
	          notinsolute=true;        
            }          
          } 
          else {
            io::messages.add("CONTACTNUMRESSPEC block: atomnum 0 does not make sense", 
                                       "In_Colvarres", io::message::error);           
          }
        } // atoms of one group read
        if (!quiet)  os << "\n" ;
      
        // don't allow non-solvent atoms
        if (notinsolute) {
          std::ostringstream msg;
          msg << "CONTACTNUMRESSPEC block: not all specified atoms are in the solute"  << std::endl;
          io::messages.add(msg.str(), "In_Colvarres", io::message::error);
        }
      } // atoms of both atom groups read      	
	
	  if(_lineStream.fail()){
	    std::ostringstream msg;
	    msg << "bad line in CONTACTNUMRESSPEC block \n";
	    io::messages.add(msg.str(),
			   "In_Colvarres",
			   io::message::error);
	  }
	
	  // create the virtual atoms
	  std::vector<std::vector<util::Virtual_Atom > >atoms(2);
	  for (int ag=0; ag<2; ag++) {
        if(atomgroup[ag].size() != 0){
          // if real atom type, one va per atom
          if (vatypes[ag] == 0) {
            for (int i=0; i<atomgroup[ag].size(); i++) {        
              std::vector<int> atomvec(1,atomgroup[ag][i]);
              util::Virtual_Atom v(util::virtual_type(vatypes[ag]), atomvec, 0.1, 0.153);
              atoms[ag].push_back(v);
            }
          }
          // else: create one virtual atom from all the atoms of one group
          else {    
            std::vector<int> atomvec(atomgroup[ag].size());
            for (int i=0; i<atomgroup[ag].size(); i++) { 
              atomvec.push_back(atomgroup[ag][i]);
            }
            util::Virtual_Atom v(util::virtual_type(vatypes[ag]), atomvec, 0.1, 0.153);
            atoms[ag].push_back(v);
          }
        } else {
          io::messages.add("CONTACTNUMRESSPEC block: no atoms in atomgroup", 
                                       "In_Colvarres", io::message::error);          
        }
	  }  
      
	  topo.contactnum_restraint().push_back
	    (topology::contactnum_restraint_struct
	      (true, atoms[0], atoms[1], cont0, w0, rcut, nn, mm, masscut));
   
      if (!quiet) os << "END\n";
  
      int check_empty;
      _lineStream >> check_empty;
      if (!_lineStream.fail()){
	    std::ostringstream msg;
	    msg << "CONTACTNUMRESSPEC block: more arguments given than required, they are being ignored\n";
	    io::messages.add(msg.str(), "In_Colvarres", io::message::warning);
	  }
	} // end if buffer size>2
  } // CONTACTNUMRES
  
  
  { // PERTURBED CONTACTNUMRES
    DEBUG(10, "PERTCONTACTNUMRESSPEC block");
    buffer = m_block["PERTCONTACTNUMRESSPEC"];
    block_read.insert("PERTCONTACTNUMRESSPEC");
    if (buffer.size()<=2){
      io::messages.add("no or empty PERTCONTACTNUMRESSPEC block in colvar restraints file",
		       "in_colvarres", io::message::warning);
    }
    else {            
      double rcut, cont0A, cont0B, w0A, w0B, masscut;
      int nn, mm;          
      
      DEBUG(10, "reading in PERTCONTACTNUMRES data");
      
      if (!quiet) {
        switch (sim.param().colvarres.colvarres) {
          case 0:
            os << "\tColvar restraints OFF\n";
            // how did you get here?
            break;
          case 1:
            os << "\tColvar restraints ON\n";
            break;
          default:
            os << "\tColvar restraints ERROR\n";
        }
      }
      
      _lineStream.clear();
      std::string s; 
      _lineStream.str(concatenate(buffer.begin() + 1, buffer.end() - 1, s));
      
             
      _lineStream >> cont0A>> w0A >> cont0B >> w0B >> rcut >> nn >> mm >> masscut;
      
      if (_lineStream.fail()) {
	    std::ostringstream msg;
	    msg << "bad line in PERTCONTACTNUMRESSPEC block: failed to read in CONT0, N and M"  << std::endl;
	    io::messages.add(msg.str(),
			 "In_Colvarres",
			 io::message::error);
      }
            
      if (!quiet){
	    os << std::setw(10) << "CONT0A"
	       << std::setw(10) << "W0A"
	       << std::setw(10) << "CONT0B"
	       << std::setw(10) << "W0B"
	       << std::setw(10) << "RCUT"
	       << std::setw(10) << "N"
	       << std::setw(10) << "M"
	       << std::setw(10) << "MASS"
	       << "\n" 
	       <<  std::setw(10)<< cont0A
	       <<  std::setw(10)<< w0A
	       <<  std::setw(10)<< cont0B
	       <<  std::setw(10)<< w0B
	       <<  std::setw(10)<< rcut
	       <<  std::setw(10)<< nn
	       <<  std::setw(10)<< mm
	       <<  std::setw(10)<< masscut
	       << "\n";
      }
      	
	  // read atom numbers for two groups of atoms
	  std::vector<std::vector<int > > atomgroup(2);
	  std::vector< int> vatypes(2);
	  bool notinsolute=false;
      for (unsigned int ag=0; ag<=1; ag++) {	
	    int natoms;
        _lineStream >> vatypes[ag] >> natoms;
        if (!quiet) {
      	  os << std::setw(8) << " type G"<< ag+1 <<" ATOMN(1..G"<< ag+1<<")\n "
          << std::setw(4) << vatypes[ag] << std::setw(8) << natoms;
        }
	
	    for (unsigned int i = 0; i < abs(natoms); ++i) {  
  	      // single atoms were specified
  	      if (natoms>0){
  	        unsigned int x;
            _lineStream >> x;
	        if(_lineStream.fail()) break;
            if (x-1 <= topo.num_solute_atoms()) {
              if (topo.mass(x-1) > masscut) {
	            // -1 because we directly convert to array indices
                DEBUG(10, "\tadding atom to contact group "<< ag+1 << " " << x - 1);
                atomgroup[ag].push_back(x - 1);
                if (!quiet)	os << "  " << x; 
              } 
            } else {
	          notinsolute=true;        
            }          
          } 
  	      // ranges of atoms were specified
          else if (natoms<0) {
  	        unsigned int start, end;
            _lineStream >> start >> end;
            if (end-1 <= topo.num_solute_atoms()) {
              if (start > end) {
                 io::messages.add("PERTCONTACTNUMRESSPEC block: atom ranges not in order", 
                                       "In_Colvarres", io::message::error);              
              } else {
                for (unsigned int x=start; x<=end; x++) {
                  if (topo.mass(x-1) > masscut) {
	                // -1 because we directly convert to array indices
                    DEBUG(10, "\tadding atom to contact group "<< ag+1 << " " << x - 1);
                    atomgroup[ag].push_back(x - 1);  
                  }        
                } 
                if (!quiet)	os << "  " << start << "  " << end;
              }  
            } else {
	          notinsolute=true;        
            }          
          } 
          else {
            io::messages.add("PERTCONTACTNUMRESSPEC block: atomnum 0 does not make sense", 
                                       "In_Colvarres", io::message::error);           
          }
        } // atoms of one group read
        if (!quiet)  os << "\n" ;
      
        // don't allow non-solvent atoms
        if (notinsolute) {
          std::ostringstream msg;
          msg << "PERTCONTACTNUMRESSPEC block: not all specified atoms are in the solute"  << std::endl;
          io::messages.add(msg.str(), "In_Colvarres", io::message::error);
        }
      } // atoms of both atom groups read      	
	
	  if(_lineStream.fail()){
	    std::ostringstream msg;
	    msg << "bad line in PERTCONTACTNUMRESSPEC block \n";
	    io::messages.add(msg.str(),
			   "In_Colvarres",
			   io::message::error);
	  }
	
	  // create the virtual atoms
	  std::vector<std::vector<util::Virtual_Atom > >atoms(2);
	  for (int ag=0; ag<2; ag++) {
        if(atomgroup[ag].size() != 0){
          // if real atom type, one va per atom
          if (vatypes[ag] == 0) {
            for (int i=0; i<atomgroup[ag].size(); i++) {        
              std::vector<int> atomvec(1,atomgroup[ag][i]);
              util::Virtual_Atom v(util::virtual_type(vatypes[ag]), atomvec, 0.1, 0.153);
              atoms[ag].push_back(v);
            }
          }
          // else: create one virtual atom from all the atoms of one group
          else {    
            std::vector<int> atomvec(atomgroup[ag].size());
            for (int i=0; i<atomgroup[ag].size(); i++) { 
              atomvec.push_back(atomgroup[ag][i]);
            }
            util::Virtual_Atom v(util::virtual_type(vatypes[ag]), atomvec, 0.1, 0.153);
            atoms[ag].push_back(v);
          }
        } else {
          io::messages.add("PERTCONTACTNUMRESSPEC block: no atoms in atomgroup", 
                                       "In_Colvarres", io::message::error);          
        }
	  }  
      
	  topo.perturbed_contactnum_restraint().push_back
	    (topology::perturbed_contactnum_restraint_struct
	      (true, atoms[0], atoms[1], cont0A, w0A, cont0B, w0B, rcut, nn, mm, masscut));
   
      if (!quiet) os << "END\n";
  
      int check_empty;
      _lineStream >> check_empty;
      if (!_lineStream.fail()){
	    std::ostringstream msg;
	    msg << "PERTCONTACTNUMRESSPEC block: more arguments given than required, they are being ignored\n";
	    io::messages.add(msg.str(), "In_Colvarres", io::message::warning);
	  }
	} // end if buffer size>2
  } // PERTURBED CONTACTNUMRES
  
}
