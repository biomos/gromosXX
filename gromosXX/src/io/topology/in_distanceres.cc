/**
 * @file in_distanceres.cc
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

#include "in_distanceres.h"

#undef MODULE
#undef SUBMODULE
#define MODULE io
#define SUBMODULE topology

static std::set<std::string> block_read;

void 
io::In_Distanceres::read(topology::Topology& topo,
		      simulation::Simulation & sim,
		      std::ostream & os){
  
  DEBUG(7, "reading in a distance restraints file");

  if (!quiet)
    os << "DISTRANCE RESTRAINTS\n";
  
  std::vector<std::string> buffer;
  std::vector<std::string>::const_iterator it;

  { // DISTANCERES
    DEBUG(10, "DISRESSPEC block");
    buffer = m_block["DISRESSPEC"];
    block_read.insert("DISRESSPEC");
    if (buffer.size()<=2){
      io::messages.add("no or empty DISRESSPEC block in distance restraints file",
		       "in_distanceres", io::message::warning);
    }
    else{      
      std::vector<std::string>::const_iterator it = buffer.begin()+1,
	to = buffer.end()-1;
      
      double dish,disc,disn;
      
      DEBUG(10, "reading in DISTANCERES data");
      
      if (!quiet){
	
	switch(sim.param().distanceres.distanceres){
	  case 0:
	    os << "\tDistance restraints OFF\n";
	    // how did you get here?
	    break;
	  case 1:
	    os << "\tDistance restraints ON\n";
	    
	    break;
	  case 2:
	    os << "\tDistance restraints ON\n"
	       << "\t\t(using force constant K*w0)\n";
	    break;
	  default:
	    os << "\tDistance restraints ERROR\n";
	}
      }
      
      
      _lineStream.clear();
      _lineStream.str(*it);
      
      _lineStream >> dish >> disc >> disn;
      
      if(_lineStream.fail()){
	std::ostringstream msg;
	msg << "bad line in DISRESSPEC block: failed to read in DISH , DISC and DISN"  << std::endl;
	io::messages.add(msg.str(),
			 "In_Distanceres",
			 io::message::error);
      }
      
      ++it;
      
      if (!quiet){
	os << std::setw(10) << "DISH"
	   << std::setw(10) << "DISC"
           << std::setw(10) << "DISN"
	   << "\n" 
	   <<  std::setw(10)<< dish 
	   <<  std::setw(10)<< disc
           <<  std::setw(10)<< disn
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
	}      
	_lineStream >> type1;
	
	for(unsigned int i = 0; i < io::In_Distanceres::MAX_ATOMS; i++) {
	  int atom;
	  _lineStream >> atom;
	  if (atom > 0) atom2.push_back(atom - 1);
	}      
	_lineStream >> type2;
	
	_lineStream >> r0>> w0 >> rah;
	
	
	if(_lineStream.fail()){
	  std::ostringstream msg;
	  msg << "bad line in DISRESSPEC block: " << line_number << std::endl
	      << "          " << *it;
	  io::messages.add(msg.str(),
			   "In_Distanceres",
			   io::message::error);
	}
	
	// g++ 3.2 fix
	util::virtual_type t1 = util::virtual_type(type1);
	util::virtual_type t2 = util::virtual_type(type2);
	
	util::Virtual_Atom v1(t1, atom1, dish, disc, disn);
	util::Virtual_Atom v2(t2, atom2, dish, disc, disn);
	
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
      std::vector<std::string>::const_iterator it = buffer.begin()+1,
	to = buffer.end()-1;
      
      double dish, disc, disn;
      
      DEBUG(10, "reading in DISTANCERES (PERTDISRESSPEC data");
      
      if (!quiet){
	switch(sim.param().distanceres.distanceres){
	  case 0:
	    os << "\tPerturbed Distance restraints OFF\n";
	    // how did you get here?
	    break;
	  case 1:
	    os << "\tPerturbed Distance restraints ON\n";
	    
	    break;
	  case 2:
	    os << "\tPerturbed Distance restraints ON\n"
	       << "\t\t(using force constant K*w0)\n";
	    break;
	  default:
	    os << "\tPerturbed Distance restraints ERROR\n";
	}
      }
      
      _lineStream.clear();
      _lineStream.str(*it);
      
      _lineStream >> dish >> disc >> disn;
      
      if(_lineStream.fail()){
	std::ostringstream msg;
	msg << "bad line in PERTDISRESSPEC block: failed to read in DISH, DISC and DISN"  << std::endl;
	io::messages.add(msg.str(),
			 "In_Distanceres",
			 io::message::error);
      }
      
      ++it;
      if (!quiet){
	os << std::setw(10) << "DISH"
	   << std::setw(10) << "DISC"
           << std::setw(10) << "DISN"
	   << "\n" 
	   <<  std::setw(10)<< dish 
	   <<  std::setw(10)<< disc
           <<  std::setw(10)<< disn
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
	}      
	_lineStream >> type1;
	
	for(unsigned int i = 0; i < io::In_Distanceres::MAX_ATOMS; i++) {
	  int atom;
	  _lineStream >> atom;
	  if (atom > 0) atom2.push_back(atom - 1);
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
	
	util::virtual_type t1 = util::virtual_type(type1);
	util::virtual_type t2 = util::virtual_type(type2);
	
	util::Virtual_Atom v1(t1, atom1, dish, disc, disn);
	util::Virtual_Atom v2(t2, atom2, dish, disc, disn);
	
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
      
    }//PERTDISRESPEC DISTANCERES
    
    //if (!quiet) os << "END\n";
    
  }
  
  { // EDSDISTANCERES
    DEBUG(10, "EDSDISRESSPEC block");
    buffer = m_block["EDSDISRESSPEC"];
    block_read.insert("EDSDISRESSPEC");
    if (buffer.size()<=2){
      io::messages.add("no or empty EDSDISRESSPEC block in distance restraints file",
		       "in_distanceres", io::message::warning);
    }
    else{
      if (!sim.param().eds.eds){
        io::messages.add("EDSDISRESSPEC block given but EDS not turned on!",
                         "in_distanceres", io::message::error);
      }
      if(sim.param().distanceres.distanceres < 0){
        io::messages.add("eds perturbed distance restraints not compatible with time averaging!",
                         "in_distanceres", io::message::error);
      }
      std::vector<std::string>::const_iterator it = buffer.begin()+1,
	to = buffer.end()-1;
      
      double dish,disc,disn;
            
      DEBUG(10, "reading in EDSDISRESSPEC data");
      
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
      
      _lineStream >> dish >> disc >> disn;
      
      if(_lineStream.fail()){
	std::ostringstream msg;
	msg << "bad line in EDSDISRESSPEC block: failed to read in DISH , DISC and DISN"  << std::endl;
	io::messages.add(msg.str(),
			 "In_Distanceres",
			 io::message::error);
      }
      
      ++it;
      
      const unsigned int numstates = sim.param().eds.numstates;
      
      if (!quiet){
	os << std::setw(10) << "DISH"
	   << std::setw(10) << "DISC"
           << std::setw(10) << "DISN"
	   << "\n" 
	   <<  std::setw(10)<< dish 
	   <<  std::setw(10)<< disc
           <<  std::setw(10)<< disn
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
	}      
	_lineStream >> type1;
	
	for(unsigned int i = 0; i < io::In_Distanceres::MAX_ATOMS; i++) {
	  int atom;
	  _lineStream >> atom;
	  if (atom > 0) atom2.push_back(atom - 1);
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
	  msg << "bad line in EDSDISRESSPEC block: " << line_number << std::endl
	      << "          " << *it;
	  io::messages.add(msg.str(),
			   "In_Distanceres",
			   io::message::error);
	}
	
	// g++ 3.2 fix
	util::virtual_type t1 = util::virtual_type(type1);
	util::virtual_type t2 = util::virtual_type(type2);
	
	util::Virtual_Atom v1(t1, atom1, dish, disc, disn);
	util::Virtual_Atom v2(t2, atom2, dish, disc, disn);
	
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
