/**
 * @file in_dihrest.cc
 * implements methods of In_Dihrest
 */


#include <stdheader.h>

#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <simulation/simulation.h>
#include <interaction/interaction_types.h>
#include <configuration/configuration.h>

#include <io/instream.h>
#include <io/blockinput.h>

#include "in_dihrest.h"

#undef MODULE
#undef SUBMODULE
#define MODULE io
#define SUBMODULE topology

static std::set<std::string> block_read;

void io::In_Dihrest::read(topology::Topology& topo,
			  configuration::Configuration & conf,
			  simulation::Simulation & sim){
  
  DEBUG(7, "reading in a dihedral restraints file");
  
  if (!quiet)
    std::cout << "DIHEDRAL RESTRAINTS\n";
  
  std::vector<std::string> buffer;
  std::vector<std::string>::const_iterator it;
  
  { // DIHRESSPEC
    DEBUG(10, "DIHRESSPEC block");
    buffer = m_block["DIHRESSPEC"];
    
    if (!buffer.size()){
      io::messages.add("no DIHRESSPEC block in dihedral restraints file",
		       "in_dihrest", io::message::error);
      return;
    }

    block_read.insert("DIHRESSPEC");

    std::vector<std::string>::const_iterator it = buffer.begin()+1,
      to = buffer.end()-1;

    int i, j, k, l;
    double delta, phi, w0;

    DEBUG(10, "reading in DIHREST data");

    if (!quiet){
      switch(sim.param().dihrest.dihrest){
	case 0:
	  std::cout << "\tDihedral restraints OFF\n";
	  // how did you get here?
	  break;
	case 1:
	  std::cout << "\tDihedral restraints ON\n"
		    << "\t\t(uniform force constant K)\n";
	  break;
	case 2:
	  std::cout << "\tDihedral restraints ON\n"
		    << "\t\t(force constant K*w0)\n";
	  break;
	default:
	  std::cout << "\tDihedral restraints: ERROR\n";
	  io::messages.add("wrong value for method in dihedral restraints block",
			   "in_dihedral", io::message::error);
	  return;
      }
    }
    
    if (!quiet){

      std::cout << std::setw(10) << "i"
		<< std::setw(8) << "j"
		<< std::setw(8) << "k"
		<< std::setw(8) << "l"
		<< std::setw(8) << "delta"
		<< std::setw(8) << "phi"
		<< std::setw(8) << "w0"
		<< "\n";

      std::cout.precision(2);
      std::cout.setf(std::ios::fixed, std::ios::floatfield);

    }
    
    for(int c=0; it != to; ++c, ++it){
      
      DEBUG(11, "\tnr " << c);
      
      _lineStream.clear();
      _lineStream.str(*it);

      _lineStream >> i >> k >> j >> l >> delta >> phi >> w0;
    
      if(_lineStream.fail()){
	io::messages.add("bad line in DIHREST block",
			 "In_Dihrest", io::message::error);
      }
      
      topo.dihedral_restraints().push_back
	(topology::dihedral_restraint_struct(i-1, j-1, k-1, l-1,
					     delta * 2 * math::Pi / 360, phi * 2 * math::Pi / 360, w0));
      
      if (!quiet){
	std::cout << std::setw(10) << i
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

    std::vector<std::string>::const_iterator it = buffer.begin()+1,
      to = buffer.end()-1;

    int i, j, k, l, m, n;
    double delta, A_phi, A_w0, B_phi, B_w0;
    
    DEBUG(10, "reading in perturbed DIHREST (PERTDIHRESSPEC data");

    if (!quiet){
      switch(sim.param().distrest.distrest){
	case 0:
	  std::cout << "\tPerturbed Dihedral restraints OFF\n";
	  // how did you get here?
	  break;
	case 1:
	  std::cout << "\tPerturbed Dihedral restraints ON\n"
		    << "\t\t(using uniform force constant K\n";
	  break;
	case 2:
	  std::cout << "\tPerturbed Dihedral restraints ON\n"
		    << "\t\t(using force constant K*w0)\n";
	  break;
	default:
	  std::cout << "\tPerturbed Dihedral restraints ERROR\n";
	  io::messages.add("wrong method for dihedral restraints",
			   "in_dihrest", io::message::error);
	  return;
      }
    }

    if (!quiet){
      std::cout << std::setw(10) << "i"
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

      std::cout.precision(2);
      std::cout.setf(std::ios::fixed, std::ios::floatfield);
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
	std::cout << std::setw(10) << i
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
    
    if (!quiet) std::cout << "END\n";
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