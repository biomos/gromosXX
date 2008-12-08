/**
 * @file in_posres.cc
 * implements methods of In_Posresspec and In_Posres.
 */


#include <stdheader.h>

#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <simulation/simulation.h>
#include <interaction/interaction_types.h>
#include <configuration/configuration.h>

#include <io/instream.h>
#include <io/blockinput.h>
#include <io/configuration/in_configuration.h>
#include <io/configuration/out_configuration.h>

#include "in_posres.h"

#undef MODULE
#undef SUBMODULE
#define MODULE io
#define SUBMODULE topology

/**
 * @section posresspec POSRESSPEC block
 * The POSRESSPEC specifies the atoms of which the position is to be restrained
 * or constrained.
 *
 * The block is read from the position restraints specification file 
 * (\@posresspec).
 *
 * @verbatim
POSRESSPEC
# RESIDUE   ATOM
    1 HEXA  CH31       1
    1 HEXA  CH22       2
    1 HEXA  CH23       3
    1 HEXA  CH24       4
    2 HEXA  CH23       9
    2 HEXA  CH24      10
    2 HEXA  CH25      11
    2 HEXA  CH36      12
END
@endverbatim
 *
 * @section refposition REFPOSITION block
 * The REFPOSITION block has the same format as the POSITION block in the
 * configuration file and can be read either from the position restraints
 * file (\@posres) (NTPORB = 1) or the configuration file (NTPORB = 0, default).
 *
 * @verbatim
REFPOSITION
# RESIDUE   ATOM            POSITION
    1 HEXA  CH31       1    3.003183500    2.472059119    3.983607307
    1 HEXA  CH22       2    2.944583705    2.441994145    4.121706839
    1 HEXA  CH23       3    3.050801324    2.493790966    4.218887665
    1 HEXA  CH24       4    3.012945845    2.500345433    4.366985611
    2 HEXA  CH23       9    2.257723206    2.570648803    4.179681483
    2 HEXA  CH24      10    2.111496670    2.612272713    4.162523150
    2 HEXA  CH25      11    2.098679480    2.760542973    4.198030366
    2 HEXA  CH36      12    1.953046837    2.805540373    4.211273000
END
@endverbatim
 *
 * @section bfactor BFACTOR block
 * The BFACTOR block is optional. If not specified the atomic B-factors @f$b@f$
 * are set to @f$1@f$. The force constant is given as
 * @f[k = k_0 / b\mathrm{,}@f]
 * where @f$k_0@f$ is specified in the input parameter file by CPOR.
 *
 * The block is read from the position restraints file (\@posres).
 *
 * @verbatim
BFACTOR
# RESIDUE   ATOM           BFACTOR
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
io::In_Posresspec::read(topology::Topology& topo,
		    simulation::Simulation & sim,
		    std::ostream & os){
  
  DEBUG(7, "reading in a position restraints file");

  if (!quiet)
    os << "POSITION RESTRAINTS\n";
  
  std::vector<std::string> buffer;
  std::vector<std::string>::const_iterator it;

  { // POSRESSPEC

    buffer = m_block["POSRESSPEC"];
    DEBUG(10, "POSRESSPEC block : " << buffer.size());
    
    if (!buffer.size()){
      io::messages.add("no POSRESSPEC block in position restraints file",
		       "in_posres", io::message::error);
      return;
    }

    std::vector<std::string>::const_iterator it = buffer.begin()+1,
      to = buffer.end()-1;

    DEBUG(10, "reading in POSRESSPEC data");

    unsigned int i, nr;
    for(i=0; it != to; ++i, ++it){
      
      DEBUG(11, "\tnr " << i);
      
      std::string line(*it);
      if (line.length() < 17) {
        io::messages.add("line too short in POSRESSPEC block", "In_Posresspec",
                io::message::error);
      }
      
      // the first 17 chars are ignored
      line.erase(line.begin(), line.begin() + 17);
      
      _lineStream.clear();
      _lineStream.str(line);

      _lineStream >> nr;
      
      DEBUG(11, "\t" << nr);
      
      if(_lineStream.fail()){
        io::messages.add("bad line in POSRESSPEC block",
                "In_Posresspec", io::message::error);
        return;
      }
      
      if (nr < 1 || nr > topo.num_atoms()){
        io::messages.add("illegal atom in POSRESSPEC block",
                "In_Posresspec", io::message::error);
        return;
      }
      
      topo.position_restraints().push_back
      (topology::position_restraint_struct(nr-1, math::Vec(0.0)));
      
    }
  } // POSRESSPEC  
  
  if (!quiet){
    switch(sim.param().posrest.posrest){
      case simulation::posrest_off :
	os << "\tPosition restraints OFF\n";
	// how did you get here?
	break;
      case simulation::posrest_on:
      case simulation::posrest_bfactor:
	os << "\tPosition restraints ON\n"
	   << "\t\trestraining to ";
        if (sim.param().posrest.read) 
          os <<  "the positions in the positions restraints file.\n";
        else
	  os << "their positions in the startup file.\n";
	break;
      case simulation::posrest_const:
	os << "\tPosition constraints ON\n"
	   << "\t\tconstraining following atoms to ";
        if (sim.param().posrest.read) 
          os <<  "the positions given in the positions restraints file.\n";
        else
	  os << "their positions in the startup file.\n";
	break;
    }
    
    if (sim.param().posrest.scale_reference_positions)
      os << "\n\tReference positions are scaled upon pressure scaling.\n";
    
    os << "END\n";
  }
}

void 
io::In_Posres::read(topology::Topology& topo,
		    simulation::Simulation & sim,
		    std::ostream & os){
  std::vector<std::string> buffer;
  std::vector<std::string>::const_iterator it;
  
  // read the reference positions from this file
  { // REFPOSITION
    buffer = m_block["REFPOSITION"];
    if (sim.param().posrest.read) {
      if (!buffer.size()) {
        io::messages.add("no REFPOSITION block in position restraints "
                         "specification file", "In_Posres", io::message::error);
        return;
      } else {
        io::In_Configuration::_read_refposition(topo.position_restraints(), 
                                                buffer, true);
      }
    } else {
      if (buffer.size()) {
        io::messages.add("REFPOSITION block ignored in position restraints "
                         "specification file", "In_Posres", io::message::warning);
      }      
    }
  } // REFPOSITION
  
  { // BFACTOR
    buffer = m_block["BFACTOR"];
    if (sim.param().posrest.read && 
        sim.param().posrest.posrest == simulation::posrest_bfactor) {
      if (!buffer.size()) {
        io::messages.add("no BFACTOR block in position restraints "
                         "specification file", "In_Posres", io::message::error);
        return;
      } else {
        io::In_Configuration::_read_bfactor(topo.position_restraints(),
                                            buffer, true);
      }
    } else {
      if (buffer.size()) {
        io::messages.add("BFACTOR block ignored in position restraints "
                         "specification file", "In_Posres", io::message::warning);
      }      
    }  
  } // BFACTOR
}
