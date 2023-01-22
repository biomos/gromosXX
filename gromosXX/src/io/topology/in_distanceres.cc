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


void
io::In_Distanceres::read(topology::Topology& topo,
              simulation::Simulation & sim,
              std::ostream & os){

  DEBUG(7, "reading in a distance restraints file");

  if (!quiet)
    os << "DISTANCE RESTRAINTS\n";

  read_DISTANCERESSPEC(topo, sim, os);
  read_PERTDISRESSPEC(topo, sim, os);
  read_DFRESSPEC(topo, sim, os);
  read_PERTDFRESSPEC(topo, sim, os);
  read_MDISRESSPEC(topo, sim, os);

  if (!quiet) os << "END\n";

  if (sim.param().distanceres.distanceres) {
    if (block_read.count("DISTANCERESSPEC") == 0
     && block_read.count("PERTDISRESSPEC") == 0
     && block_read.count("MDISRESSPEC") == 0)
      io::messages.add("Distance restraints are on but neither DISTANCERESSPEC nor PERTDISRESSPEC nor MDISRESSPEC found",
               "In_Distanceres",
               io::message::error);
  }
  if (sim.param().distancefield.distancefield) {
    if (block_read.count("DFRESSPEC") == 0
     && block_read.count("PERTDFRESSPEC") == 0)
      io::messages.add("DF restraints on but neither DFRESSPEC nor PERTDFRESSPEC found",
               "In_Distanceres",
               io::message::error);
  }

  for(std::map<std::string, std::vector<std::string> >::const_iterator
    it = m_block.begin(),
    to = m_block.end();
      it != to;
      ++it){

    if (block_read.count(it->first) == 0 && it->second.size()){
      io::messages.add(" block " + it->first + " not supported!",
               "In_Distanceres",
               io::message::warning);
    }
  }

}

/**
 * @section distanceresspec DISTANCERESSPEC block
 * The DISTANCERESSPEC block is read from the distance restraints specification
 * file.
 *
 * \c DISH is the carbon-hydrogen, \c DISC the carbon-carbon distance.
 * See @ref util::virtual_type for valid virtual atom types.
 * \c r0 is the restraint distance, \c w0 a weight factor (multiplied by the force 
 * constant specified in the input file, \c CDIR) and rah specifies the type of 
 * restraint (half harmonic repulsive, full harmonic, half harmonic attractive). 
 * The restraint may be applied in a reduced set of dimensions, which is also set
 * by the value of \c rah. Allowed values of rah are \c dim - 1, \c dim or 
 * \c dim + 1, where dim can take the following values:
 * - dim = 0  : dimensions to apply distance restraint: X, Y, Z
 * - dim = 10 : dimensions to apply distance restraint: X, Y
 * - dim = 20 : dimensions to apply distance restraint: X, Z
 * - dim = 30 : dimensions to apply distance restraint: Y, Z
 * - dim = 40 : dimension to apply distance restraint: X
 * - dim = 50 : dimension to apply distance restraint: Y
 * - dim = 60 : dimension to apply distance restraint: Z
 *
 * The type of restraint is determined as follows:
 * - rah = dim - 1: half harmonic repulsive distance restraint
 * - rah = dim: full harmonic distance restraint
 * - rah = dim + 1: half harmonic attractive distance restraint
 *
 * @todo add documentation for further rah values for restraining only x,y,z or combinations
 *
 * @snippet snippets/snippets.cc DISTANCERESSPEC
 * @sa util::virtual_type util::Virtual_Atom
 */
void io::In_Distanceres::read_DISTANCERESSPEC(topology::Topology &topo,
                      simulation::Simulation &sim,
                      std::ostream & os)  { // DISTANCERES
  DEBUG(10, "DISTANCERESSPEC block");

  std::stringstream exampleblock;
  // lines starting with 'exampleblock<<"' and ending with '\n";' (spaces don't matter)
  // will be used to generate snippets that can be included in the doxygen doc;
  // the first line is the tag
  exampleblock << "DISTANCERESSPEC\n";
  exampleblock << "# DISH, DISC carbon-hydrogen/carbon-carbon distance\n";
  exampleblock << "# i,j,k,l  atoms comprising the virtual atom (put 0 if less than four atoms in use)\n";
  exampleblock << "# type   virtual atom type\n";
  exampleblock << "# r0, w0  target distance and force constant weighting factor\n";
  exampleblock << "# rah    form and dimension of the potential\n";
  exampleblock << "# full harmonic:\n";
  exampleblock << "#     0: x,y,z\n";
  exampleblock << "#    10: x,y\n";
  exampleblock << "#    20: x,z\n";
  exampleblock << "#    30: y,z\n";
  exampleblock << "#    40: x\n";
  exampleblock << "#    50: y\n";
  exampleblock << "#    60: z\n";
  exampleblock << "#  subtract or add 1 from these numbers to select a half harmonic\n";
  exampleblock << "#  repulsive or attractive potential\n";
  exampleblock << "# DISH  DISC\n";
  exampleblock << "  0.1   0.153\n";
  exampleblock << "# i  j  k  l  type    i  j  k  l  type    r0    w0    rah\n";
  exampleblock << "  1  0  0  0  0       10 12 11 13 3       0.2   1.0   0\n";
  exampleblock << "END\n";

  std::string blockname = "DISTANCERESSPEC";
  Block block(blockname, exampleblock.str());

  if (block.read_buffer(m_block[blockname], false) == 0) {
    block_read.insert(blockname);

    double dish = 0.0,disc = 0.0;

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

    block.get_next_parameter("DISH", dish, ">0", "");
    block.get_next_parameter("DISC", disc, ">0", "");

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

    unsigned int num = block.numlines()-3;
    for(unsigned int line_number=0; line_number < num; line_number++){

      DEBUG(11, "\tnr " << line_number);

      int type1 = 0, type2 = 0;
      std::vector<int> atom1, atom2;
      double r0 = 0.0,w0 = 0.0;
      int rah = 0;

      for(unsigned int i = 0; i < io::In_Distanceres::MAX_ATOMS; i++) {
        unsigned int atom = 0;
        std::string str_i = io::to_string(i);
        block.get_next_parameter("ATOM["+str_i+"]", atom, ">=0", "");
        if (atom > topo.num_atoms()) {
          std::ostringstream msg;
          msg << blockname << " block: atom number out of range: " << atom << ", last atom is "<< topo.num_atoms();
          io::messages.add(msg.str(), "In_Distanceres", io::message::error);
        }

        // -1 because we directly convert to array indices
        if (atom > 0) atom1.push_back(atom - 1);
      }
      block.get_next_parameter("TYPE", type1, "", "-2,-1,0,1,2,3,4,5,6,7");

      for(unsigned int i = 0; i < io::In_Distanceres::MAX_ATOMS; i++) {
        unsigned int atom = 0;
        std::string str_i = io::to_string(i);
        block.get_next_parameter("ATOM["+str_i+"]", atom, ">=0", "");
        if (atom > topo.num_atoms()) {
          std::ostringstream msg;
          msg << blockname << " block: atom number out of range: " << atom << ", last atom is "<< topo.num_atoms();
          io::messages.add(msg.str(), "In_Distanceres", io::message::error);
        }

        // -1 because we directly convert to array indices
        if (atom > 0) atom2.push_back(atom - 1);
      }
      block.get_next_parameter("TYPE", type2, "", "-2,-1,0,1,2,3,4,5,6,7");

      block.get_next_parameter("R0", r0, ">=0", "");
      block.get_next_parameter("W0", w0, ">=0", "");
      block.get_next_parameter("RAH", rah, "", "");

    // g++ 3.2 fix
      if(!block.error()){
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
    } // for restraints
    block.get_final_messages();
  } // if block content
} // DISTANCERES

/**
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
 * @snippet snippets/snippets.cc PERTDISRESSPEC
 */
void io::In_Distanceres::read_PERTDISRESSPEC(topology::Topology &topo,
                      simulation::Simulation &sim,
                      std::ostream & os)   { // PERTDISRESPEC
  DEBUG(10, "PERTDISRESSPEC block");

  std::stringstream exampleblock;
  // lines starting with 'exampleblock<<"' and ending with '\n";' (spaces don't matter)
  // will be used to generate snippets that can be included in the doxygen doc;
  // the first line is the tag
  exampleblock << "PERTDISRESSPEC\n";
  exampleblock << "# DISH  DISC\n";
  exampleblock << "  0.1   0.153\n";
  exampleblock << "# i  j  k  l  type    i  j  k  l  type n m   A_r0  A_w0  B_r0   B_w0  rah\n";
  exampleblock << "  1  0  0  0  0       10 12 11 13 3    1 1    0.2   1.0   0.5    2.0   0\n";
  exampleblock << "END\n";


  std::string blockname = "PERTDISRESSPEC";
  Block block(blockname, exampleblock.str());

  if (block.read_buffer(m_block[blockname], false) == 0) {
    block_read.insert(blockname);

    if (!sim.param().perturbation.perturbation) {
      io::messages.add("Perturbation is off but found a "+blockname+" block",
             "in_distanceres", io::message::warning);
      return;
    }
    sim.param().perturbation.perturbed_par=true;

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

    double dish = 0.0, disc = 0.0;
    block.get_next_parameter("DISH", dish, ">0", "");
    block.get_next_parameter("DISC", disc, ">0", "");

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

    unsigned int num = block.numlines()-3;
    for(unsigned int line_number=0; line_number < num; ++line_number){

      DEBUG(11, "\tnr " << line_number-2);

      int type1 = 0, type2 = 0;
      int n = 0,m = 0;
      std::vector<int> atom1, atom2;
      double A_r0 = 0.0, B_r0 = 0.0, A_w0 = 0.0, B_w0 = 0.0;
      int rah = 0;

      for(unsigned int i = 0; i < io::In_Distanceres::MAX_ATOMS; i++) {
        unsigned int atom = 0;
        std::string str_atom = io::to_string(i);
        block.get_next_parameter("ATOM["+str_atom+"]", atom, ">=0", "");
        if (atom > topo.num_atoms()) {
          std::ostringstream msg;
          msg << blockname << " block: atom number out of range: " << atom << ", last atom is "<< topo.num_atoms();
          io::messages.add(msg.str(), "In_Distanceres", io::message::error);
        }

        if (atom > 0) atom1.push_back(atom - 1);
      }
      block.get_next_parameter("TYPE", type1, "", "-2,-1,0,1,2,3,4,5,6,7");

      for(unsigned int i = 0; i < io::In_Distanceres::MAX_ATOMS; i++) {
        unsigned int atom = 0;
        std::string str_atom = io::to_string(i);
        block.get_next_parameter("ATOM["+str_atom+"]", atom, ">=0", "");
        if (atom > topo.num_atoms()) {
          std::ostringstream msg;
          msg << blockname << " block: atom number out of range: " << atom << ", last atom is "<< topo.num_atoms();
          io::messages.add(msg.str(), "In_Distanceres", io::message::error);
        }

        if (atom > 0) atom2.push_back(atom - 1);
      }
      block.get_next_parameter("TYPE", type2, "", "-2,-1,0,1,2,3,4,5,6,7");

      block.get_next_parameter("N", n, "", "");
      block.get_next_parameter("M", m, "", "");
      block.get_next_parameter("A_R0", A_r0, ">=0", "");
      block.get_next_parameter("A_W0", A_w0, ">=0", "");
      block.get_next_parameter("B_R0", B_r0, ">=0", "");
      block.get_next_parameter("B_W0", B_w0, ">=0", "");
      block.get_next_parameter("RAH", rah, "", "");


      if( !block.error()){
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
    block.get_final_messages();
  }
}//PERTDISRESPEC


/**
 * @section dfresspec DFRESSPEC block
 * The DFRESSPEC block is read from the distance restraints specification
 * file and used for distancefield restraints.
 *
 * See:
 * - A. de Ruiter and C. Oostenbrink, Protein-ligand binding from distancefield
 *   distances and Hamiltonian replica exchange simulations, J. Chem. Theory
 *   Comp. 9 (2013) 883 - 892, doi: 10.1021/ct300967a
 *
 * @snippet snippets/snippets.cc DFRESSPEC
 */
void io::In_Distanceres::read_DFRESSPEC(topology::Topology &topo,
                      simulation::Simulation &sim,
                      std::ostream & os) { // DISTANCEFIELD RES
  DEBUG(10, "DFRESSPEC block");

  std::stringstream exampleblock;
  // lines starting with 'exampleblock<<"' and ending with '\n";' (spaces don't matter)
  // will be used to generate snippets that can be included in the doxygen doc;
  // the first line is the tag
  exampleblock << "DFRESSPEC\n";
  exampleblock << "#   DISH H-C bond length for virtual atoms\n";
  exampleblock << "#   DISC C-C bond length for virtual atoms\n";
  exampleblock << "#   PROTEINATOMS > 0 last atom of the host\n";
  exampleblock << "#   K >= 0.0 Force constant\n";
  exampleblock << "#   r0 >=0 zero energy distance\n";
  exampleblock << "#   TYPE_I Virtual atom type for interaction site I\n";
  exampleblock << "#   NUM_I  Number of atoms defining interaction site I\n";
  exampleblock << "#   ATOM_I[0..NUM_I] Index numbers of atoms defining interaction site I\n";
  exampleblock << "#   TYPE_J Virtual atom type for interaction site J\n";
  exampleblock << "#   NUM_J  Number of atoms defining interaction site J\n";
  exampleblock << "#   ATOM_J[0..NUM_J] Index numbers of atoms defining interaction site J\n";
  exampleblock << "# DISH  DISC\n";
  exampleblock << "  0.1   0.153\n";
  exampleblock << "# PROTEINATOMS  K    r0\n";
  exampleblock << "  1190          500  0.0\n";
  exampleblock << "# TYPE_I  NUM_I  ATOM_I[0] .. ATOM_I[NUM_I]\n";
  exampleblock << "  -1      7        16  190  249  312  486  632 1208\n";
  exampleblock << "# TYPE_J  NUM_J  ATOM_J[0] .. ATOM_J[NUM_J]\n";
  exampleblock << "  -1      2      1194 1203\n";
  exampleblock << "END\n";

  std::string blockname = "DFRESSPEC";
  Block block(blockname, exampleblock.str());

  if (block.read_buffer(m_block[blockname], false) == 0) {
    block_read.insert(blockname);

    DEBUG(10, "reading in DFRES data");

    if (!quiet) {
    switch (sim.param().distancefield.distancefield) {
      case 0:
        os << "\tDistancefield restraints OFF\n";
        io::messages.add("Distancefield is off but found a "+blockname+" block",
             "in_distanceres", io::message::warning);
        return;
      case 1:
        os << "\tDistancefield restraints ON\n";
        break;
      default:
        os << "\tDistancefield restraints ERROR\n";
    }
    }

    double dish = 0.0,disc = 0.0;
    block.get_next_parameter("DISH", dish, ">0", "");
    block.get_next_parameter("DISC", disc, ">0", "");

    if (!quiet){
          os << std::setw(10) << "DISH"
          << std::setw(10) << "DISC"
          << "\n"
          <<  std::setw(10)<< dish
          <<  std::setw(10)<< disc
          << "\n";
    }

    int vtype_i = 0, vtype_j = 0;
    unsigned int num = 0, atom = 0;
    std::vector<int> atomi, atomj;

    topo.disfield_restraints().on = true;

    block.get_next_parameter("PROTEINATOMS", topo.disfield_restraints().proteinatoms, ">=0", "");
    if (topo.disfield_restraints().proteinatoms > topo.num_solute_atoms()) {
          std::ostringstream msg;
          msg << blockname << " block: PROTEINATOMS out of range: "
              << topo.disfield_restraints().proteinatoms
              << ", last solute atom is "<< topo.num_solute_atoms();
          io::messages.add(msg.str(), "In_Distanceres", io::message::error);
    }
    topo.disfield_restraints().proteinatoms--;
    block.get_next_parameter("K", topo.disfield_restraints().K, ">=0", "");
    block.get_next_parameter("r0", topo.disfield_restraints().r0, ">=0", "");

    block.get_next_parameter("TYPE_I", vtype_i, "", "");
    block.get_next_parameter("NUM_I", num, ">0", "");
    for(unsigned int i=0; i< num; i++){
      block.get_next_parameter("ATOM", atom, ">0", "");
      if (atom > topo.num_atoms()) {
          std::ostringstream msg;
          msg << blockname << " block: atom number out of range: " << atom << ", last atom is "<< topo.num_atoms();
          io::messages.add(msg.str(), "In_Distanceres", io::message::error);
      }
      atomi.push_back(atom-1);
    }
    block.get_next_parameter("TYPE_J", vtype_j, "", "");
    block.get_next_parameter("NUM_J", num, ">0", "");
    for(unsigned int i=0; i< num; i++){
      block.get_next_parameter("ATOM", atom, ">0", "");
      if (atom > topo.num_atoms()) {
          std::ostringstream msg;
          msg << blockname << " block: atom number out of range: " << atom << ", last atom is "<< topo.num_atoms();
          io::messages.add(msg.str(), "In_Distanceres", io::message::error);
      }
      atomj.push_back(atom-1);
    }

    if (!block.error()) {
      util::virtual_type t1 = util::virtual_type(vtype_i);
      util::virtual_type t2 = util::virtual_type(vtype_j);

      util::Virtual_Atom v1(t1, atomi, dish, disc);
      util::Virtual_Atom v2(t2, atomj, dish, disc);

      topo.disfield_restraints().v1 = v1;
      topo.disfield_restraints().v2 = v2;
    }

    block.get_final_messages();
  }

} // DISTANCEFIELD

/**
 * @section pertdfresspec PERTDFRESSPEC block
 * The PERTDFRESSPEC block is read from the distance restraints specification
 * file and used for perturbed distancefield restraints.
 *
 * See:
 * - A. de Ruiter and C. Oostenbrink, Protein-ligand binding from distancefield
 *   distances and Hamiltonian replica exchange simulations, J. Chem. Theory
 *   Comp. 9 (2013) 883 - 892, doi: 10.1021/ct300967a
 *
 * @snippet snippets/snippets.cc PERTDFRESSPEC
 */
void io::In_Distanceres::read_PERTDFRESSPEC(topology::Topology &topo,
                      simulation::Simulation &sim,
                      std::ostream & os) { // PERTDFRESPEC DISTANCERES
  DEBUG(10, "PERTDFRESSPEC block");

  std::stringstream exampleblock;
  // lines starting with 'exampleblock<<"' and ending with '\n";' (spaces don't matter)
  // will be used to generate snippets that can be included in the doxygen doc;
  // the first line is the tag
  exampleblock << "PERTDFRESSPEC\n";
  exampleblock << "#   DISH H-C bond length for virtual atoms\n";
  exampleblock << "#   DISC C-C bond length for virtual atoms\n";
  exampleblock << "#   PROTEINATOMS > 0 last atom of the host\n";
  exampleblock << "#   A_r0 >=0 reference distance for state A\n";
  exampleblock << "#   B_r0 >=0 reference distance for state B\n";
  exampleblock << "#   K_A >= 0 force constant state A\n";
  exampleblock << "#   K_B >= 0 force constant state B\n";
  exampleblock << "#   n >= 0 hidden restraint parameter n\n";
  exampleblock << "#   m >= 0 hidden restraint parameter m\n";
  exampleblock << "#   TYPE_I Virtual atom type for interaction site I\n";
  exampleblock << "#   NUM_I  Number of atoms defining interaction site I\n";
  exampleblock << "#   ATOM_I[0..NUM_I] Index numbers of atoms defining interaction site I\n";
  exampleblock << "#   TYPE_J Virtual atom type for interaction site J\n";
  exampleblock << "#   NUM_J  Number of atoms defining interaction site J\n";
  exampleblock << "#   ATOM_J[0..NUM_J] Index numbers of atoms defining interaction site J\n";
  exampleblock << "# DISH  DISC\n";
  exampleblock << "  0.1   0.153\n";
  exampleblock << "# PROTEINATOMS  A_r0  K_A  B_r0  K_B  n  m\n";
  exampleblock << "  1190          4.5   500  0.0   500  0  0\n";
  exampleblock << "# TYPE_I  NUM_I  ATOM_I[0] .. ATOM_I[NUM_I]\n";
  exampleblock << "  -1      7        16  190  249  312  486  632 1208\n";
  exampleblock << "# TYPE_J  NUM_J  ATOM_J[0] .. ATOM_J[NUM_J]\n";
  exampleblock << "  -1      2      1194 1203\n";
  exampleblock << "END\n";

  std::string blockname = "PERTDFRESSPEC";
  Block block(blockname, exampleblock.str());

  if (block.read_buffer(m_block[blockname], false) == 0) {
    block_read.insert(blockname);

    DEBUG(10, "reading in DISTANCERES (PERTDFRESSPEC) data");

    if (!quiet){
    switch(sim.param().distancefield.distancefield * sim.param().perturbation.perturbation){
      case 0:
        os << "\tPerturbed Distancefield restraints OFF\n";
        io::messages.add("No perturbation or no distancefield, but found a "+blockname+" block",
             "in_distanceres", io::message::warning);
        return;
      case 1:
        os << "\tPerturbed Distancefield restraints ON\n";
        break;
      default:
        os << "\tPerturbed Distancefield restraints ERROR\n";
    }
    }
    sim.param().perturbation.perturbed_par=true;

    double dish = 0.0,disc = 0.0;
    block.get_next_parameter("DISH", dish, ">0", "");
    block.get_next_parameter("DISC", disc, ">0", "");

    if (!quiet){
          os << std::setw(10) << "DISH"
          << std::setw(10) << "DISC"
          << "\n"
          <<  std::setw(10)<< dish
          <<  std::setw(10)<< disc
          << "\n";
    }

    int vtype_i = 0, vtype_j = 0;
    unsigned int num = 0, atom = 0;
    std::vector<int> atomi, atomj;

    topo.perturbed_disfield_restraints().on = true;

    block.get_next_parameter("PROTEINATOMS", topo.perturbed_disfield_restraints().proteinatoms, ">=0", "");
    if (topo.perturbed_disfield_restraints().proteinatoms > topo.num_solute_atoms()) {
          std::ostringstream msg;
          msg << blockname << " block: PROTEINATOMS out of range: "
              << topo.perturbed_disfield_restraints().proteinatoms
              << ", last solute atom is "<< topo.num_solute_atoms();
          io::messages.add(msg.str(), "In_Distanceres", io::message::error);
    }
    topo.perturbed_disfield_restraints().proteinatoms--;
    block.get_next_parameter("A_r0", topo.perturbed_disfield_restraints().A_r0, ">=0", "");
    block.get_next_parameter("K_A", topo.perturbed_disfield_restraints().K_A, ">=0", "");
    block.get_next_parameter("B_r0", topo.perturbed_disfield_restraints().B_r0, ">=0", "");
    block.get_next_parameter("K_B", topo.perturbed_disfield_restraints().K_B, ">=0", "");
    block.get_next_parameter("n", topo.perturbed_disfield_restraints().n, ">=0", "");
    block.get_next_parameter("m", topo.perturbed_disfield_restraints().m, ">=0", "");

    block.get_next_parameter("TYPE_I", vtype_i, "", "");
    block.get_next_parameter("NUM_I", num, ">0", "");
    for(unsigned int i=0; i< num; i++){
      block.get_next_parameter("ATOM", atom, ">0", "");
      if (atom > topo.num_atoms()) {
          std::ostringstream msg;
          msg << blockname << " block: atom number out of range: " << atom << ", last atom is "<< topo.num_atoms();
          io::messages.add(msg.str(), "In_Distanceres", io::message::error);
      }
      atomi.push_back(atom-1);
    }
    block.get_next_parameter("TYPE_J", vtype_j, "", "");
    block.get_next_parameter("NUM_J", num, ">0", "");
    for(unsigned int i=0; i< num; i++){
      block.get_next_parameter("ATOM", atom, ">0", "");
      if (atom > topo.num_atoms()) {
          std::ostringstream msg;
          msg << blockname << " block: atom number out of range: " << atom << ", last atom is "<< topo.num_atoms();
          io::messages.add(msg.str(), "In_Distanceres", io::message::error);
      }
      atomj.push_back(atom-1);
    }

    if (!block.error()) {
      util::virtual_type t1 = util::virtual_type(vtype_i);
      util::virtual_type t2 = util::virtual_type(vtype_j);

      util::Virtual_Atom v1(t1, atomi, dish, disc);
      util::Virtual_Atom v2(t2, atomj, dish, disc);

      topo.perturbed_disfield_restraints().v1 = v1;
      topo.perturbed_disfield_restraints().v2 = v2;
    }
    block.get_final_messages();
  }
}//PERTDFRESPEC DISTANCERES


/**
 * @section mdisresspec MDISRESSPEC block
 * The MDISRESSPEC block is read from the distance restraints specification
 * file and used for EDS restraints.
 *
 * The format is very similar to the @ref distanceresspec with the difference that
 * one may give values for multipde states.
 * @snippet snippets/snippets.cc MDISRESSPEC
 */
void io::In_Distanceres::read_MDISRESSPEC(topology::Topology &topo,
                      simulation::Simulation &sim,
                      std::ostream & os) { // MDISRESSPEC
  DEBUG(10, "MDISRESSPEC block");

  std::stringstream exampleblock;
  // lines starting with 'exampleblock<<"' and ending with '\n";' (spaces don't matter)
  // will be used to generate snippets that can be included in the doxygen doc;
  // the first line is the tag
  exampleblock << "MDISRESSPEC\n";
  exampleblock << "# DISH  DISC\n";
  exampleblock << "  0.1   0.153\n";
  exampleblock << "# N: number of eds states (3 in this example)\n";
  exampleblock << "# i  j  k  l  type    i  j  k  l  type    r0[1 ... N]    w0[1 ... N]    rah\n";
  exampleblock << "  1  0  0  0  0       10 12 11 13 3       0.2  0.2  0.2  1.0  0.0 0.0    0\n";
  exampleblock << "  5  0  0  0  0       10 12 11 13 3       0.2  0.2  0.2  0.0  1.0 0.0    0\n";
  exampleblock << "  8  0  0  0  0       10 12 11 13 3       0.2  0.2  0.2  0.0  0.0 1.0    0\n";
  exampleblock << "END\n";

  std::string blockname = "MDISRESSPEC";
  Block block(blockname, exampleblock.str());

  if (block.read_buffer(m_block[blockname], false) == 0) {
    block_read.insert(blockname);

    if (!sim.param().eds.eds){
        io::messages.add("MDISRESSPEC block given but EDS not turned on!",
                         "in_distanceres", io::message::error);
        return;
    }

    if (sim.param().distanceres.distanceres < 0){
        io::messages.add("eds perturbed distance restraints not compatible with time averaging!",
                         "in_distanceres", io::message::error);
    }

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

    double dish = 0.0,disc = 0.0;
    block.get_next_parameter("DISH", dish, ">0", "");
    block.get_next_parameter("DISC", disc, ">0", "");

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

    unsigned int num = block.numlines()-3;
    for(unsigned int line_number=0; line_number < num; line_number++){

      DEBUG(11, "\tnr " << line_number - 2);

      int type1 = 0, type2 = 0;
      std::vector<int> atom1, atom2;
      std::vector<double> r0(numstates),w0(numstates);
      int rah = 0;

      for(unsigned int i = 0; i < io::In_Distanceres::MAX_ATOMS; i++) {
        unsigned int atom = 0;
        std::string str_atom = io::to_string(i);
        block.get_next_parameter("ATOM["+str_atom+"]", atom, ">=0", "");
        if (atom > topo.num_atoms()) {
          std::ostringstream msg;
          msg << blockname << " block: atom number out of range: " << atom << ", last atom is "<< topo.num_atoms();
          io::messages.add(msg.str(), "In_Distanceres", io::message::error);
        }

        // -1 because we directly convert to array indices
        if (atom > 0) atom1.push_back(atom - 1);
      }
      block.get_next_parameter("TYPE", type1, "", "-2,-1,0,1,2,3,4,5,6,7");

      for(unsigned int i = 0; i < io::In_Distanceres::MAX_ATOMS; i++) {
        unsigned int atom = 0;
        std::string str_atom = io::to_string(i);
        block.get_next_parameter("ATOM["+str_atom+"]", atom, ">=0", "");
        if (atom > topo.num_atoms()) {
          std::ostringstream msg;
          msg << blockname << " block: atom number out of range: " << atom << ", last atom is "<< topo.num_atoms();
          io::messages.add(msg.str(), "In_Distanceres", io::message::error);
        }

        // -1 because we directly convert to array indices
        if (atom > 0) atom2.push_back(atom - 1);
      }
      block.get_next_parameter("TYPE", type2, "", "-2,-1,0,1,2,3,4,5,6,7");

      for(unsigned int i = 0; i < numstates; i++){
          std::string str_i = io::to_string(i);
          block.get_next_parameter("r0["+str_i+"]", r0[i], ">=0", "");
      }
      for(unsigned int i = 0; i < numstates; i++){
          std::string str_i = io::to_string(i);
          block.get_next_parameter("w0["+str_i+"]", w0[i], ">=0", "");
      }
      block.get_next_parameter("RAH", rah, "", "");

    // g++ 3.2 fix
      if(!block.error()){
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

    block.get_final_messages();
  }
} // MDISRESSPEC
