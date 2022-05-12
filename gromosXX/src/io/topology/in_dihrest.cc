/**
 * @file in_dihrest.cc
 * implements methods of In_Dihrest
 */

#include <sstream>
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
void io::In_Dihrest::read(topology::Topology& topo,
    simulation::Simulation & sim,
    std::ostream & os){

  DEBUG(7, "reading in a dihedral restraints file");

  std::ostringstream oss;
  oss << "the dihedral (specified in the DIHEDRALRESSPEC block) \n"
       << "         is calculated from the nearest-image\n"
       << "         vectors, you have to make sure that none of the vectors defining it become \n"
       << "         longer than half the box length at any point during the simulations!\n";
  io::messages.add(oss.str(), "in_dihedralres", io::message::warning);

  if (!quiet)
    os << "DIHEDRAL RESTRAINTS\n";

  read_DIHEDRALRESSPEC(topo, sim, os);
  read_PERTDIHRESSPEC(topo, sim, os);

  if (!quiet) os << "END\n";

  if (sim.param().dihrest.dihrest) {
    if (block_read.count("DIHEDRALRESSPEC") == 0
     && block_read.count("PERTDIHRESSPEC") == 0)
      io::messages.add("DIHEDRAL restraints are on but neither DIHEDRALRESSPEC nor PERTDIHRESSPEC found",
               "In_Dihrest",
               io::message::error);
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

/**
  * @section dihedralresspec DIHEDRALRESSPEC block
  * This block is for unperturbed dihedrals restraints/constraints. It is read
  * from the dihedral restraint specifcation file.
 * @snippet snippets/snippets.cc DIHEDRALRESSPEC
 */

void io::In_Dihrest::read_DIHEDRALRESSPEC(topology::Topology &topo,
  simulation::Simulation &sim,
  std::ostream & os) { // DIHEDRALRESSPEC
  DEBUG(10, "DIHEDRALRESSPEC block");

    std::stringstream exampleblock;
    // lines starting with 'exampleblock<<"' and ending with '\n";' (spaces don't matter)
    // will be used to generate snippets that can be included in the doxygen doc;
    // the first line is the tag
   exampleblock << "DIHEDRALRESSPEC\n";
   exampleblock << "# IPLR, JPLR, KPLR, LPLR   atom sequence numbers of the atoms defining the\n";
   exampleblock << "#                          restrained dihedral i-j-k-l\n";
   exampleblock << "# WDLR                     individual dihedral restraint weight factor by which\n";
   exampleblock << "#                          the harmonic dihedral restraining term may be multiplied.\n";
   exampleblock << "# PDLR                     dihedral angle value (in degrees) at minimum energy\n";
   exampleblock << "#                          of the harmonic dihedral restraining term.\n";
   exampleblock << "# DELTA                    deviation of the zero-energy dihedral angle value after\n";
   exampleblock << "#                          which the restraint switches to the next periodic\n";
   exampleblock << "#                          value. The dihedral angle is put in the interval\n";
   exampleblock << "#                          [ PDLR + DELTA - 360 , PDLR + DELTA ]\n";
   exampleblock << "#  IPLR  JPLR  KPLR  LPLR  WDLR  PDLR  DELTA\n";
   exampleblock << "     1     2     3     4    1.0   120.0 180.0\n";
   exampleblock << "END\n";

  std::string blockname = "DIHEDRALRESSPEC";
  Block block(blockname, exampleblock.str());

  if (block.read_buffer(m_block[blockname], false) == 0) {
    block_read.insert(blockname);

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
          os << "\tDihedral constraints (sine-and-cosine alg.) ON\n";
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
     << std::setw(8) << "w0"
     << std::setw(8) << "phi"
     << std::setw(8) << "delta"
     << "\n";

   os.precision(2);
   os.setf(std::ios::fixed, std::ios::floatfield);
  }


  int i = 0, j = 0, k = 0, l = 0;
  double delta = 0.0, phi = 0.0, w0 = 0.0;

  unsigned int num = block.numlines()-2;
  for(unsigned int line_number=0; line_number < num; ++line_number){
    DEBUG(11, "\tnr " << line_number);
    block.get_next_parameter("IPLR", i, ">0", "");
    block.get_next_parameter("JPLR", j, ">0", "");
    block.get_next_parameter("KPLR", k, ">0", "");
    block.get_next_parameter("LPLR", l, ">0", "");
    block.get_next_parameter("WDLR", w0, ">=0", "");
    block.get_next_parameter("PDLR", phi, "", "");
    block.get_next_parameter("DELTA", delta, "", "");

    // move phi into range -180 to 180 degrees
    double phi_input = phi;
    while (phi > 180) phi -= 360;
    while (phi < -180) phi += 360;
    

  topo.dihedral_restraints().push_back
(topology::dihedral_restraint_struct(i-1, j-1, k-1, l-1,
   delta * 2 * math::Pi / 360, phi * 2 * math::Pi / 360, w0));

  if (!quiet){
      os << std::setw(10) << i
         << std::setw(8) << j
         << std::setw(8) << k
         << std::setw(8) << l
         << std::setw(8) <<  w0
         << std::setw(8) << phi
         << std::setw(8) << delta;
      if (phi != phi_input) os << " # WARNING: angle was mapped to between -180 to 180 degrees";
      os << "\n";
  }
    }
    block.get_final_messages();
  } // if block empty or not there
} // DIHREST


/**
* @section pertdihresspec PERTDIHRESSPEC block
* This block is for perturbed dihedral restraints/constraints. It is read from
* the dihedral restraints specification file.
* @snippet snippets/snippets.cc PERTDIHRESSPEC
*/
void io::In_Dihrest::read_PERTDIHRESSPEC(topology::Topology &topo,
  simulation::Simulation &sim,
  std::ostream & os) { // PERTDIHRESPEC
  DEBUG(10, "PERTDIHRESSPEC block");

  std::stringstream exampleblock;
  // lines starting with 'exampleblock<<"' and ending with '\n";' (spaces don't matter)
  // will be used to generate snippets that can be included in the doxygen doc;
  // the first line is the tag
  exampleblock << "PERTDIHRESSPEC\n";
  exampleblock << "# IPLR, JPLR, KPLR, LPLR   atom sequence numbers of the atoms defining the\n";
  exampleblock << "#                          restrained dihedral i-j-k-l\n";
  exampleblock << "# APDLR    dihedral angle value (in degrees) at minimum energy of the harmonic\n";
  exampleblock << "#          dihedral restraining term in state A.\n";
  exampleblock << "# AWDLR    individual dihedral restraint weight factor by which the harmonic\n";
  exampleblock << "#          dihedral restraining term may be multiplied in state A.\n";
  exampleblock << "# BPDLR    as APDLR, but for state B.\n";
  exampleblock << "# BWDLR    as AWDLR, but for state B.\n";
  exampleblock << "# M        hidden restraint parameter m and\n";
  exampleblock << "# N        hidden restraint parameter n of\n";
  exampleblock << "#          hidden restraint prefactor l^n*(1-l)^m.\n";
  exampleblock << "# DELTA    deviation of the zero-energy dihedral angle value after which the\n";
  exampleblock << "#          restraint switches to the next periodic value. The dihedral angle\n";
  exampleblock << "#          is put in the interval\n";
  exampleblock << "#          [ (1-RLAM)*APDLR + RLAM*BPDLR + DELTA - 360 ,\n";
  exampleblock << "#                                       (1-RLAM)*APDLR + RLAM*BPDLR + DELTA ]\n";
  exampleblock << "#\n";
  exampleblock << "# IPLR  JPLR  KPLR  LPLR  M   N  DELTA  APDLR  AWDLR  BPDLR  BWDLR\n";
  exampleblock << "     1     2     3     4  2   1  180.0  120.0    1.0  160.0    1.0\n";
  exampleblock << "END\n";

  std::string blockname = "PERTDIHRESSPEC";
  Block block(blockname, exampleblock.str());

  if (block.read_buffer(m_block[blockname], false) == 0) {
    block_read.insert(blockname);

    sim.param().perturbation.perturbed_par=true;

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

  int i = 0, j = 0, k = 0, l = 0, m = 0, n = 0;
  double delta = 0.0, A_phi = 0.0, A_w0 = 0.0, B_phi = 0.0, B_w0 = 0.0;

  unsigned int num = block.numlines()-2;
  for(unsigned int line_number=0; line_number < num; ++line_number){

    DEBUG(11, "\tnr " << line_number);
    block.get_next_parameter("IPLR", i, ">0", "");
    block.get_next_parameter("JPLR", j, ">0", "");
    block.get_next_parameter("KPLR", k, ">0", "");
    block.get_next_parameter("LPLR", l, ">0", "");
    block.get_next_parameter("M", m, ">0", "");
    block.get_next_parameter("N", n, ">0", "");
    block.get_next_parameter("DELTA", delta, "", "");
    block.get_next_parameter("APDLR", A_phi, "", "");
    block.get_next_parameter("AWDLR", A_w0, ">=0", "");
    block.get_next_parameter("BPDLR", B_phi, "", "");
    block.get_next_parameter("BWDLR", B_w0, ">=0", "");

    topo.perturbed_dihedral_restraints().push_back
(topology::perturbed_dihedral_restraint_struct(i-1, j-1, k-1, l-1, m, n, delta * 2 * math::Pi / 360,
   A_phi * 2 * math::Pi / 360, A_w0,
   B_phi * 2 * math::Pi / 360, B_w0 ));

    // move phi into range -180 to 180 degrees
    double A_phi_input = A_phi;
    while (A_phi > 180) A_phi -= 360;
    while (A_phi < -180) A_phi += 360;
    double B_phi_input = B_phi;
    while (B_phi > 180) B_phi -= 360;
    while (B_phi < -180) B_phi += 360;

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
          ;
      if (B_phi != B_phi_input || A_phi != A_phi_input) os << " # WARNING: angle was mapped to between -180 to 180 degrees";
      os << "\n";
      }
    }

    block.get_final_messages();
  }
} // PERTDIHRESSPEC
