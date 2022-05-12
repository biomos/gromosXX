/**
 * @file in_angrest.cc
 * implements methods of In_Angrest
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

#include "in_angrest.h"

#undef MODULE
#undef SUBMODULE
#define MODULE io
#define SUBMODULE topology

static std::set<std::string> block_read;
void io::In_Angrest::read(topology::Topology& topo,
    simulation::Simulation & sim,
    std::ostream & os){

  DEBUG(7, "reading in an angle restraints file");

  std::ostringstream oss;
  oss << "the angle is calculated from the nearest-image\n"
       << "         vectors, you have to make sure that none of the vectors defining it become \n"
       << "         longer than half the box length at any point during the simulations!\n";
  io::messages.add(oss.str(), "in_angleres", io::message::warning);

  if (!quiet)
    os << "ANGLE RESTRAINTS\n";

  read_ANGRESSPEC(topo, sim, os);
  read_PERTANGRESSPEC(topo, sim, os);

  if (!quiet) os << "END\n";

  if (sim.param().angrest.angrest) {
    if (block_read.count("ANGRESSPEC") == 0
     && block_read.count("PERTANGRESSPEC") == 0)
      io::messages.add("ANGLE restraints are on but neither ANGRESSPEC nor PERTANGRESSPEC found",
               "In_Angrest",
               io::message::error);
  }

  for(std::map<std::string, std::vector<std::string> >::const_iterator
  it = m_block.begin(),
  to = m_block.end();
      it != to;
      ++it){

    if (block_read.count(it->first) == 0 && it->second.size()){
      io::messages.add("block " + it->first + " not supported!",
         "In_Angrest",
         io::message::warning);
    }
  }

}

/**
  * @section angresspec ANGRESSPEC block
  * This block is for unperturbed angle restraints/constraints. It is read
  * from the angle restraint specification file.
 * @snippet snippets/snippets.cc ANGRESSPEC
 */

void io::In_Angrest::read_ANGRESSPEC(topology::Topology &topo,
  simulation::Simulation &sim,
  std::ostream & os) { // ANGRESSPEC
  DEBUG(10, "ANGRESSPEC block");

    std::stringstream exampleblock;
    // lines starting with 'exampleblock<<"' and ending with '\n";' (spaces don't matter)
    // will be used to generate snippets that can be included in the doxygen doc;
    // the first line is the tag
   exampleblock << "ANGRESSPEC\n";
   exampleblock << "# IPLR, JPLR, KPLR         atom sequence numbers of the atoms defining the\n";
   exampleblock << "#                          restrained angle i-j-k\n";
   exampleblock << "# WALR                     individual angle restraint weight factor by which\n";
   exampleblock << "#                          the harmonic angle restraining term may be multiplied.\n";
   exampleblock << "# PALR                     angle value (in degrees) at minimum energy\n";
   exampleblock << "#                          of the harmonic angle restraining term.\n";
   exampleblock << "#  IPLR  JPLR  KPLR  WALR  PALR\n";
   exampleblock << "     1     2     3   1.0   120.0\n";
   exampleblock << "END\n";

  std::string blockname = "ANGRESSPEC";
  Block block(blockname, exampleblock.str());

  if (block.read_buffer(m_block[blockname], false) == 0) {
    block_read.insert(blockname);

    DEBUG(10, "reading in ANGREST data");

    if (!quiet){
      switch(sim.param().angrest.angrest){
        case simulation::angle_restr_off:
          os << "\tAngle restraints OFF\n";
          // how did you get here?
          break;
        case simulation::angle_restr_inst:
          os << "\tAngle restraints ON\n"
             << "\t\t(uniform force constant K)\n";
          break;
        case simulation::angle_restr_inst_weighted:
     	    os << "\tAngle restraints ON\n"
     		     << "\t\t(force constant K*w0)\n";
          break;
        case simulation::angle_constr:
          os << "\tAngle constraints ON\n";
          break;
        default:
          os << "\tAngle restraints: ERROR\n";
          io::messages.add("wrong value for method in angle restraints block",
     		    "in_angle", io::message::error);
          return;
      }
    }

  if (!quiet){
  os << std::setw(10) << "i"
     << std::setw(8) << "j"
     << std::setw(8) << "k"
     << std::setw(8) << "w0"
     << std::setw(8) << "theta"
     << "\n";

   os.precision(2);
   os.setf(std::ios::fixed, std::ios::floatfield);
  }


  int i = 0, j = 0, k = 0;
  double theta = 0.0, w0 = 0.0;

  unsigned int num = block.numlines()-2;
  for(unsigned int line_number=0; line_number < num; ++line_number){
    DEBUG(11, "\tnr " << line_number);
    block.get_next_parameter("IPLR", i, ">0", "");
    block.get_next_parameter("JPLR", j, ">0", "");
    block.get_next_parameter("KPLR", k, ">0", "");
    block.get_next_parameter("WALR", w0, ">=0", "");
    block.get_next_parameter("PALR", theta, "", "");

    // move theta into range 0 to 180 degrees
    double theta_input = theta;
    if (theta < 0) theta *= -1;
    while (theta > 360) theta -= 360;
    

  topo.angle_restraints().push_back
(topology::angle_restraint_struct(i-1, j-1, k-1, theta * math::Pi / 180, w0));

  if (sim.param().angrest.angrest == simulation::angle_constr) {
    // remove bond angle potentials for constrained angles
     std::vector<topology::three_body_term_struct>::iterator it=topo.solute().angles().begin(), 
     to=topo.solute().angles().end();
     for (;it != to; ++it) {
       if (it->i == i-1 && it->j == j-1 && it->k == k-1) {
         topo.solute().angles().erase(it);
         DEBUG(9, "removed angle " <<  i << ", number of angles " << topo.solute().angles().size());
         break;
       }
     }
  }

  if (!quiet){
      os << std::setw(10) << i
         << std::setw(8) << j
         << std::setw(8) << k
         << std::setw(8) << w0
         << std::setw(8) << theta;
      if (theta != theta_input) os << " # WARNING: angle was mapped to between 0 and 180 degrees";
      os << "\n";
  }
    }
    block.get_final_messages();
  } // if block empty or not there
} // ANGREST


/**
* @section pertangresspec PERTANGRESSPEC block
* This block is for perturbed angle restraints/constraints. It is read from
* the angle restraints specification file.
* @snippet snippets/snippets.cc PERTANGRESSPEC
*/
void io::In_Angrest::read_PERTANGRESSPEC(topology::Topology &topo,
  simulation::Simulation &sim,
  std::ostream & os) { // PERTANGRESPEC
  DEBUG(10, "PERTANGRESSPEC block");

  std::stringstream exampleblock;
  // lines starting with 'exampleblock<<"' and ending with '\n";' (spaces don't matter)
  // will be used to generate snippets that can be included in the doxygen doc;
  // the first line is the tag
  exampleblock << "PERTANGRESSPEC\n";
  exampleblock << "# IPLR, JPLR, KPLR       atom sequence numbers of the atoms defining the\n";
  exampleblock << "#                          restrained angle i-j-k-l\n";
  exampleblock << "# APALR    angle angle value (in degrees) at minimum energy of the harmonic\n";
  exampleblock << "#          angle restraining term in state A.\n";
  exampleblock << "# AWALR    individual angle restraint weight factor by which the harmonic\n";
  exampleblock << "#          angle restraining term may be multiplied in state A.\n";
  exampleblock << "# BPALR    as APALR, but for state B.\n";
  exampleblock << "# BWALR    as AWALR, but for state B.\n";
  exampleblock << "# M        hidden restraint parameter m and\n";
  exampleblock << "# N        hidden restraint parameter n of\n";
  exampleblock << "#          hidden restraint prefactor l^n*(1-l)^m.\n";
  exampleblock << "# IPLR  JPLR  KPLR   M   N  APALR  AWALR  BPALR  BWALR\n";
  exampleblock << "     1     2     3   2   1  120.0    1.0  160.0    1.0\n";
  exampleblock << "END\n";

  std::string blockname = "PERTANGRESSPEC";
  Block block(blockname, exampleblock.str());

  if (block.read_buffer(m_block[blockname], false) == 0) {
    block_read.insert(blockname);

    sim.param().perturbation.perturbed_par=true;

  DEBUG(10, "reading in perturbed ANGREST (PERTANGRESSPEC data");

  if (!quiet){
    switch(sim.param().angrest.angrest){
      case 0:
        os << "\tPerturbed Angle restraints OFF\n";
        // how did you get here?
        break;
      case 1:
        os << "\tPerturbed Angle restraints ON\n"
         << "\t\t(using uniform force constant K\n";
        break;
      case 2:
        os << "\tPerturbed Angle restraints ON\n"
         << "\t\t(using force constant K*w0)\n";
        break;
      case 3:
        os << "\tPerturbed Angle constraints ON\n";
        break;
      default:
        os << "\tPerturbed Angle restraints ERROR\n";
        io::messages.add("wrong method for angle restraints",
         "in_angrest", io::message::error);
        return;
    }
  }

  if (!quiet){
    os << std::setw(10) << "i"
      << std::setw(8) << "j"
      << std::setw(8) << "k"
      << std::setw(8) << "m"
      << std::setw(8) << "n"
      << std::setw(8) << "A_theta"
      << std::setw(8) << "A_w0"
      << std::setw(8) << "B_theta"
      << std::setw(8) << "B_w0"
      << "\n";

      os.precision(2);
      os.setf(std::ios::fixed, std::ios::floatfield);
  }

  int i = 0, j = 0, k = 0, m = 0, n = 0;
  double A_theta = 0.0, A_w0 = 0.0, B_theta = 0.0, B_w0 = 0.0;

  unsigned int num = block.numlines()-2;
  for(unsigned int line_number=0; line_number < num; ++line_number){

    DEBUG(11, "\tnr " << line_number);
    block.get_next_parameter("IPLR", i, ">0", "");
    block.get_next_parameter("JPLR", j, ">0", "");
    block.get_next_parameter("KPLR", k, ">0", "");
    block.get_next_parameter("M", m, ">0", "");
    block.get_next_parameter("N", n, ">0", "");
    block.get_next_parameter("APALR", A_theta, "", "");
    block.get_next_parameter("AWALR", A_w0, ">=0", "");
    block.get_next_parameter("BPALR", B_theta, "", "");
    block.get_next_parameter("BWALR", B_w0, ">=0", "");

    topo.perturbed_angle_restraints().push_back
(topology::perturbed_angle_restraint_struct(i-1, j-1, k-1, m, n,
   A_theta * math::Pi / 180, A_w0,
   B_theta * math::Pi / 180, B_w0 ));

    // move theta into range 0 to 180 degrees
    double A_theta_input = A_theta;
    if (A_theta < 0) A_theta*=-1;
    while (A_theta > 360) A_theta -= 360;
    double B_theta_input = B_theta;
    if (B_theta < 0) B_theta*=-1;
    while (B_theta > 360) B_theta -= 360;

      if (!quiet){
        os << std::setw(10) << i
          << std::setw(8) << j
          << std::setw(8) << k
          << std::setw(8) << m
          << std::setw(8) << n
          << std::setw(8) << A_theta
          << std::setw(8) << A_w0
          << std::setw(8) << B_theta
          << std::setw(8) << B_w0
          ;
      if (B_theta != B_theta_input || A_theta != A_theta_input) os << " # WARNING: angle was mapped to between 0 and 180 degrees";
      os << "\n";
      }
    }

    block.get_final_messages();
  }
} // PERTANGRESSPEC
