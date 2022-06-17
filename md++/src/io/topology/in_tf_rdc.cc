/**
 * @file in_tf_rdc.cc
 * implements methods of In_Tfrdcresspec
 */

#include <stdheader.h>

#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <simulation/simulation.h>
#include <interaction/interaction_types.h>
#include <configuration/configuration.h>

#include "../../io/instream.h"
#include "../../io/blockinput.h"
#include "in_tf_rdc.h"

#undef MODULE
#undef SUBMODULE
#define MODULE io
#define SUBMODULE topology

static std::set<std::string> block_read;

void
io::In_Tfrdcresspec::read(topology::Topology& topo,
        simulation::Simulation & sim,
        std::ostream & os) {
  std::set<std::string> block_read;
  DEBUG(7, "reading in a tensor-free RDC restraints file");

  read_TFRDCCONV(topo, sim, os);
  read_TFRDCRESSPEC(topo, sim, os);
  read_TFRDCMOLAXIS(topo, sim, os);

}

/**
 * @section tfrdcconversion TFRDCCONV block
 * The TFRDCCONV block is read from the TF-RDC restraint specification file.
 *
 * - Frequency unit conversion factor to ps-1.
 * - Gyr.magn. ratio unit conversion factor to (e/u).
 * @snippet snippets/snippets.cc TFRDCCONV
 */
void io::In_Tfrdcresspec::read_TFRDCCONV(topology::Topology &topo,
                                 simulation::Simulation &sim,
                                 std::ostream &os)
{ // TFRDCCONV

  DEBUG(10, "TFRDCCONV block");
  std::stringstream exampleblock;
  // lines starting with 'exampleblock<<"' and ending with '\n";' (spaces don't matter)
  // will be used to generate snippets that can be included in the doxygen doc;
  // the first line is the tag
  exampleblock << "TFRDCCONV\n";
  exampleblock << "# FACFREQ   frequency conversion factor from input units to ps-1\n";
  exampleblock << "#           typically input is Hz and the factor is 10^-12\n";
  exampleblock << "# FACGYR    rgyr conversion factor from input units to e/u\n";
  exampleblock << "#           typically input is 10^6*rad/T*s and the factor is 0.010375\n";
  exampleblock << "#  FACFREQ        FACGYR\n";
  exampleblock << "   0.000000000001 0.010364272\n";
  exampleblock << "END\n";

  std::string blockname = "TFRDCCONV";
  Block block(blockname, exampleblock.str());

  if (block.read_buffer(m_block[blockname], true) == 0)
  {
    block_read.insert(blockname);

    DEBUG(10, "reading in TF-RDC TFRDCCONV block")
    block.get_next_parameter("FACFREQ", sim.param().tfrdc.factorFreq, ">0", "");
    block.get_next_parameter("FACGYR", sim.param().tfrdc.factorGyr, ">0", "");

    os.precision(12);
    os.setf(std::ios_base::fixed, std::ios_base::floatfield);

    DEBUG(10, scientific << setprecision(6) << setw(14) << "factorFreq: " << sim.param().tfrdc.factorFreq)
    DEBUG(10, scientific << setprecision(6) << setw(14) << "factorGyr: " << sim.param().tfrdc.factorGyr)

    block.get_final_messages();
  } // if block empty or not there
} // TFRDCCONV

/**
 * @section tfrdcresspec TFRDCRESSPEC block
 * The TFRDCRESSPEC block is read from the tensor-free RDC restraints
 *  specification file .tfr.
 * @snippet snippets/snippets.cc TFRDCRESSPEC
 */
void io::In_Tfrdcresspec::read_TFRDCRESSPEC(topology::Topology &topo,
                                 simulation::Simulation &sim,
                                 std::ostream &os)
{ // TFRDCRESSPEC
    DEBUG(10, "TFRDCRESSPEC block");
  std::stringstream exampleblock;
  // lines starting with 'exampleblock<<"' and ending with '\n";' (spaces don't matter)
  // will be used to generate snippets that can be included in the doxygen doc;
  // the first line is the tag
  exampleblock << "TFRDCRESSPEC\n";
  exampleblock << "# DISH, DISC carbon-hydrogen/carbon-carbon distance\n";
  exampleblock << "# i,j,k,l  atoms comprising the virtual atom (fill with 0 if less than 4 atoms used)\n";
  exampleblock << "# type     virtual atom type\n";
  exampleblock << "# R0       normalisation distance\n";
  exampleblock << "# G1, G2   gyromagnetic ratio for nucleus 1/nucleus 2\n";
  exampleblock << "# D0       target RDC\n";
  exampleblock << "# DD0      half the width of the flatbottom potential\n";
  exampleblock << "# WTFRDC   a weight factor (multiplied by the force constant \n";
  exampleblock << "#          specified in the input file)\n";
  exampleblock << "# DISH  DISC\n";
  exampleblock << "  0.1   0.153\n";
  exampleblock << "# i   j   k   l   type  i    j    k   l  type  R0      G1      G2      D0      DD0   WTFRDC\n";
  exampleblock << "  2   0   0   0   0     4    0    0   0    0   0.1  267.513  19.331   0.136    0.05     1.0\n";
  exampleblock << "  3   0   0   0   0     5    0    0   0    0   0.1  267.513  67.262   0.136    0.05     1.0\n";
  exampleblock << "END\n";

  std::string blockname = "TFRDCRESSPEC";
  Block block(blockname, exampleblock.str());

  if (block.read_buffer(m_block[blockname], true) == 0)
  {
    block_read.insert(blockname);

    block.get_next_parameter("DISH", dish, ">0", "");
    block.get_next_parameter("DISC", disc, ">0", "");

    DEBUG(10, "reading in TFRDCRESSPEC data")

    os.setf(std::ios_base::fixed, std::ios_base::floatfield);

    DEBUG(10, setw(6) << "i"
                      << setw(6) << "j"
                      << setw(6) << "k"
                      << setw(6) << "l"
                      << setw(6) << "type"
                      << setw(6) << "i"
                      << setw(6) << "j"
                      << setw(6) << "k"
                      << setw(6) << "l"
                      << setw(6) << "type"
                      << setw(19) << "R0"
                      << setw(12) << "G1 [e/u]"
                      << setw(12) << "G2 [e/u]"
                      << setw(8) << "D0 [1/ps]"
                      << setw(8) << "DD0 [1/ps]"
                      << setw(8) << "WTFRDC");

    unsigned int num = block.numlines() - 3;
    for (unsigned int line_number = 0; line_number < num; line_number++)
    {
      std::vector<int> atom1, atom2;
      int type1, type2;
      double weight, D0, DD0, gyr1, gyr2, r0;

      DEBUG(11, "\ttfrdc line " << line_number);

      for (unsigned int i = 0; i < 4; i++)
      {
        unsigned int atom;
        std::string str_i = io::to_string(i);
        block.get_next_parameter("ATOM[" + str_i + "]", atom, ">=0", "");
        DEBUG(11, "\tat " << i << " " << atom);
        if (atom > topo.num_atoms())
        {
          std::ostringstream msg;
          msg << blockname << " block: atom number out of range: " << atom << ", last atom is " << topo.num_atoms();
          io::messages.add(msg.str(), "In_TFRDC", io::message::error);
        }

        // -1 because we directly convert to array indices
        if (atom > 0)
          atom1.push_back(atom - 1);
      }
      block.get_next_parameter("TYPE", type1, "", "-2,-1,0,1,2,3,4,5,6,7");

      for (unsigned int i = 0; i < 4; i++)
      {
        unsigned int atom;
        std::string str_i = io::to_string(i);
        block.get_next_parameter("ATOM[" + str_i + "]", atom, ">=0", "");
        if (atom > topo.num_atoms())
        {
          std::ostringstream msg;
          msg << blockname << " block: atom number out of range: " << atom << ", last atom is " << topo.num_atoms();
          io::messages.add(msg.str(), "In_TFRDC", io::message::error);
        }

        // -1 because we directly convert to array indices
        if (atom > 0)
          atom2.push_back(atom - 1);
      }
      block.get_next_parameter("TYPE", type2, "", "-2,-1,0,1,2,3,4,5,6,7");

      block.get_next_parameter("R0", r0, ">0", "");
      block.get_next_parameter("GYR1", gyr1, "", "");
      block.get_next_parameter("GYR2", gyr2, "", "");
      block.get_next_parameter("D0", D0, "", "");
      block.get_next_parameter("DD0", DD0, ">=0", "");
      block.get_next_parameter("W0", weight, ">=0", "");

      if (!block.error())
      {

        D0 *= sim.param().tfrdc.factorFreq;
        DD0 *= sim.param().tfrdc.factorFreq;
        gyr1 *= sim.param().tfrdc.factorGyr;
        gyr2 *= sim.param().tfrdc.factorGyr;

        // check for sensible choice of RDC
        if (abs(-(math::eps0_i * math::h_bar * gyr1 * gyr2) / (pow(math::spd_l, 2) * 4.0 * pow(math::Pi, 2)) * 1000) < abs(D0))
        {
          io::messages.add("The chosen RDC is larger in magnitude than RDC_max.  This is probably a mistake and may result in strange behaviour.",
                           "In_Tfrdcresspec", io::message::warning);
        }
        if(fabs(-( math::eps0_i * math::h_bar * gyr1 * gyr2 )/( pow(math::spd_l,2) * 4.0 * pow(math::Pi,2) * pow(r0,3.0))) < fabs(D0)){
          io::messages.add("The chosen RDC is larger in magnitude than RDC_max. This is probably a mistake and may result in strange behaviour.",
           "In_Tfrdcresspec", io::message::warning);
        }

        util::virtual_type t1 = util::virtual_type(type1);
        util::virtual_type t2 = util::virtual_type(type2);

        util::Virtual_Atom v1(t1, atom1, dish, disc);
        util::Virtual_Atom v2(t2, atom2, dish, disc);

        topo.tf_rdc_restraints().push_back
            (topology::tf_rdc_restraint_struct(v1, v2, r0, gyr1, gyr2, D0, DD0, weight));

        DEBUG(10, setw(6) << atom1[0]
                          << setw(6) << atom1[1]
                          << setw(6) << atom1[2]
                          << setw(6) << atom1[3]
                          << setw(6) << type1
                          << setw(6) << atom2[0]
                          << setw(6) << atom2[1]
                          << setw(6) << atom2[2]
                          << setw(6) << atom2[3]
                          << setw(6) << type2
                          << setw(8) << r0
                          << setprecision(4) << setw(12) << gyr1
                          << setw(12) << gyr2
                          << setprecision(14) << setw(19) << D0
                          << setprecision(14) << setw(19) << DD0
                          << setprecision(2) << setw(8) << weight);
      }
    } // for restraint-lines
    block.get_final_messages();

  } // if block content
} // TFRDCRESSPEC



/**
 * @section tfrdcmolaxis TFRDCMOLAXIS block
 * The TFRDCMOLAXIS block is read from the TFRDC restraint specification file.
 * It defines a molecular axis, the angle of which with the magn. field vector is calculated to write out its distribution at the end of a run.
 * @snippet snippets/snippets.cc TFRDCMOLAXIS
 */
void io::In_Tfrdcresspec::read_TFRDCMOLAXIS(topology::Topology &topo,
                                 simulation::Simulation &sim,
                                 std::ostream &os)
{ // TFRDCMOLAXIS
    DEBUG(10, "TFRDCMOLAXIS block");
  std::stringstream exampleblock;
  // lines starting with 'exampleblock<<"' and ending with '\n";' (spaces don't matter)
  // will be used to generate snippets that can be included in the doxygen doc;
  // the first line is the tag
  exampleblock << "TFRDCMOLAXIS\n";
  exampleblock << "# i,j,k,l  atoms comprising the virtual atom (put 0 if less than four atoms in use)\n";
  exampleblock << "# type   virtual atom type \n";
  exampleblock << "# i   j   k   l   type  i    j    k   l  type\n";
  exampleblock << "  16  0   0   0    0    17   0    0   0  0\n";
  exampleblock << "END\n";
  std::string blockname = "TFRDCMOLAXIS";
  Block block(blockname, exampleblock.str());

  if (block.read_buffer(m_block[blockname], true) == 0)
  {
    block_read.insert(blockname);

    DEBUG(10, "reading in TFRDCMOLAXIS data");

    int type1, type2;
    std::vector<int> atom1, atom2;

    for (unsigned int i = 0; i < io::In_Tfrdcresspec::MAX_ATOMS; i++)
    {
      unsigned int atom;
      std::string str_i = io::to_string(i);
      std::string condition = i == 0 ? ">0" : ">=0";
      block.get_next_parameter("ATOM[" + str_i + "]", atom, condition, "");
      if (atom > topo.num_atoms())
      {
        std::ostringstream msg;
        msg << blockname << " block: atom number out of range: " << atom << ", last atom is " << topo.num_atoms();
        io::messages.add(msg.str(), "In_Tfrdcresspec", io::message::error);
      }

      // -1 because we directly convert to array indices
      if (atom > 0)
        atom1.push_back(atom - 1);
    }
    block.get_next_parameter("type1", type1, "", "");

    for (unsigned int i = 0; i < io::In_Tfrdcresspec::MAX_ATOMS; i++)
    {
      unsigned int atom;
      std::string str_i = io::to_string(i);
      std::string condition = i == 0 ? ">0" : ">=0";
      block.get_next_parameter("ATOM[" + str_i + "]", atom, condition, "");
      if (atom > topo.num_atoms())
      {
        std::ostringstream msg;
        msg << blockname << " block: atom number out of range: " << atom << ", last atom is " << topo.num_atoms();
        io::messages.add(msg.str(), "In_Tfrdcresspec", io::message::error);
      }

      // -1 because we directly convert to array indices
      if (atom > 0)
        atom2.push_back(atom - 1);
    }
    block.get_next_parameter("type2", type2, "", "");

    util::virtual_type t1 = util::virtual_type(type1);
    util::virtual_type t2 = util::virtual_type(type2);

    util::Virtual_Atom v1(t1, atom1, dish, disc);
    util::Virtual_Atom v2(t2, atom2, dish, disc);


    topo.tf_rdc_molaxis().push_back(v1);
    topo.tf_rdc_molaxis().push_back(v2);
    block.get_final_messages();
  } // if block empty or not there
} // TFRDCMOLAXIS