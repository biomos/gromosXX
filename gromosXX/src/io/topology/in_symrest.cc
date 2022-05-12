/**
 * @file in_symrest.cc
 * implements methods of In_Symrest
 */

#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../interaction/interaction_types.h"
#include "../../configuration/configuration.h"

#include "../../io/instream.h"
#include "../../io/blockinput.h"
#include "../../io/configuration/in_configuration.h"
#include "../../io/configuration/out_configuration.h"

#include "in_symrest.h"

#undef MODULE
#undef SUBMODULE
#define MODULE io
#define SUBMODULE topology

/**
 * @section transform TRANSFORM
@verbatim
TRANSFORM
# Number of symmetry transformations 
# SYMNUMSYM
     1
# Transformations
# M                                   V
  1.0         0.0         0.0         0.0  
  0.0         1.0         0.0         0.0  
  0.0         0.0         1.0         0.0  
END
@endverbatim
 *
 * @section symresspec SYMRESSPEC block
 * The SYMRESSPEC block specifies parameters of the method. On an additional
 * line data about asymmetric units have to be given. The rest is
 * just the specification of the atoms.
 *
 * @verbatim
SYMRESSPEC
# ASUDEF     : pointer to first atom in every ASU
# SYMATOMS   : the atoms to restrain
#
# ASUDEF[1..SYMNUMSYM]
  1  101  201  301
# SYMATOMS
    2 HEXA  CH23       9
    2 HEXA  CH24      10
    2 HEXA  CH25      11
    2 HEXA  CH36      12
END
@endverbatim
 */
void
io::In_Symrest::read(topology::Topology& topo,
        simulation::Simulation & sim,
        std::ostream & os) {

  DEBUG(7, "reading in a xray restraints file");

  if (!quiet)
    os << "SYMMETRY RESTRAINTS\n";

  std::vector<std::string> buffer;
  { // TRANSFORM
    buffer = m_block["TRANSFORM"];
    DEBUG(10, "TRANSFORM block : " << buffer.size());

    if (buffer.empty()) {
      io::messages.add("TRANSFORM block missing",
              "In_Symrest", io::message::error);
      return;
    }

    if (buffer.size() < 3) {
      io::messages.add("TRANSFORM block: Not enough lines given",
              "In_Symrest", io::message::error);
      return;
    }

    unsigned int num_sym = 0;
    _lineStream.clear();
    _lineStream.str(buffer[1]);
    _lineStream >> num_sym;
    if (_lineStream.fail()) {
      io::messages.add("TRANSFORM block: Cannot read number of symmetry transformations",
              "In_Symrest", io::message::error);
      return;
    }

    math::Matrix rot;
    math::Vec v;

    std::vector<std::string>::const_iterator it = buffer.begin() + 2,
            to = buffer.end() - 1;

    for (unsigned int i = 0; i < num_sym; i++) {
      for (unsigned int j = 0; j < 3; j++, ++it) {
        if (it == to) {
          io::messages.add("TRANSFORM block: Not enough lines given",
                  "In_Symrest", io::message::error);
          return;
        }

        _lineStream.clear();
        _lineStream.str(*it);
        for (int k = 0; k < 3; ++k) {
          if (!(_lineStream >> rot(j, k))) {
            io::messages.add("TRANSFORM block: Cannot read matrix.",
                    "In_Symrest", io::message::error);
            return;
          }
        }
        if (!(_lineStream >> v[j])) {
          io::messages.add("TRANSFORM block: Cannot read vector.",
                  "In_Symrest", io::message::error);
          return;
        }
      }
      sim.param().symrest.symmetry_operations.push_back(std::pair<math::Matrix, math::Vec > (rot, v));
    }
  }

  { // SYMRESSPEC
    buffer = m_block["SYMRESSPEC"];
    DEBUG(10, "SYMRESSPEC block : " << buffer.size());

    if (buffer.size() < 3) {
      io::messages.add("SYMRESSPEC block: not enough lines given",
              "In_Symrest", io::message::error);
      return;
    }

    _lineStream.clear();
    _lineStream.str(buffer[1]);

    topo.sym_asu().clear();
    for (unsigned int i = 0; i < sim.param().symrest.symmetry_operations.size(); ++i) {
      int atom_pointer = 0;
      _lineStream >> atom_pointer;
      if (_lineStream.fail()) {
        std::ostringstream msg;
        msg << "SYMRESSPEC block: Cannot read ASU pointer " << (i + 1) << ".";
        io::messages.add("In_Symrest", msg.str(), io::message::error);
        return;
      }
      --atom_pointer;
      if (atom_pointer < 0 || atom_pointer >= int(topo.num_atoms())) {
        std::ostringstream msg;
        msg << "SYMRESSPEC block: ASU pointer " << (i + 1) << " is out of range.";
        io::messages.add("In_Symrest", msg.str(), io::message::error);
        return;
      }
      DEBUG(9, "\t\tASU: " << atom_pointer);
      topo.sym_asu().push_back(atom_pointer);
    }

    topo.sym_restraints().clear();
    std::vector<std::string>::const_iterator it = buffer.begin() + 2,
            to = buffer.end() - 1;
    for (unsigned int line_nr = 6; it != to; ++it, ++line_nr) {
      std::string line(*it);
      if (line.length() < 17) {
        std::ostringstream msg;
        msg << "SYMRESSPEC block: Line " << line_nr << " is too short.";
        io::messages.add("In_Symrest", msg.str(), io::message::error);
        return;
      }

      // the first 17 chars are ignored
      line.erase(line.begin(), line.begin() + 17);

      _lineStream.clear();
      _lineStream.str(line);

      int atom = 0;
      _lineStream >> atom;

      DEBUG(11, "\t" << atom);

      if (_lineStream.fail()) {
        std::ostringstream msg;
        msg << "SYMRESSPEC block: Line " << line_nr << ": Cannot read atom.";
        io::messages.add("In_Symrest", msg.str(), io::message::error);
        return;
      }

      --atom;
      if (atom < 0 || atom >= int(topo.num_atoms())) {
        std::ostringstream msg;
        msg << "SYMRESSPEC block: Line " << line_nr << ": Atom out of range.";
        io::messages.add("In_Symrest", msg.str(), io::message::error);
        return;
      }
      if (atom < int(topo.sym_asu()[0])) {
        std::ostringstream msg;
        msg << "SYMRESSPEC block: Line " << line_nr << ": Atom not in first ASU.";
        io::messages.add("In_Symrest", msg.str(), io::message::error);
        return;
      }
      const unsigned int atom_p = atom - topo.sym_asu()[0];
      for (unsigned int i = 1; i < topo.sym_asu().size(); ++i) {
        const unsigned int atom_img = topo.sym_asu()[i] + atom_p;
        if (atom_img >= topo.num_atoms()) {
          std::ostringstream msg;
          msg << "SYMRESSPEC block: Line " << line_nr << ": The image nr. "
                  << i << " of atom " << (atom + 1) << " is out of range.";
          io::messages.add("In_Symrest", msg.str(), io::message::error);
          return;
        }
      }
      topo.sym_restraints().push_back(atom);
    } // for atoms
  } // SYMRESSPEC


  if (!quiet) {

    switch (sim.param().symrest.symrest) {
      case simulation::xray_symrest_off:
        os << "\tsymmetry restraints OFF\n";
        break;
      case simulation::xray_symrest_ind:
        os << "\tsymmetry restraints on individual atom positions\n";
        break;
      case simulation::xray_symrest_constr:
        os << "\tsymmetry constraints on individual atom positions\n";
        break;
    }
    os << "END\n";
  }
}


