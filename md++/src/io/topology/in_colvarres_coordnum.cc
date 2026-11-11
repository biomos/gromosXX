/*
 * This file is part of GROMOS.
 *
 * Copyright (c) 2011, 2012, 2016, 2018, 2021, 2023 Biomos b.v.
 * See <https://www.gromos.net> for details.
 */

/**
 * @file in_colvarres_coordnum.cc
 * implements methods of In_Topology for COLVAR geometry specification files.
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

#include "in_colvarres_helpers.h"

using namespace io::colvar_input;

void io::In_Colvarres::read_coordnumresspec(topology::Topology &topo, simulation::Simulation &sim, std::ostream &os)
{
  std::vector<std::string> buffer;
  std::set<std::string> &block_read = m_colvar_blocks_read;

  // -------------------------------------------------------------------------
  // COORDNUMRESSPEC
  // -------------------------------------------------------------------------
  {
    const std::string blockname = "COORDNUMRESSPEC";
    DEBUG(10, "COORDNUMRESSPEC block");
    buffer = m_block[blockname];

    if (buffer.size() <= 2) {
      if (count_colvar_specs(sim, "COORDNUM") > 0) {
        io::messages.add("RESTRAINTS contains COORDNUM entries but no non-empty "
                         "COORDNUMRESSPEC block was found",
                         "In_Colvarres", io::message::error);
      }
    } else {
      block_read.insert(blockname);

      size_t coordnum_index = 0;
      unsigned int geom_index = 0;
      double rcut = 0.0;
      int nn = 0, mm = 0;

      DEBUG(10, "reading in COORDNUM COLVAR geometry");

      std::vector<std::string>::const_iterator data_it = buffer.begin() + 1;
      const std::vector<std::string>::const_iterator data_to = buffer.end() - 1;
      for (; data_it != data_to && !line_is_data(*data_it); ++data_it) {}

      if (data_it != data_to) {
        _lineStream.clear();
        _lineStream.str(*data_it);
        std::vector<std::string> tokens = tokens_from_line(*data_it);
        if (tokens.size() == 4) {
          _lineStream >> geom_index;
        }
        _lineStream >> rcut >> nn >> mm;
      } else {
        io::messages.add("COORDNUMRESSPEC block: missing RCUT N M data line",
                         "In_Colvarres", io::message::error);
      }

      if (_lineStream.fail()) {
        io::messages.add("bad first data line in COORDNUMRESSPEC block: "
                         "expected RCUT N M",
                         "In_Colvarres", io::message::error);
      }

      if (!quiet) {
        os << std::setw(10) << "RCUT"
           << std::setw(10) << "N"
           << std::setw(10) << "M"
           << "\n"
           << std::setw(10) << rcut
           << std::setw(10) << nn
           << std::setw(10) << mm
           << "\n";
        os << std::setw(8) << "TYPE"
           << std::setw(8) << "GROUP"
           << std::setw(10) << "ATOM"
           << "\n";
      }

      std::vector<std::vector<int> > atomgroup(2);
      std::vector<int> vatypes(2, 0);
      std::vector<bool> seen_type(2, false);

      for (std::vector<std::string>::const_iterator it =
             data_it == data_to ? data_to : data_it + 1,
           to = data_to; it != to; ++it) {
        std::string line(*it);
        if (!line_is_data(line)) continue;

        _lineStream.clear();
        _lineStream.str(line);

        int type = 0;
        int group = 0;
        unsigned int atom_no = 0;

        _lineStream >> type >> group >> atom_no;

        if (_lineStream.fail()) {
          io::messages.add("bad line in COORDNUMRESSPEC block: expected TYPE GROUP ATOM",
                           "In_Colvarres", io::message::error);
          continue;
        }
        if (group != 1 && group != 2) {
          io::messages.add("COORDNUMRESSPEC block: GROUP must be 1 or 2",
                           "In_Colvarres", io::message::error);
          continue;
        }
        if (atom_no == 0 || atom_no > topo.num_atoms()) {
          std::ostringstream msg;
          msg << "COORDNUMRESSPEC block: atom number out of range: "
              << atom_no << ", last atom is " << topo.num_atoms();
          io::messages.add(msg.str(), "In_Colvarres", io::message::error);
          continue;
        }

        const int ag = group - 1;
        if (!seen_type[ag]) {
          vatypes[ag] = type;
          seen_type[ag] = true;
        } else if (vatypes[ag] != type) {
          io::messages.add("COORDNUMRESSPEC block: all atoms within one group "
                           "must use the same virtual atom type",
                           "In_Colvarres", io::message::error);
        }

        atomgroup[ag].push_back(atom_no - 1);

        if (!quiet) {
          os << std::setw(8) << type
             << std::setw(8) << group
             << std::setw(10) << atom_no
             << "\n";
        }
      }

      std::vector<std::vector<util::Virtual_Atom> > atoms;
      make_coordnum_virtual_atoms(atomgroup, vatypes, atoms);

      const simulation::Parameter::colvar_bias_spec *spec =
        geom_index > 0
          ? indexed_colvar_spec(sim.param().colvarres.bias_specs,
                                "COORDNUM", geom_index)
          : nth_colvar_spec(sim, "COORDNUM", coordnum_index);
      if (spec == NULL) {
        if (geom_index > 0) {
          std::ostringstream msg;
          msg << "COORDNUMRESSPEC index " << geom_index
              << " has no matching RESTRAINTS COORDNUM entry; geometry is ignored.";
          io::messages.add(msg.str(), "In_Colvarres", io::message::warning);
        }
        else {
          io::messages.add("COORDNUMRESSPEC contains more geometries than RESTRAINTS COORDNUM entries; extra geometries are ignored.",
                           "In_Colvarres", io::message::warning);
        }
      }
      else {
      topo.coordnum_restraint().push_back
        (topology::coordnum_restraint_struct(true, atoms[0], atoms[1],
                                             0.0, 1.0, rcut, nn, mm));
      }
      ++coordnum_index;

      if (!quiet) os << "END\n";
    }
  }
}

void io::In_Colvarres::read_pertcoordnumresspec(topology::Topology &topo, simulation::Simulation &sim, std::ostream &os)
{
  std::vector<std::string> buffer;
  std::set<std::string> &block_read = m_colvar_blocks_read;

  // -------------------------------------------------------------------------
  // PERTCOORDNUMRESSPEC
  // -------------------------------------------------------------------------
  {
    const std::string blockname = "PERTCOORDNUMRESSPEC";
    DEBUG(10, "PERTCOORDNUMRESSPEC block");
    buffer = m_block[blockname];

    if (buffer.size() <= 2) {
      if (count_colvar_specs(sim.param().colvarres.pert_bias_specs, "COORDNUM") > 0) {
        io::messages.add("PERTRESTRAINTS contains COORDNUM entries but no non-empty "
                         "PERTCOORDNUMRESSPEC block was found",
                         "In_Colvarres", io::message::error);
      }
    } else {
      block_read.insert(blockname);

      size_t coordnum_index = 0;
      unsigned int geom_index = 0;
      double rcut = 0.0;
      int nn = 0, mm = 0;
      std::vector<std::string>::const_iterator data_it = buffer.begin() + 1;
      const std::vector<std::string>::const_iterator data_to = buffer.end() - 1;
      for (; data_it != data_to && !line_is_data(*data_it); ++data_it) {}

      if (data_it != data_to) {
        _lineStream.clear();
        _lineStream.str(*data_it);
        std::vector<std::string> tokens = tokens_from_line(*data_it);
        if (tokens.size() == 4) {
          _lineStream >> geom_index;
        }
        _lineStream >> rcut >> nn >> mm;
      }
      if (_lineStream.fail()) {
        io::messages.add("bad first data line in PERTCOORDNUMRESSPEC block: expected RCUT N M",
                         "In_Colvarres", io::message::error);
      }

      std::vector<std::vector<int> > atomgroup(2);
      std::vector<int> vatypes(2, 0);
      std::vector<bool> seen_type(2, false);

      for (std::vector<std::string>::const_iterator it =
             data_it == data_to ? data_to : data_it + 1,
           to = data_to; it != to; ++it) {
        if (!line_is_data(*it)) continue;

        int type = 0, group = 0;
        unsigned int atom_no = 0;
        _lineStream.clear();
        _lineStream.str(*it);
        _lineStream >> type >> group >> atom_no;

        if (_lineStream.fail() || (group != 1 && group != 2)
            || atom_no == 0 || atom_no > topo.num_atoms()) {
          io::messages.add("bad line in PERTCOORDNUMRESSPEC block: expected valid TYPE GROUP ATOM",
                           "In_Colvarres", io::message::error);
          continue;
        }

        const int ag = group - 1;
        if (!seen_type[ag]) {
          vatypes[ag] = type;
          seen_type[ag] = true;
        }
        atomgroup[ag].push_back(atom_no - 1);
      }

      std::vector<std::vector<util::Virtual_Atom> > atoms;
      make_coordnum_virtual_atoms(atomgroup, vatypes, atoms);

      const simulation::Parameter::colvar_bias_spec *spec =
        geom_index > 0
          ? indexed_colvar_spec(sim.param().colvarres.pert_bias_specs,
                                "COORDNUM", geom_index)
          : nth_colvar_spec(sim.param().colvarres.pert_bias_specs,
                            "COORDNUM", 2 * coordnum_index);
      if (spec == NULL) {
        if (geom_index > 0) {
          std::ostringstream msg;
          msg << "PERTCOORDNUMRESSPEC index " << geom_index
              << " has no matching PERTRESTRAINTS COORDNUM entry; geometry is ignored.";
          io::messages.add(msg.str(), "In_Colvarres", io::message::warning);
        }
        else {
          io::messages.add("PERTCOORDNUMRESSPEC contains more geometries than PERTRESTRAINTS COORDNUM entries; extra geometries are ignored.",
                           "In_Colvarres", io::message::warning);
        }
      }
      else {
      topo.coordnum_restraint().push_back
        (topology::coordnum_restraint_struct(true, atoms[0], atoms[1],
                                             0.0, 1.0, rcut, nn, mm));
      }
      ++coordnum_index;
    }
  }
}

namespace io {
namespace colvar_input {

void make_coordnum_virtual_atoms(
    const std::vector<std::vector<int> > &atomgroup,
    const std::vector<int> &vatypes,
    std::vector<std::vector<util::Virtual_Atom> > &atoms) {

  atoms.clear();
  atoms.resize(2);

  for (int ag = 0; ag < 2; ++ag) {
    if (atomgroup[ag].empty()) {
      io::messages.add("COORDNUMRESSPEC block: no atoms in atom group",
                       "In_Colvarres", io::message::error);
      continue;
    }

    if (vatypes[ag] == 0) {
      // Real atoms: create one virtual atom per atom.
      for (size_t i = 0; i < atomgroup[ag].size(); ++i) {
        std::vector<int> atomvec(1, atomgroup[ag][i]);
        util::Virtual_Atom v(util::virtual_type(vatypes[ag]), atomvec,
                             0.1, 0.153);
        atoms[ag].push_back(v);
      }
    } else {
      // Non-zero virtual atom type: create one virtual atom from the whole group.
      std::vector<int> atomvec = atomgroup[ag];
      util::Virtual_Atom v(util::virtual_type(vatypes[ag]), atomvec,
                           0.1, 0.153);
      atoms[ag].push_back(v);
    }
  }
}

} // namespace colvar_input
} // namespace io
