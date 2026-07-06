/*
 * This file is part of GROMOS.
 *
 * Copyright (c) 2011, 2012, 2016, 2018, 2021, 2023 Biomos b.v.
 * See <https://www.gromos.net> for details.
 */

/**
 * @file in_colvarres.cc
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

static std::set<std::string> block_read;

namespace {

bool line_is_data(const std::string &line) {
  const std::string::size_type pos = line.find_first_not_of(" \t\r\n");
  return pos != std::string::npos && line[pos] != '#';
}

std::string colvar_type_from_int(int type) {
  switch (type) {
    case 0: return "DISTANCE";
    case 1: return "ANGLE";
    case 2: return "DIHEDRAL";
    case 3: return "COORDNUM";
    default: return "";
  }
}

bool decode_average_mode(int raw,
                         int &averaging,
                         int &force_scale,
                         const std::string &type) {
  averaging = raw % 10;
  force_scale = raw / 10;

  if (raw < 0 || averaging < 0 || averaging > 2
      || force_scale < 0 || force_scale > 2) {
    io::messages.add("RESTRAINTS block: AVG must be one of 0, 1, 2, 11, 12, or 22.",
                     "In_Colvarres", io::message::error);
    return false;
  }

  if (force_scale > 0 && averaging == 0) {
    io::messages.add("RESTRAINTS block: force scaling in AVG requires time averaging.",
                     "In_Colvarres", io::message::error);
    return false;
  }

  if (force_scale == 2 && averaging != 2) {
    io::messages.add("RESTRAINTS block: chain-rule force scaling is only meaningful for AVG 2.",
                     "In_Colvarres", io::message::error);
    return false;
  }

  if (force_scale > 0 && type != "DISTANCE") {
    io::messages.add("RESTRAINTS block: force-scaled averaging is currently only supported for DISTANCE colvars.",
                     "In_Colvarres", io::message::error);
    return false;
  }

  return true;
}

const simulation::Parameter::colvar_bias_spec *nth_colvar_spec(
    const simulation::Simulation &sim,
    const std::string &type,
    size_t index) {
  size_t n = 0;
  for (std::vector<simulation::Parameter::colvar_bias_spec>::const_iterator
       it = sim.param().colvarres.bias_specs.begin(),
       to = sim.param().colvarres.bias_specs.end(); it != to; ++it) {
    if (it->type == type) {
      if (n == index) return &(*it);
      ++n;
    }
  }
  return NULL;
}

unsigned int count_colvar_specs(const simulation::Simulation &sim,
                                const std::string &type) {
  unsigned int n = 0;
  for (std::vector<simulation::Parameter::colvar_bias_spec>::const_iterator
       it = sim.param().colvarres.bias_specs.begin(),
       to = sim.param().colvarres.bias_specs.end(); it != to; ++it) {
    if (it->type == type) ++n;
  }
  return n;
}

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

} // anonymous namespace

/**
 * @section coordnumresspec COORDNUMRESSPEC block
 * Geometry and switching-function parameters for one COORDNUM collective
 * variable. The target, force constant, mode, averaging, virial, and
 * linear-tail options are read from the RESTRAINTS block in this file.
 * AVG encodes both the averaging transform and, optionally, the force scaling:
 * 0 no averaging, 1 exponential average, 2 inverse-cubic distance average,
 * 11/12 with relaxation force scaling, 22 with chain-rule force scaling.
 *
 * @verbatim
COORDNUMRESSPEC
# RCUT    N      M
  0.3    20     40
# type  GROUP[1,2]   ATOM[i]
     0          1      22
     0          2      96
     0          2      99
END
@endverbatim
 */
void io::In_Colvarres::read(topology::Topology& topo,
                            simulation::Simulation & sim,
                            std::ostream & os) {

  DEBUG(7, "reading in a colvar restraints file");
  block_read.clear();
  sim.param().colvarres.bias_specs.clear();

  if (!quiet)
    os << "COLVAR RESTRAINTS\n";

  std::vector<std::string> buffer;

  // -------------------------------------------------------------------------
  // RESTRAINTS
  // -------------------------------------------------------------------------
  {
    const std::string blockname = "RESTRAINTS";
    DEBUG(10, "RESTRAINTS block");
    buffer = m_block[blockname];

    if (buffer.size() <= 2) {
      io::messages.add("COLVARRES is on but no non-empty RESTRAINTS block was found",
                       "In_Colvarres", io::message::error);
    } else {
      block_read.insert(blockname);

      std::vector<std::string>::const_iterator it = buffer.begin() + 1;
      const std::vector<std::string>::const_iterator to = buffer.end() - 1;

      for (; it != to && !line_is_data(*it); ++it) {}

      unsigned int nres = 0;
      if (it != to) {
        _lineStream.clear();
        _lineStream.str(*it);
        _lineStream >> nres;
      }
      if (_lineStream.fail() || nres == 0) {
        io::messages.add("bad NRES line in RESTRAINTS block",
                         "In_Colvarres", io::message::error);
      }

      if (!quiet) {
        os << std::setw(8) << "TYPE"
           << std::setw(12) << "TARGET"
           << std::setw(12) << "CVK"
           << std::setw(8) << "RAH"
           << std::setw(8) << "AVG"
           << std::setw(10) << "FSCALE"
           << std::setw(10) << "TAU"
           << std::setw(10) << "VIRIAL"
           << std::setw(10) << "CVLIN"
           << "\n";
      }

      unsigned int parsed = 0;
      if (it != to) ++it;
      for (; it != to && parsed < nres; ++it) {
        if (!line_is_data(*it)) continue;

        int type = -1;
        int rah = 0;
        simulation::Parameter::colvar_bias_spec spec;

        _lineStream.clear();
        _lineStream.str(*it);
        _lineStream >> type
                    >> spec.target
                    >> spec.k
                    >> rah
                    >> spec.averaging
                    >> spec.tau
                    >> spec.virial
                    >> spec.linear_tail;

        if (_lineStream.fail()) {
          io::messages.add("bad line in RESTRAINTS block: expected TYPE TARGET CVK RAH AVG TAU VIRIAL CVLIN",
                           "In_Colvarres", io::message::error);
          continue;
        }

        spec.type = colvar_type_from_int(type);
        if (spec.type.empty()) {
          io::messages.add("RESTRAINTS block: TYPE must be 0 DISTANCE, 1 ANGLE, 2 DIHEDRAL, or 3 COORDNUM",
                           "In_Colvarres", io::message::error);
          continue;
        }
        const double input_target = spec.target;
        if (spec.type == "ANGLE" || spec.type == "DIHEDRAL") {
          spec.target *= math::Pi / 180.0;
          spec.linear_tail *= math::Pi / 180.0;
        }

        spec.rah = rah;

        if (spec.k < 0.0) {
          io::messages.add("RESTRAINTS block: CVK must be >= 0.0.",
                           "In_Colvarres", io::message::error);
        }
        decode_average_mode(spec.averaging, spec.averaging,
                            spec.force_scale, spec.type);
        if (spec.tau < 0.0) {
          io::messages.add("RESTRAINTS block: TAU must be >= 0.0.",
                           "In_Colvarres", io::message::error);
        }
        if (spec.virial != 0 && spec.virial != 1) {
          io::messages.add("RESTRAINTS block: VIRIAL must be 0 or 1.",
                           "In_Colvarres", io::message::error);
        }

        sim.param().colvarres.bias_specs.push_back(spec);
        ++parsed;

        if (!quiet) {
          os << std::setw(8) << spec.type
             << std::setw(12) << input_target
             << std::setw(12) << spec.k
             << std::setw(8) << rah
             << std::setw(8) << spec.averaging
             << std::setw(10) << spec.force_scale
             << std::setw(10) << spec.tau
             << std::setw(10) << spec.virial
             << std::setw(10) << spec.linear_tail
             << "\n";
        }
      }

      if (parsed != nres) {
        std::ostringstream msg;
        msg << "RESTRAINTS block: expected " << nres
            << " restraint line(s), parsed " << parsed << ".";
        io::messages.add(msg.str(), "In_Colvarres", io::message::error);
      }
    }
  }

  // -------------------------------------------------------------------------
  // DISTANCERESSPEC
  // -------------------------------------------------------------------------
  {
    const std::string blockname = "DISTANCERESSPEC";
    DEBUG(10, "DISTANCERESSPEC block");
    buffer = m_block[blockname];

    if (buffer.size() <= 2) {
      if (count_colvar_specs(sim, "DISTANCE") > 0) {
        io::messages.add("RESTRAINTS contains DISTANCE entries but no non-empty "
                         "DISTANCERESSPEC block was found",
                         "In_Colvarres", io::message::error);
      }
    } else {
      block_read.insert(blockname);

      std::vector<std::string>::const_iterator it = buffer.begin() + 1;
      const std::vector<std::string>::const_iterator to = buffer.end() - 1;

      for (; it != to && !line_is_data(*it); ++it) {}

      double dish = 0.0, disc = 0.0;
      if (it != to) {
        _lineStream.clear();
        _lineStream.str(*it);
        _lineStream >> dish >> disc;
      }
      if (_lineStream.fail()) {
        io::messages.add("bad first data line in DISTANCERESSPEC block: expected DISH DISC",
                         "In_Colvarres", io::message::error);
      }

      if (!quiet) {
        os << std::setw(10) << "DISH"
           << std::setw(10) << "DISC"
           << "\n"
           << std::setw(10) << dish
           << std::setw(10) << disc
           << "\n";
      }

      size_t distance_index = 0;
      if (it != to) ++it;
      for (; it != to; ++it) {
        if (!line_is_data(*it)) continue;

        int type1 = 0, type2 = 0;
        std::vector<int> atom1, atom2;

        _lineStream.clear();
        _lineStream.str(*it);

        for (unsigned int i = 0; i < io::In_Colvarres::MAX_ATOMS; ++i) {
          unsigned int atom = 0;
          _lineStream >> atom;
          if (atom > topo.num_atoms()) {
            std::ostringstream msg;
            msg << blockname << " block: atom number out of range: "
                << atom << ", last atom is " << topo.num_atoms();
            io::messages.add(msg.str(), "In_Colvarres", io::message::error);
          }
          if (atom > 0) atom1.push_back(atom - 1);
        }
        _lineStream >> type1;

        for (unsigned int i = 0; i < io::In_Colvarres::MAX_ATOMS; ++i) {
          unsigned int atom = 0;
          _lineStream >> atom;
          if (atom > topo.num_atoms()) {
            std::ostringstream msg;
            msg << blockname << " block: atom number out of range: "
                << atom << ", last atom is " << topo.num_atoms();
            io::messages.add(msg.str(), "In_Colvarres", io::message::error);
          }
          if (atom > 0) atom2.push_back(atom - 1);
        }
        _lineStream >> type2;

        if (_lineStream.fail()) {
          io::messages.add("bad line in DISTANCERESSPEC block: expected two virtual atom definitions",
                           "In_Colvarres", io::message::error);
          continue;
        }

        const simulation::Parameter::colvar_bias_spec *spec =
          nth_colvar_spec(sim, "DISTANCE", distance_index);
        if (spec == NULL) {
          io::messages.add("DISTANCERESSPEC contains more geometries than RESTRAINTS DISTANCE entries; extra geometries are ignored.",
                           "In_Colvarres", io::message::warning);
          ++distance_index;
          continue;
        }

        util::Virtual_Atom v1(util::virtual_type(type1), atom1, dish, disc);
        util::Virtual_Atom v2(util::virtual_type(type2), atom2, dish, disc);

        topo.distance_restraints().push_back
          (topology::distance_restraint_struct(v1, v2,
                                               spec->target, 1.0,
                                               spec->rah));

        ++distance_index;
      }
    }
  }

  // -------------------------------------------------------------------------
  // ANGRESSPEC
  // -------------------------------------------------------------------------
  {
    const std::string blockname = "ANGRESSPEC";
    DEBUG(10, "ANGRESSPEC block");
    buffer = m_block[blockname];

    if (buffer.size() <= 2) {
      if (count_colvar_specs(sim, "ANGLE") > 0) {
        io::messages.add("RESTRAINTS contains ANGLE entries but no non-empty "
                         "ANGRESSPEC block was found",
                         "In_Colvarres", io::message::error);
      }
    } else {
      block_read.insert(blockname);

      if (!quiet) {
        os << std::setw(10) << "IPLR"
           << std::setw(10) << "JPLR"
           << std::setw(10) << "KPLR"
           << "\n";
      }

      size_t angle_index = 0;
      for (std::vector<std::string>::const_iterator it = buffer.begin() + 1,
           to = buffer.end() - 1; it != to; ++it) {
        if (!line_is_data(*it)) continue;

        unsigned int i = 0, j = 0, k = 0;
        _lineStream.clear();
        _lineStream.str(*it);
        _lineStream >> i >> j >> k;

        if (_lineStream.fail()) {
          io::messages.add("bad line in ANGRESSPEC block: expected IPLR JPLR KPLR",
                           "In_Colvarres", io::message::error);
          continue;
        }
        if (i == 0 || j == 0 || k == 0
            || i > topo.num_atoms() || j > topo.num_atoms() || k > topo.num_atoms()) {
          io::messages.add("ANGRESSPEC block: atom number out of range",
                           "In_Colvarres", io::message::error);
          continue;
        }

        const simulation::Parameter::colvar_bias_spec *spec =
          nth_colvar_spec(sim, "ANGLE", angle_index);
        if (spec == NULL) {
          io::messages.add("ANGRESSPEC contains more geometries than RESTRAINTS ANGLE entries; extra geometries are ignored.",
                           "In_Colvarres", io::message::warning);
          ++angle_index;
          continue;
        }

        topo.angle_restraints().push_back
          (topology::angle_restraint_struct(i - 1, j - 1, k - 1,
                                            spec->target, 1.0));
        ++angle_index;

        if (!quiet) {
          os << std::setw(10) << i
             << std::setw(10) << j
             << std::setw(10) << k
             << "\n";
        }
      }
    }
  }

  // -------------------------------------------------------------------------
  // DIHEDRALRESSPEC
  // -------------------------------------------------------------------------
  {
    const std::string blockname = "DIHEDRALRESSPEC";
    DEBUG(10, "DIHEDRALRESSPEC block");
    buffer = m_block[blockname];

    if (buffer.size() <= 2) {
      if (count_colvar_specs(sim, "DIHEDRAL") > 0) {
        io::messages.add("RESTRAINTS contains DIHEDRAL entries but no non-empty "
                         "DIHEDRALRESSPEC block was found",
                         "In_Colvarres", io::message::error);
      }
    } else {
      block_read.insert(blockname);

      if (!quiet) {
        os << std::setw(10) << "IPLR"
           << std::setw(10) << "JPLR"
           << std::setw(10) << "KPLR"
           << std::setw(10) << "LPLR"
           << "\n";
      }

      size_t dihedral_index = 0;
      for (std::vector<std::string>::const_iterator it = buffer.begin() + 1,
           to = buffer.end() - 1; it != to; ++it) {
        if (!line_is_data(*it)) continue;

        unsigned int i = 0, j = 0, k = 0, l = 0;
        _lineStream.clear();
        _lineStream.str(*it);
        _lineStream >> i >> j >> k >> l;

        if (_lineStream.fail()) {
          io::messages.add("bad line in DIHEDRALRESSPEC block: expected IPLR JPLR KPLR LPLR",
                           "In_Colvarres", io::message::error);
          continue;
        }
        if (i == 0 || j == 0 || k == 0 || l == 0
            || i > topo.num_atoms() || j > topo.num_atoms()
            || k > topo.num_atoms() || l > topo.num_atoms()) {
          io::messages.add("DIHEDRALRESSPEC block: atom number out of range",
                           "In_Colvarres", io::message::error);
          continue;
        }

        const simulation::Parameter::colvar_bias_spec *spec =
          nth_colvar_spec(sim, "DIHEDRAL", dihedral_index);
        if (spec == NULL) {
          io::messages.add("DIHEDRALRESSPEC contains more geometries than RESTRAINTS DIHEDRAL entries; extra geometries are ignored.",
                           "In_Colvarres", io::message::warning);
          ++dihedral_index;
          continue;
        }

        topo.dihedral_restraints().push_back
          (topology::dihedral_restraint_struct(i - 1, j - 1, k - 1, l - 1,
                                               math::Pi, spec->target, 1.0));
        ++dihedral_index;

        if (!quiet) {
          os << std::setw(10) << i
             << std::setw(10) << j
             << std::setw(10) << k
             << std::setw(10) << l
             << "\n";
        }
      }
    }
  }

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

      double rcut = 0.0;
      int nn = 0, mm = 0;

      DEBUG(10, "reading in COORDNUM COLVAR geometry");

      std::vector<std::string>::const_iterator data_it = buffer.begin() + 1;
      const std::vector<std::string>::const_iterator data_to = buffer.end() - 1;
      for (; data_it != data_to && !line_is_data(*data_it); ++data_it) {}

      if (data_it != data_to) {
        _lineStream.clear();
        _lineStream.str(*data_it);
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

      topo.coordnum_restraint().push_back
        (topology::coordnum_restraint_struct(true, atoms[0], atoms[1],
                                             0.0, 1.0, rcut, nn, mm));

      if (!quiet) os << "END\n";
    }
  }

  if (m_block["PERTCOORDNUMRESSPEC"].size()) {
    io::messages.add("PERTCOORDNUMRESSPEC is not handled by the new generic "
                     "COLVAR reader yet.",
                     "In_Colvarres", io::message::warning);
  }

  for (std::map<std::string, std::vector<std::string> >::const_iterator
       it = m_block.begin(), to = m_block.end(); it != to; ++it) {
    if (block_read.count(it->first) == 0 && it->second.size()) {
      io::messages.add(" block " + it->first + " not supported!",
                       "In_Colvarres", io::message::warning);
    }
  }

  const unsigned int n_distance_specs = count_colvar_specs(sim, "DISTANCE");
  if (n_distance_specs != topo.distance_restraints().size()) {
    std::ostringstream msg;
    msg << "RESTRAINTS/DISTANCERESSPEC mismatch: " << n_distance_specs
        << " DISTANCE bias specification(s), but "
        << topo.distance_restraints().size()
        << " DISTANCE geometry specification(s).";
    io::messages.add(msg.str(), "In_Colvarres", io::message::warning);
  }

  const unsigned int n_coord_specs = count_colvar_specs(sim, "COORDNUM");
  if (n_coord_specs != topo.coordnum_restraint().size()) {
    std::ostringstream msg;
    msg << "RESTRAINTS/COORDNUMRESSPEC mismatch: " << n_coord_specs
        << " COORDNUM bias specification(s), but "
        << topo.coordnum_restraint().size()
        << " COORDNUM geometry specification(s).";
    io::messages.add(msg.str(), "In_Colvarres", io::message::warning);
  }

  const unsigned int n_angle_specs = count_colvar_specs(sim, "ANGLE");
  if (n_angle_specs != topo.angle_restraints().size()) {
    std::ostringstream msg;
    msg << "RESTRAINTS/ANGRESSPEC mismatch: " << n_angle_specs
        << " ANGLE bias specification(s), but "
        << topo.angle_restraints().size()
        << " ANGLE geometry specification(s).";
    io::messages.add(msg.str(), "In_Colvarres", io::message::warning);
  }

  const unsigned int n_dihedral_specs = count_colvar_specs(sim, "DIHEDRAL");
  if (n_dihedral_specs != topo.dihedral_restraints().size()) {
    std::ostringstream msg;
    msg << "RESTRAINTS/DIHEDRALRESSPEC mismatch: " << n_dihedral_specs
        << " DIHEDRAL bias specification(s), but "
        << topo.dihedral_restraints().size()
        << " DIHEDRAL geometry specification(s).";
    io::messages.add(msg.str(), "In_Colvarres", io::message::warning);
  }
}
