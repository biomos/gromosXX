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

std::vector<std::string> tokens_from_line(const std::string &line) {
  std::vector<std::string> tokens;
  std::istringstream stream(line);
  std::string token;
  while (stream >> token) tokens.push_back(token);
  return tokens;
}

bool decode_method(int raw,
                   int &averaging,
                   int &force_scale,
                   bool &constraining,
                   const std::string &type,
                   const std::string &blockname) {
  constraining = false;

  if (raw == -1) {
    averaging = 0;
    force_scale = 0;
    constraining = true;
    io::messages.add(blockname + " block: METHOD -1 colvar constraints are "
                     "recognised but not implemented yet.",
                     "In_Colvarres", io::message::error);
    return false;
  }

  averaging = raw % 10;
  force_scale = raw / 10;

  if (raw < 0 || averaging < 0 || averaging > 2
      || force_scale < 0 || force_scale > 2) {
    io::messages.add(blockname + " block: METHOD must be one of -1, 0, 1, 2, 11, 12, or 22.",
                     "In_Colvarres", io::message::error);
    return false;
  }

  if (force_scale > 0 && averaging == 0) {
    io::messages.add(blockname + " block: force scaling in METHOD requires time averaging.",
                     "In_Colvarres", io::message::error);
    return false;
  }

  if (force_scale == 2 && averaging != 2) {
    io::messages.add(blockname + " block: chain-rule force scaling is only meaningful for METHOD 22.",
                     "In_Colvarres", io::message::error);
    return false;
  }

  if (force_scale > 0 && type != "DISTANCE") {
    io::messages.add(blockname + " block: force-scaled averaging is currently only supported for DISTANCE colvars.",
                     "In_Colvarres", io::message::error);
    return false;
  }

  return true;
}

const simulation::Parameter::colvar_bias_spec *indexed_colvar_spec(
    const std::vector<simulation::Parameter::colvar_bias_spec> &specs,
    const std::string &type,
    unsigned int index) {
  for (std::vector<simulation::Parameter::colvar_bias_spec>::const_iterator
       it = specs.begin(), to = specs.end(); it != to; ++it) {
    if (it->type == type && it->index == index) return &(*it);
  }
  return NULL;
}

bool colvar_index_exists(
    const std::vector<simulation::Parameter::colvar_bias_spec> &specs,
    unsigned int index) {
  if (index == 0) return false;
  for (std::vector<simulation::Parameter::colvar_bias_spec>::const_iterator
       it = specs.begin(), to = specs.end(); it != to; ++it) {
    if (it->index == index) return true;
  }
  return false;
}

const simulation::Parameter::colvar_bias_spec *nth_colvar_spec(
    const std::vector<simulation::Parameter::colvar_bias_spec> &specs,
    const std::string &type,
    size_t index) {
  size_t n = 0;
  for (std::vector<simulation::Parameter::colvar_bias_spec>::const_iterator
       it = specs.begin(), to = specs.end(); it != to; ++it) {
    if (it->type == type) {
      if (n == index) return &(*it);
      ++n;
    }
  }
  return NULL;
}

const simulation::Parameter::colvar_bias_spec *nth_colvar_spec(
    const simulation::Simulation &sim,
    const std::string &type,
    size_t index) {
  return nth_colvar_spec(sim.param().colvarres.bias_specs, type, index);
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

unsigned int count_colvar_specs(
    const std::vector<simulation::Parameter::colvar_bias_spec> &specs,
    const std::string &type) {
  unsigned int n = 0;
  for (std::vector<simulation::Parameter::colvar_bias_spec>::const_iterator
       it = specs.begin(), to = specs.end(); it != to; ++it) {
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
 * @section colvarresfile COLVAR restraint specification file
 *
 * The RESTRAINTS block stores the bias parameters. Geometry is read from the
 * matching *RESSPEC block. If a non-zero index is present, the geometry row is
 * matched to the RESTRAINTS row with the same index. Old files without indices
 * are still matched by order within each TYPE.
 *
 * TYPE values:
 * 0 DISTANCE, 1 ANGLE, 2 DIHEDRAL, 3 COORDNUM
 *
 * ANGLE and DIHEDRAL TARGET/CVLIN values are read in degrees and converted
 * internally to radians. METHOD encodes restraining/averaging behavior:
 * -1 colvar constraint (recognised but not implemented yet),
 * 0 no averaging, 1 exponential average, 2 inverse-cubic distance average,
 * 11/12 with relaxation force scaling, 22 with chain-rule force scaling.
 *
 * @verbatim
TITLE
collective variable restraints
END

RESTRAINTS
# NRES
4
# INDX TYPE METHOD TARGET CVK RAH TAU VIRIAL CVLIN
  1    0      0    0.24   10000   0   0.0   1   0.0
  2    1      0     104     490   0   0.0   0   0.0
  3    2      0      -5   15000   0   0.0   0   0.0
  4    3      0       1   10000   0   0.0   1   0.0
END

DISTANCERESSPEC
# DISH DISC
  0.1  0.153
# INDX i j k l TYPE   i j k l TYPE
  1    95 0 0 0 0     14 0 0 0 0
END

ANGRESSPEC
# INDX IPLR JPLR KPLR
  2    14   95   23
END

DIHEDRALRESSPEC
# INDX IPLR JPLR KPLR LPLR
  3    23   24   25   95
END

COORDNUMRESSPEC
# INDX RCUT N M
  4    0.3  20 40
# TYPE GROUP ATOM
  0    1     95
  0    2     96
  0    2     99
END
@endverbatim
 *
 * Perturbed restraints use PERTRESTRAINTS plus matching PERT*RESSPEC geometry
 * blocks. Each PERTRESTRAINTS row defines one perturbed term with A and B
 * state bias parameters. Geometry is not perturbed and can be matched by index
 * in the same way as normal restraints.
 *
 * @verbatim
PERTRESTRAINTS
# NRES
4
# INDX TYPE METHOD TARGETA CVKA TARGETB CVKB RAH TAU VIRIAL CVLIN
  5    0      0    0.24  10000   0.30  20000   0   0.0   1   0.0
  6    1      0     104    490    109    700   0   0.0   0   0.0
  7    2      0      -5  15000     10   9000   0   0.0   0   0.0
  8    3      0       1  10000      2  12000   0   0.0   1   0.0
END

PERTDISTANCERESSPEC
# DISH DISC
  0.1  0.153
# INDX i j k l TYPE   i j k l TYPE
  5    95 0 0 0 0     14 0 0 0 0
END

PERTANGRESSPEC
# INDX IPLR JPLR KPLR
  6    14   95   23
END

PERTDIHEDRALRESSPEC
# INDX IPLR JPLR KPLR LPLR
  7    23   24   25   95
END

PERTCOORDNUMRESSPEC
# INDX RCUT N M
  8    0.3  20 40
# TYPE GROUP ATOM
  0    1     95
  0    2     96
  0    2     99
END
@endverbatim
 */
void io::In_Colvarres::read(topology::Topology& topo,
                            simulation::Simulation & sim,
                            std::ostream & os) {

  DEBUG(7, "reading in a colvar restraints file");
  block_read.clear();
  sim.param().colvarres.bias_specs.clear();
  sim.param().colvarres.pert_bias_specs.clear();

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
      const bool has_pert_restraints = m_block["PERTRESTRAINTS"].size() > 2;
      if (has_pert_restraints) {
        // A perturbed-only colvar file is valid.
      }
      else {
      io::messages.add("COLVARRES is on but no non-empty RESTRAINTS block was found",
                       "In_Colvarres", io::message::error);
      }
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
        os << std::setw(8) << "INDX"
           << std::setw(8) << "TYPE"
           << std::setw(8) << "METHOD"
           << std::setw(12) << "TARGET"
           << std::setw(12) << "CVK"
           << std::setw(8) << "RAH"
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
        std::vector<std::string> tokens = tokens_from_line(*it);

        _lineStream.clear();
        _lineStream.str(*it);
        if (tokens.size() == 9) {
          _lineStream >> spec.index
                      >> type
                      >> spec.method
                      >> spec.target
                      >> spec.k
                      >> rah
                      >> spec.tau
                      >> spec.virial
                      >> spec.linear_tail;
        }
        else {
          _lineStream >> type
                      >> spec.target
                      >> spec.k
                      >> rah
                      >> spec.method
                      >> spec.tau
                      >> spec.virial
                      >> spec.linear_tail;
        }

        if (_lineStream.fail()) {
          io::messages.add("bad line in RESTRAINTS block: expected INDX TYPE METHOD TARGET CVK RAH TAU VIRIAL CVLIN",
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
        decode_method(spec.method, spec.averaging, spec.force_scale,
                      spec.constraining, spec.type, blockname);
        if (spec.tau < 0.0) {
          io::messages.add("RESTRAINTS block: TAU must be >= 0.0.",
                           "In_Colvarres", io::message::error);
        }
        if (spec.virial != 0 && spec.virial != 1) {
          io::messages.add("RESTRAINTS block: VIRIAL must be 0 or 1.",
                           "In_Colvarres", io::message::error);
        }
        if (colvar_index_exists(sim.param().colvarres.bias_specs,
                                spec.index)) {
          std::ostringstream msg;
          msg << "RESTRAINTS block: duplicate index " << spec.index << ".";
          io::messages.add(msg.str(), "In_Colvarres", io::message::error);
        }

        sim.param().colvarres.bias_specs.push_back(spec);
        ++parsed;

        if (!quiet) {
          os << std::setw(8) << spec.index
             << std::setw(8) << spec.type
             << std::setw(8) << spec.method
             << std::setw(12) << input_target
             << std::setw(12) << spec.k
             << std::setw(8) << rah
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
  // PERTRESTRAINTS
  // -------------------------------------------------------------------------
  {
    const std::string blockname = "PERTRESTRAINTS";
    DEBUG(10, "PERTRESTRAINTS block");
    buffer = m_block[blockname];

    if (buffer.size() > 2) {
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
        io::messages.add("bad NRES line in PERTRESTRAINTS block",
                         "In_Colvarres", io::message::error);
      }

      if (!quiet) {
        os << "PERTRESTRAINTS\n";
        os << std::setw(8) << "INDX"
           << std::setw(8) << "TYPE"
           << std::setw(8) << "METHOD"
           << std::setw(12) << "TARGETA"
           << std::setw(12) << "CVKA"
           << std::setw(12) << "TARGETB"
           << std::setw(12) << "CVKB"
           << std::setw(8) << "RAH"
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
        simulation::Parameter::colvar_bias_spec specA;
        simulation::Parameter::colvar_bias_spec specB;
        std::vector<std::string> tokens = tokens_from_line(*it);

        _lineStream.clear();
        _lineStream.str(*it);
        if (tokens.size() == 11) {
          _lineStream >> specA.index
                      >> type
                      >> specA.method
                      >> specA.target
                      >> specA.k
                      >> specB.target
                      >> specB.k
                      >> rah
                      >> specA.tau
                      >> specA.virial
                      >> specA.linear_tail;
        }
        else {
          _lineStream >> type
                      >> specA.target
                      >> specA.k
                      >> specB.target
                      >> specB.k
                      >> rah
                      >> specA.method
                      >> specA.tau
                      >> specA.virial
                      >> specA.linear_tail;
        }

        if (_lineStream.fail()) {
          io::messages.add("bad line in PERTRESTRAINTS block: expected INDX TYPE METHOD TARGETA CVKA TARGETB CVKB RAH TAU VIRIAL CVLIN",
                           "In_Colvarres", io::message::error);
          continue;
        }

        specB.index = specA.index;
        specA.type = colvar_type_from_int(type);
        specB.type = specA.type;
        if (specA.type.empty()) {
          io::messages.add("PERTRESTRAINTS block: TYPE must be 0 DISTANCE, 1 ANGLE, 2 DIHEDRAL, or 3 COORDNUM",
                           "In_Colvarres", io::message::error);
          continue;
        }
        const double input_targetA = specA.target;
        const double input_targetB = specB.target;
        if (specA.type == "ANGLE" || specA.type == "DIHEDRAL") {
          specA.target *= math::Pi / 180.0;
          specB.target *= math::Pi / 180.0;
          specA.linear_tail *= math::Pi / 180.0;
        }

        specA.rah = rah;
        specB.rah = rah;
        specB.method = specA.method;
        specB.averaging = specA.averaging;
        specB.tau = specA.tau;
        specB.virial = specA.virial;
        specB.linear_tail = specA.linear_tail;

        if (specA.k < 0.0 || specB.k < 0.0) {
          io::messages.add("PERTRESTRAINTS block: CVK must be >= 0.0.",
                           "In_Colvarres", io::message::error);
        }
        decode_method(specA.method, specA.averaging, specA.force_scale,
                      specA.constraining, specA.type, blockname);
        specB.averaging = specA.averaging;
        specB.force_scale = specA.force_scale;
        specB.constraining = specA.constraining;
        if (specA.tau < 0.0) {
          io::messages.add("PERTRESTRAINTS block: TAU must be >= 0.0.",
                           "In_Colvarres", io::message::error);
        }
        if (specA.virial != 0 && specA.virial != 1) {
          io::messages.add("PERTRESTRAINTS block: VIRIAL must be 0 or 1.",
                           "In_Colvarres", io::message::error);
        }
        if (colvar_index_exists(sim.param().colvarres.pert_bias_specs,
                                specA.index)) {
          std::ostringstream msg;
          msg << "PERTRESTRAINTS block: duplicate index " << specA.index << ".";
          io::messages.add(msg.str(), "In_Colvarres", io::message::error);
        }

        sim.param().colvarres.pert_bias_specs.push_back(specA);
        sim.param().colvarres.pert_bias_specs.push_back(specB);
        sim.param().perturbation.perturbed_par = true;
        ++parsed;

        if (!quiet) {
          os << std::setw(8) << specA.index
             << std::setw(8) << specA.type
             << std::setw(8) << specA.method
             << std::setw(12) << input_targetA
             << std::setw(12) << specA.k
             << std::setw(12) << input_targetB
             << std::setw(12) << specB.k
             << std::setw(8) << rah
             << std::setw(10) << specA.force_scale
             << std::setw(10) << specA.tau
             << std::setw(10) << specA.virial
             << std::setw(10) << specA.linear_tail
             << "\n";
        }
      }

      if (parsed != nres) {
        std::ostringstream msg;
        msg << "PERTRESTRAINTS block: expected " << nres
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
      bool indexed_distance_geometries = false;
      std::vector<std::pair<unsigned int, topology::distance_restraint_struct> >
        indexed_distances;
      if (it != to) ++it;
      for (; it != to; ++it) {
        if (!line_is_data(*it)) continue;

        int type1 = 0, type2 = 0;
        unsigned int geom_index = 0;
        std::vector<int> atom1, atom2;
        std::vector<std::string> tokens = tokens_from_line(*it);

        _lineStream.clear();
        _lineStream.str(*it);
        if (tokens.size() == 11) {
          _lineStream >> geom_index;
        }

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
          geom_index > 0
            ? indexed_colvar_spec(sim.param().colvarres.bias_specs,
                                  "DISTANCE", geom_index)
            : nth_colvar_spec(sim, "DISTANCE", distance_index);
        if (spec == NULL) {
          if (geom_index > 0) {
            std::ostringstream msg;
            msg << "DISTANCERESSPEC index " << geom_index
                << " has no matching RESTRAINTS DISTANCE entry; geometry is ignored.";
            io::messages.add(msg.str(), "In_Colvarres", io::message::warning);
          }
          else {
            io::messages.add("DISTANCERESSPEC contains more geometries than RESTRAINTS DISTANCE entries; extra geometries are ignored.",
                             "In_Colvarres", io::message::warning);
          }
          ++distance_index;
          continue;
        }

        util::Virtual_Atom v1(util::virtual_type(type1), atom1, dish, disc);
        util::Virtual_Atom v2(util::virtual_type(type2), atom2, dish, disc);

        topology::distance_restraint_struct restraint(v1, v2,
                                                      spec->target, 1.0,
                                                      spec->rah);
        if (geom_index > 0) {
          indexed_distance_geometries = true;
          indexed_distances.push_back(std::make_pair(geom_index, restraint));
        }
        else {
          topo.distance_restraints().push_back(restraint);
        }

        ++distance_index;
      }

      if (indexed_distance_geometries) {
        for (std::vector<simulation::Parameter::colvar_bias_spec>::const_iterator
             spec_it = sim.param().colvarres.bias_specs.begin(),
             spec_to = sim.param().colvarres.bias_specs.end();
             spec_it != spec_to; ++spec_it) {
          if (spec_it->type != "DISTANCE" || spec_it->index == 0) continue;
          bool found = false;
          for (std::vector<std::pair<unsigned int, topology::distance_restraint_struct> >::const_iterator
               geom_it = indexed_distances.begin(),
               geom_to = indexed_distances.end(); geom_it != geom_to; ++geom_it) {
            if (geom_it->first != spec_it->index) continue;
            topo.distance_restraints().push_back(geom_it->second);
            found = true;
            break;
          }
          if (!found) {
            std::ostringstream msg;
            msg << "RESTRAINTS DISTANCE index " << spec_it->index
                << " has no matching DISTANCERESSPEC geometry.";
            io::messages.add(msg.str(), "In_Colvarres", io::message::warning);
          }
        }
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
      bool indexed_angle_geometries = false;
      std::vector<std::pair<unsigned int, topology::angle_restraint_struct> >
        indexed_angles;
      for (std::vector<std::string>::const_iterator it = buffer.begin() + 1,
           to = buffer.end() - 1; it != to; ++it) {
        if (!line_is_data(*it)) continue;

        unsigned int geom_index = 0;
        unsigned int i = 0, j = 0, k = 0;
        std::vector<std::string> tokens = tokens_from_line(*it);
        _lineStream.clear();
        _lineStream.str(*it);
        if (tokens.size() == 4) {
          _lineStream >> geom_index;
        }
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
          geom_index > 0
            ? indexed_colvar_spec(sim.param().colvarres.bias_specs,
                                  "ANGLE", geom_index)
            : nth_colvar_spec(sim, "ANGLE", angle_index);
        if (spec == NULL) {
          if (geom_index > 0) {
            std::ostringstream msg;
            msg << "ANGRESSPEC index " << geom_index
                << " has no matching RESTRAINTS ANGLE entry; geometry is ignored.";
            io::messages.add(msg.str(), "In_Colvarres", io::message::warning);
          }
          else {
            io::messages.add("ANGRESSPEC contains more geometries than RESTRAINTS ANGLE entries; extra geometries are ignored.",
                             "In_Colvarres", io::message::warning);
          }
          ++angle_index;
          continue;
        }

        topology::angle_restraint_struct restraint(i - 1, j - 1, k - 1,
                                                   spec->target, 1.0);
        if (geom_index > 0) {
          indexed_angle_geometries = true;
          indexed_angles.push_back(std::make_pair(geom_index, restraint));
        }
        else {
          topo.angle_restraints().push_back(restraint);
        }
        ++angle_index;

        if (!quiet) {
          os << std::setw(10) << i
             << std::setw(10) << j
             << std::setw(10) << k
             << "\n";
        }
      }

      if (indexed_angle_geometries) {
        for (std::vector<simulation::Parameter::colvar_bias_spec>::const_iterator
             spec_it = sim.param().colvarres.bias_specs.begin(),
             spec_to = sim.param().colvarres.bias_specs.end();
             spec_it != spec_to; ++spec_it) {
          if (spec_it->type != "ANGLE" || spec_it->index == 0) continue;
          bool found = false;
          for (std::vector<std::pair<unsigned int, topology::angle_restraint_struct> >::const_iterator
               geom_it = indexed_angles.begin(),
               geom_to = indexed_angles.end(); geom_it != geom_to; ++geom_it) {
            if (geom_it->first != spec_it->index) continue;
            topo.angle_restraints().push_back(geom_it->second);
            found = true;
            break;
          }
          if (!found) {
            std::ostringstream msg;
            msg << "RESTRAINTS ANGLE index " << spec_it->index
                << " has no matching ANGRESSPEC geometry.";
            io::messages.add(msg.str(), "In_Colvarres", io::message::warning);
          }
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
      bool indexed_dihedral_geometries = false;
      std::vector<std::pair<unsigned int, topology::dihedral_restraint_struct> >
        indexed_dihedrals;
      for (std::vector<std::string>::const_iterator it = buffer.begin() + 1,
           to = buffer.end() - 1; it != to; ++it) {
        if (!line_is_data(*it)) continue;

        unsigned int geom_index = 0;
        unsigned int i = 0, j = 0, k = 0, l = 0;
        std::vector<std::string> tokens = tokens_from_line(*it);
        _lineStream.clear();
        _lineStream.str(*it);
        if (tokens.size() == 5) {
          _lineStream >> geom_index;
        }
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
          geom_index > 0
            ? indexed_colvar_spec(sim.param().colvarres.bias_specs,
                                  "DIHEDRAL", geom_index)
            : nth_colvar_spec(sim, "DIHEDRAL", dihedral_index);
        if (spec == NULL) {
          if (geom_index > 0) {
            std::ostringstream msg;
            msg << "DIHEDRALRESSPEC index " << geom_index
                << " has no matching RESTRAINTS DIHEDRAL entry; geometry is ignored.";
            io::messages.add(msg.str(), "In_Colvarres", io::message::warning);
          }
          else {
            io::messages.add("DIHEDRALRESSPEC contains more geometries than RESTRAINTS DIHEDRAL entries; extra geometries are ignored.",
                             "In_Colvarres", io::message::warning);
          }
          ++dihedral_index;
          continue;
        }

        topology::dihedral_restraint_struct restraint(i - 1, j - 1, k - 1, l - 1,
                                                      math::Pi, spec->target, 1.0);
        if (geom_index > 0) {
          indexed_dihedral_geometries = true;
          indexed_dihedrals.push_back(std::make_pair(geom_index, restraint));
        }
        else {
          topo.dihedral_restraints().push_back(restraint);
        }
        ++dihedral_index;

        if (!quiet) {
          os << std::setw(10) << i
             << std::setw(10) << j
             << std::setw(10) << k
             << std::setw(10) << l
             << "\n";
        }
      }

      if (indexed_dihedral_geometries) {
        for (std::vector<simulation::Parameter::colvar_bias_spec>::const_iterator
             spec_it = sim.param().colvarres.bias_specs.begin(),
             spec_to = sim.param().colvarres.bias_specs.end();
             spec_it != spec_to; ++spec_it) {
          if (spec_it->type != "DIHEDRAL" || spec_it->index == 0) continue;
          bool found = false;
          for (std::vector<std::pair<unsigned int, topology::dihedral_restraint_struct> >::const_iterator
               geom_it = indexed_dihedrals.begin(),
               geom_to = indexed_dihedrals.end(); geom_it != geom_to; ++geom_it) {
            if (geom_it->first != spec_it->index) continue;
            topo.dihedral_restraints().push_back(geom_it->second);
            found = true;
            break;
          }
          if (!found) {
            std::ostringstream msg;
            msg << "RESTRAINTS DIHEDRAL index " << spec_it->index
                << " has no matching DIHEDRALRESSPEC geometry.";
            io::messages.add(msg.str(), "In_Colvarres", io::message::warning);
          }
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

  // -------------------------------------------------------------------------
  // PERTDISTANCERESSPEC
  // -------------------------------------------------------------------------
  {
    const std::string blockname = "PERTDISTANCERESSPEC";
    DEBUG(10, "PERTDISTANCERESSPEC block");
    buffer = m_block[blockname];

    if (buffer.size() <= 2) {
      if (count_colvar_specs(sim.param().colvarres.pert_bias_specs, "DISTANCE") > 0) {
        io::messages.add("PERTRESTRAINTS contains DISTANCE entries but no non-empty "
                         "PERTDISTANCERESSPEC block was found",
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
        io::messages.add("bad first data line in PERTDISTANCERESSPEC block: expected DISH DISC",
                         "In_Colvarres", io::message::error);
      }

      size_t distance_index = 0;
      bool indexed_distance_geometries = false;
      std::vector<std::pair<unsigned int, topology::distance_restraint_struct> >
        indexed_distances;
      if (it != to) ++it;
      for (; it != to; ++it) {
        if (!line_is_data(*it)) continue;

        int type1 = 0, type2 = 0;
        unsigned int geom_index = 0;
        std::vector<int> atom1, atom2;
        std::vector<std::string> tokens = tokens_from_line(*it);
        _lineStream.clear();
        _lineStream.str(*it);
        if (tokens.size() == 11) {
          _lineStream >> geom_index;
        }

        for (unsigned int i = 0; i < io::In_Colvarres::MAX_ATOMS; ++i) {
          unsigned int atom = 0;
          _lineStream >> atom;
          if (atom > topo.num_atoms()) {
            io::messages.add("PERTDISTANCERESSPEC block: atom number out of range",
                             "In_Colvarres", io::message::error);
          }
          if (atom > 0) atom1.push_back(atom - 1);
        }
        _lineStream >> type1;

        for (unsigned int i = 0; i < io::In_Colvarres::MAX_ATOMS; ++i) {
          unsigned int atom = 0;
          _lineStream >> atom;
          if (atom > topo.num_atoms()) {
            io::messages.add("PERTDISTANCERESSPEC block: atom number out of range",
                             "In_Colvarres", io::message::error);
          }
          if (atom > 0) atom2.push_back(atom - 1);
        }
        _lineStream >> type2;

        if (_lineStream.fail()) {
          io::messages.add("bad line in PERTDISTANCERESSPEC block: expected two virtual atom definitions",
                           "In_Colvarres", io::message::error);
          continue;
        }

        const simulation::Parameter::colvar_bias_spec *spec =
          geom_index > 0
            ? indexed_colvar_spec(sim.param().colvarres.pert_bias_specs,
                                  "DISTANCE", geom_index)
            : nth_colvar_spec(sim.param().colvarres.pert_bias_specs,
                              "DISTANCE", 2 * distance_index);
        if (spec == NULL) {
          if (geom_index > 0) {
            std::ostringstream msg;
            msg << "PERTDISTANCERESSPEC index " << geom_index
                << " has no matching PERTRESTRAINTS DISTANCE entry; geometry is ignored.";
            io::messages.add(msg.str(), "In_Colvarres", io::message::warning);
          }
          else {
            io::messages.add("PERTDISTANCERESSPEC contains more geometries than PERTRESTRAINTS DISTANCE entries; extra geometries are ignored.",
                             "In_Colvarres", io::message::warning);
          }
          ++distance_index;
          continue;
        }

        util::Virtual_Atom v1(util::virtual_type(type1), atom1, dish, disc);
        util::Virtual_Atom v2(util::virtual_type(type2), atom2, dish, disc);
        topology::distance_restraint_struct restraint(v1, v2,
                                                      spec->target, 1.0,
                                                      spec->rah);
        if (geom_index > 0) {
          indexed_distance_geometries = true;
          indexed_distances.push_back(std::make_pair(geom_index, restraint));
        }
        else {
          topo.distance_restraints().push_back(restraint);
        }
        ++distance_index;
      }

      if (indexed_distance_geometries) {
        std::set<unsigned int> appended;
        for (std::vector<simulation::Parameter::colvar_bias_spec>::const_iterator
             spec_it = sim.param().colvarres.pert_bias_specs.begin(),
             spec_to = sim.param().colvarres.pert_bias_specs.end();
             spec_it != spec_to; ++spec_it) {
          if (spec_it->type != "DISTANCE" || spec_it->index == 0
              || appended.count(spec_it->index)) continue;
          bool found = false;
          for (std::vector<std::pair<unsigned int, topology::distance_restraint_struct> >::const_iterator
               geom_it = indexed_distances.begin(),
               geom_to = indexed_distances.end(); geom_it != geom_to; ++geom_it) {
            if (geom_it->first != spec_it->index) continue;
            topo.distance_restraints().push_back(geom_it->second);
            appended.insert(spec_it->index);
            found = true;
            break;
          }
          if (!found) {
            std::ostringstream msg;
            msg << "PERTRESTRAINTS DISTANCE index " << spec_it->index
                << " has no matching PERTDISTANCERESSPEC geometry.";
            io::messages.add(msg.str(), "In_Colvarres", io::message::warning);
          }
        }
      }
    }
  }

  // -------------------------------------------------------------------------
  // PERTANGRESSPEC
  // -------------------------------------------------------------------------
  {
    const std::string blockname = "PERTANGRESSPEC";
    DEBUG(10, "PERTANGRESSPEC block");
    buffer = m_block[blockname];

    if (buffer.size() <= 2) {
      if (count_colvar_specs(sim.param().colvarres.pert_bias_specs, "ANGLE") > 0) {
        io::messages.add("PERTRESTRAINTS contains ANGLE entries but no non-empty "
                         "PERTANGRESSPEC block was found",
                         "In_Colvarres", io::message::error);
      }
    } else {
      block_read.insert(blockname);
      size_t angle_index = 0;
      bool indexed_angle_geometries = false;
      std::vector<std::pair<unsigned int, topology::angle_restraint_struct> >
        indexed_angles;
      for (std::vector<std::string>::const_iterator it = buffer.begin() + 1,
           to = buffer.end() - 1; it != to; ++it) {
        if (!line_is_data(*it)) continue;

        unsigned int geom_index = 0;
        unsigned int i = 0, j = 0, k = 0;
        std::vector<std::string> tokens = tokens_from_line(*it);
        _lineStream.clear();
        _lineStream.str(*it);
        if (tokens.size() == 4) {
          _lineStream >> geom_index;
        }
        _lineStream >> i >> j >> k;

        if (_lineStream.fail() || i == 0 || j == 0 || k == 0
            || i > topo.num_atoms() || j > topo.num_atoms() || k > topo.num_atoms()) {
          io::messages.add("bad line in PERTANGRESSPEC block: expected valid IPLR JPLR KPLR",
                           "In_Colvarres", io::message::error);
          continue;
        }

        const simulation::Parameter::colvar_bias_spec *spec =
          geom_index > 0
            ? indexed_colvar_spec(sim.param().colvarres.pert_bias_specs,
                                  "ANGLE", geom_index)
            : nth_colvar_spec(sim.param().colvarres.pert_bias_specs,
                              "ANGLE", 2 * angle_index);
        if (spec == NULL) {
          if (geom_index > 0) {
            std::ostringstream msg;
            msg << "PERTANGRESSPEC index " << geom_index
                << " has no matching PERTRESTRAINTS ANGLE entry; geometry is ignored.";
            io::messages.add(msg.str(), "In_Colvarres", io::message::warning);
          }
          else {
            io::messages.add("PERTANGRESSPEC contains more geometries than PERTRESTRAINTS ANGLE entries; extra geometries are ignored.",
                             "In_Colvarres", io::message::warning);
          }
          ++angle_index;
          continue;
        }

        topology::angle_restraint_struct restraint(i - 1, j - 1, k - 1,
                                                   spec->target, 1.0);
        if (geom_index > 0) {
          indexed_angle_geometries = true;
          indexed_angles.push_back(std::make_pair(geom_index, restraint));
        }
        else {
          topo.angle_restraints().push_back(restraint);
        }
        ++angle_index;
      }

      if (indexed_angle_geometries) {
        std::set<unsigned int> appended;
        for (std::vector<simulation::Parameter::colvar_bias_spec>::const_iterator
             spec_it = sim.param().colvarres.pert_bias_specs.begin(),
             spec_to = sim.param().colvarres.pert_bias_specs.end();
             spec_it != spec_to; ++spec_it) {
          if (spec_it->type != "ANGLE" || spec_it->index == 0
              || appended.count(spec_it->index)) continue;
          bool found = false;
          for (std::vector<std::pair<unsigned int, topology::angle_restraint_struct> >::const_iterator
               geom_it = indexed_angles.begin(),
               geom_to = indexed_angles.end(); geom_it != geom_to; ++geom_it) {
            if (geom_it->first != spec_it->index) continue;
            topo.angle_restraints().push_back(geom_it->second);
            appended.insert(spec_it->index);
            found = true;
            break;
          }
          if (!found) {
            std::ostringstream msg;
            msg << "PERTRESTRAINTS ANGLE index " << spec_it->index
                << " has no matching PERTANGRESSPEC geometry.";
            io::messages.add(msg.str(), "In_Colvarres", io::message::warning);
          }
        }
      }
    }
  }

  // -------------------------------------------------------------------------
  // PERTDIHEDRALRESSPEC
  // -------------------------------------------------------------------------
  {
    const std::string blockname = "PERTDIHEDRALRESSPEC";
    DEBUG(10, "PERTDIHEDRALRESSPEC block");
    buffer = m_block[blockname];

    if (buffer.size() <= 2) {
      if (count_colvar_specs(sim.param().colvarres.pert_bias_specs, "DIHEDRAL") > 0) {
        io::messages.add("PERTRESTRAINTS contains DIHEDRAL entries but no non-empty "
                         "PERTDIHEDRALRESSPEC block was found",
                         "In_Colvarres", io::message::error);
      }
    } else {
      block_read.insert(blockname);
      size_t dihedral_index = 0;
      bool indexed_dihedral_geometries = false;
      std::vector<std::pair<unsigned int, topology::dihedral_restraint_struct> >
        indexed_dihedrals;
      for (std::vector<std::string>::const_iterator it = buffer.begin() + 1,
           to = buffer.end() - 1; it != to; ++it) {
        if (!line_is_data(*it)) continue;

        unsigned int geom_index = 0;
        unsigned int i = 0, j = 0, k = 0, l = 0;
        std::vector<std::string> tokens = tokens_from_line(*it);
        _lineStream.clear();
        _lineStream.str(*it);
        if (tokens.size() == 5) {
          _lineStream >> geom_index;
        }
        _lineStream >> i >> j >> k >> l;

        if (_lineStream.fail() || i == 0 || j == 0 || k == 0 || l == 0
            || i > topo.num_atoms() || j > topo.num_atoms()
            || k > topo.num_atoms() || l > topo.num_atoms()) {
          io::messages.add("bad line in PERTDIHEDRALRESSPEC block: expected valid IPLR JPLR KPLR LPLR",
                           "In_Colvarres", io::message::error);
          continue;
        }

        const simulation::Parameter::colvar_bias_spec *spec =
          geom_index > 0
            ? indexed_colvar_spec(sim.param().colvarres.pert_bias_specs,
                                  "DIHEDRAL", geom_index)
            : nth_colvar_spec(sim.param().colvarres.pert_bias_specs,
                              "DIHEDRAL", 2 * dihedral_index);
        if (spec == NULL) {
          if (geom_index > 0) {
            std::ostringstream msg;
            msg << "PERTDIHEDRALRESSPEC index " << geom_index
                << " has no matching PERTRESTRAINTS DIHEDRAL entry; geometry is ignored.";
            io::messages.add(msg.str(), "In_Colvarres", io::message::warning);
          }
          else {
            io::messages.add("PERTDIHEDRALRESSPEC contains more geometries than PERTRESTRAINTS DIHEDRAL entries; extra geometries are ignored.",
                             "In_Colvarres", io::message::warning);
          }
          ++dihedral_index;
          continue;
        }

        topology::dihedral_restraint_struct restraint(i - 1, j - 1, k - 1, l - 1,
                                                      math::Pi, spec->target, 1.0);
        if (geom_index > 0) {
          indexed_dihedral_geometries = true;
          indexed_dihedrals.push_back(std::make_pair(geom_index, restraint));
        }
        else {
          topo.dihedral_restraints().push_back(restraint);
        }
        ++dihedral_index;
      }

      if (indexed_dihedral_geometries) {
        std::set<unsigned int> appended;
        for (std::vector<simulation::Parameter::colvar_bias_spec>::const_iterator
             spec_it = sim.param().colvarres.pert_bias_specs.begin(),
             spec_to = sim.param().colvarres.pert_bias_specs.end();
             spec_it != spec_to; ++spec_it) {
          if (spec_it->type != "DIHEDRAL" || spec_it->index == 0
              || appended.count(spec_it->index)) continue;
          bool found = false;
          for (std::vector<std::pair<unsigned int, topology::dihedral_restraint_struct> >::const_iterator
               geom_it = indexed_dihedrals.begin(),
               geom_to = indexed_dihedrals.end(); geom_it != geom_to; ++geom_it) {
            if (geom_it->first != spec_it->index) continue;
            topo.dihedral_restraints().push_back(geom_it->second);
            appended.insert(spec_it->index);
            found = true;
            break;
          }
          if (!found) {
            std::ostringstream msg;
            msg << "PERTRESTRAINTS DIHEDRAL index " << spec_it->index
                << " has no matching PERTDIHEDRALRESSPEC geometry.";
            io::messages.add(msg.str(), "In_Colvarres", io::message::warning);
          }
        }
      }
    }
  }

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

  for (std::map<std::string, std::vector<std::string> >::const_iterator
       it = m_block.begin(), to = m_block.end(); it != to; ++it) {
    if (block_read.count(it->first) == 0 && it->second.size()) {
      io::messages.add(" block " + it->first + " not supported!",
                       "In_Colvarres", io::message::warning);
    }
  }

  const unsigned int n_distance_specs = count_colvar_specs(sim, "DISTANCE");
  const unsigned int n_pertdistance_specs =
    count_colvar_specs(sim.param().colvarres.pert_bias_specs, "DISTANCE") / 2;
  const unsigned int n_distance_specs_total =
    n_distance_specs + n_pertdistance_specs;
  if (n_distance_specs_total > 0
      && n_distance_specs_total != topo.distance_restraints().size()) {
    std::ostringstream msg;
    msg << "RESTRAINTS+PERTRESTRAINTS/DISTANCERESSPEC+PERTDISTANCERESSPEC mismatch: "
        << n_distance_specs_total
        << " DISTANCE bias specification(s), but "
        << topo.distance_restraints().size()
        << " DISTANCE geometry specification(s).";
    io::messages.add(msg.str(), "In_Colvarres", io::message::warning);
  }

  const unsigned int n_coord_specs = count_colvar_specs(sim, "COORDNUM");
  const unsigned int n_pertcoord_specs =
    count_colvar_specs(sim.param().colvarres.pert_bias_specs, "COORDNUM") / 2;
  const unsigned int n_coord_specs_total =
    n_coord_specs + n_pertcoord_specs;
  if (n_coord_specs_total > 0
      && n_coord_specs_total != topo.coordnum_restraint().size()) {
    std::ostringstream msg;
    msg << "RESTRAINTS+PERTRESTRAINTS/COORDNUMRESSPEC+PERTCOORDNUMRESSPEC mismatch: "
        << n_coord_specs_total
        << " COORDNUM bias specification(s), but "
        << topo.coordnum_restraint().size()
        << " COORDNUM geometry specification(s).";
    io::messages.add(msg.str(), "In_Colvarres", io::message::warning);
  }

  const unsigned int n_angle_specs = count_colvar_specs(sim, "ANGLE");
  const unsigned int n_pertangle_specs =
    count_colvar_specs(sim.param().colvarres.pert_bias_specs, "ANGLE") / 2;
  const unsigned int n_angle_specs_total =
    n_angle_specs + n_pertangle_specs;
  if (n_angle_specs_total > 0
      && n_angle_specs_total != topo.angle_restraints().size()) {
    std::ostringstream msg;
    msg << "RESTRAINTS+PERTRESTRAINTS/ANGRESSPEC+PERTANGRESSPEC mismatch: "
        << n_angle_specs_total
        << " ANGLE bias specification(s), but "
        << topo.angle_restraints().size()
        << " ANGLE geometry specification(s).";
    io::messages.add(msg.str(), "In_Colvarres", io::message::warning);
  }

  const unsigned int n_dihedral_specs = count_colvar_specs(sim, "DIHEDRAL");
  const unsigned int n_pertdihedral_specs =
    count_colvar_specs(sim.param().colvarres.pert_bias_specs, "DIHEDRAL") / 2;
  const unsigned int n_dihedral_specs_total =
    n_dihedral_specs + n_pertdihedral_specs;
  if (n_dihedral_specs_total > 0
      && n_dihedral_specs_total != topo.dihedral_restraints().size()) {
    std::ostringstream msg;
    msg << "RESTRAINTS+PERTRESTRAINTS/DIHEDRALRESSPEC+PERTDIHEDRALRESSPEC mismatch: "
        << n_dihedral_specs_total
        << " DIHEDRAL bias specification(s), but "
        << topo.dihedral_restraints().size()
        << " DIHEDRAL geometry specification(s).";
    io::messages.add(msg.str(), "In_Colvarres", io::message::warning);
  }
}
