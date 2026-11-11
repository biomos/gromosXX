/*
 * This file is part of GROMOS.
 *
 * Copyright (c) 2011, 2012, 2016, 2018, 2021, 2023 Biomos b.v.
 * See <https://www.gromos.net> for details.
 */

/**
 * @file in_colvarres_angular.cc
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

void io::In_Colvarres::read_angresspec(topology::Topology &topo, simulation::Simulation &sim, std::ostream &os)
{
  std::vector<std::string> buffer;
  std::set<std::string> &block_read = m_colvar_blocks_read;

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
          topo.colvar_angle_geometries().push_back(restraint);
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
            topo.colvar_angle_geometries().push_back(geom_it->second);
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
}

void io::In_Colvarres::read_dihedralresspec(topology::Topology &topo, simulation::Simulation &sim, std::ostream &os)
{
  std::vector<std::string> buffer;
  std::set<std::string> &block_read = m_colvar_blocks_read;

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
          topo.colvar_dihedral_geometries().push_back(restraint);
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
            topo.colvar_dihedral_geometries().push_back(geom_it->second);
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
}

void io::In_Colvarres::read_pertangresspec(topology::Topology &topo, simulation::Simulation &sim, std::ostream &os)
{
  std::vector<std::string> buffer;
  std::set<std::string> &block_read = m_colvar_blocks_read;

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
          topo.colvar_angle_geometries().push_back(restraint);
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
            topo.colvar_angle_geometries().push_back(geom_it->second);
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
}

void io::In_Colvarres::read_pertdihedralresspec(topology::Topology &topo, simulation::Simulation &sim, std::ostream &os)
{
  std::vector<std::string> buffer;
  std::set<std::string> &block_read = m_colvar_blocks_read;

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
          topo.colvar_dihedral_geometries().push_back(restraint);
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
            topo.colvar_dihedral_geometries().push_back(geom_it->second);
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
}

namespace io {
namespace colvar_input {

bool matching_wall_settings(const simulation::Parameter::colvar_bias_spec &a,
                            const simulation::Parameter::colvar_bias_spec &b)
{
  return a.k == b.k && a.method == b.method && a.tau == b.tau
      && a.virial == b.virial && a.linear_tail == b.linear_tail
      && a.hidden_n == b.hidden_n && a.hidden_m == b.hidden_m;
}

void pair_angular_walls(std::vector<simulation::Parameter::colvar_bias_spec> &specs,
    const std::vector<simulation::Parameter::colvar_bias_spec> &ordinary_specs,
    const topology::Topology &topo, bool perturbed)
{
  const std::string block = perturbed ? "PERTRESTRAINTS" : "RESTRAINTS";
  const size_t states = perturbed ? 2 : 1;
  size_t angle = perturbed ? count_colvar_specs(ordinary_specs, "ANGLE") : 0;
  size_t dihedral = perturbed ? count_colvar_specs(ordinary_specs, "DIHEDRAL") : 0;
  std::map<std::vector<int>, std::vector<size_t> > groups;
  for (size_t i = 0; i < specs.size(); i += states) {
    const auto &spec = specs[i];
    std::vector<int> geometry;
    if (spec.type == "ANGLE") {
      const size_t index = angle++;
      if (index < topo.colvar_angle_geometries().size()) {
        const auto &g = topo.colvar_angle_geometries()[index];
        geometry = {1, g.i, g.j, g.k};
      }
    } else if (spec.type == "DIHEDRAL") {
      const size_t index = dihedral++;
      if (index < topo.colvar_dihedral_geometries().size()) {
        const auto &g = topo.colvar_dihedral_geometries()[index];
        geometry = {2, g.i, g.j, g.k, g.l};
      }
    } else {
      continue;
    }
    if (spec.rah != -1 && spec.rah != 1) continue;
    if (geometry.empty()) {
      io::messages.add(block + ": angular wall has no matching geometry.",
                       "In_Colvarres", io::message::error);
      continue;
    }
    groups[geometry].push_back(i);
  }
  for (const auto &group : groups) {
    const auto &rows = group.second;
    if (rows.size() != 2 || specs[rows[0]].rah == specs[rows[1]].rah) {
      io::messages.add(block + ": angular walls require exactly one RAH -1 and one RAH +1 on the same ordered geometry.",
                       "In_Colvarres", io::message::error);
      continue;
    }
    const size_t lower = specs[rows[0]].rah == -1 ? rows[0] : rows[1];
    const size_t upper = specs[rows[0]].rah == 1 ? rows[0] : rows[1];
    if (lower + states > specs.size() || upper + states > specs.size()) {
      io::messages.add(block + ": incomplete angular wall A/B states.",
                       "In_Colvarres", io::message::error);
      continue;
    }
    bool valid = true;
    for (size_t state = 0; state < states; ++state) {
      const auto &lo = specs[lower + state];
      const auto &hi = specs[upper + state];
      if (!matching_wall_settings(lo, hi)) {
        io::messages.add(block + ": paired angular walls must have identical CVK, METHOD, TAU, VIRIAL, CVLIN and hidden exponents.",
                         "In_Colvarres", io::message::error);
        valid = false;
      }
      if (!std::isfinite(lo.target) || !std::isfinite(hi.target)) {
        valid = false;
        io::messages.add(block + ": angular wall targets must be finite.",
                         "In_Colvarres", io::message::error);
      } else if (lo.type == "ANGLE") {
        if (lo.target < 0 || hi.target > math::Pi || lo.target >= hi.target) {
          io::messages.add(block + ": ANGLE walls require 0 <= lower TARGET < upper TARGET <= 180 degrees.",
                           "In_Colvarres", io::message::error);
          valid = false;
        }
      } else {
        const double period = 2.0 * math::Pi;
        double width = hi.target - lo.target;
        width -= period * std::floor(width / period);
        if (width < 1e-12 || period - width < 1e-12) {
          io::messages.add(block + ": DIHEDRAL wall targets must define a nonzero arc shorter than 360 degrees.",
                           "In_Colvarres", io::message::error);
          valid = false;
        }
        if (lo.averaging != 0) {
          io::messages.add(block + ": paired DIHEDRAL walls currently require METHOD 0; periodic averaging is not implemented.",
                           "In_Colvarres", io::message::error);
          valid = false;
        }
      }
    }
    if (perturbed && specs[lower].type == "DIHEDRAL") {
      const double period = 2.0 * math::Pi;
      const double winding_a = std::floor((specs[upper].target - specs[lower].target) / period);
      const double winding_b = std::floor((specs[upper + 1].target - specs[lower + 1].target) / period);
      if (winding_a != winding_b) {
        io::messages.add(block + ": use consistent unwrapped DIHEDRAL wall targets across A/B states so the interpolated arc does not collapse.",
                         "In_Colvarres", io::message::error);
        valid = false;
      }
    }
    if (!valid) continue;
    for (size_t state = 0; state < states; ++state) {
      specs[lower + state].paired_wall = specs[upper + state].paired_wall = true;
      specs[lower + state].wall_partner_target = specs[upper + state].target;
      specs[upper + state].wall_partner_target = specs[lower + state].target;
    }
  }
}

} // namespace colvar_input
} // namespace io
