/*
 * This file is part of GROMOS.
 *
 * Copyright (c) 2011, 2012, 2016, 2018, 2021, 2023 Biomos b.v.
 * See <https://www.gromos.net> for details.
 */

/**
 * @file in_colvarres_bias.cc
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

void io::In_Colvarres::read_restraints(topology::Topology &topo, simulation::Simulation &sim, std::ostream &os)
{
  std::vector<std::string> buffer;
  std::set<std::string> &block_read = m_colvar_blocks_read;

  // -------------------------------------------------------------------------
  // RESTRAINTS
  // -------------------------------------------------------------------------
  {
    const std::string blockname = "RESTRAINTS";
    DEBUG(10, "RESTRAINTS block");
    buffer = m_block[blockname];

    if (buffer.size() <= 2) {
      const bool has_pert_restraints = m_block["PERTRESTRAINTS"].size() > 2
                                    || m_block["EDSRESTRAINTS"].size() > 2;
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
          io::messages.add("RESTRAINTS block: TYPE must be 0 DISTANCE, 1 ANGLE, 2 DIHEDRAL, 3 COORDNUM, or 4 DISTANCE_DIFFERENCE",
                           "In_Colvarres", io::message::error);
          continue;
        }
        if (spec.type == "DISTANCE_DIFFERENCE" && spec.index == 0) {
          io::messages.add(blockname + ": DISTANCE_DIFFERENCE requires a positive INDX.",
                           "In_Colvarres", io::message::error);
        }
        const double input_target = spec.target;
        if (spec.type == "ANGLE" || spec.type == "DIHEDRAL") {
          spec.target *= math::Pi / 180.0;
          spec.linear_tail *= math::Pi / 180.0;
        }

        spec.rah = rah;
        validate_rah(rah, spec.type, blockname);

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
}

void io::In_Colvarres::read_pertrestraints(topology::Topology &topo, simulation::Simulation &sim, std::ostream &os)
{
  std::vector<std::string> buffer;
  std::set<std::string> &block_read = m_colvar_blocks_read;

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
           << std::setw(8) << "N" << std::setw(8) << "M"
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
        if (tokens.size() == 11 || tokens.size() == 13) {
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
        else if (tokens.size() == 10) {
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

        int hidden_n = 0, hidden_m = 0;
        if (tokens.size() == 13)
          _lineStream >> hidden_n >> hidden_m;
        const bool parse_error = _lineStream.fail();
        std::string trailing;
        const bool extra = bool(_lineStream >> trailing);
        // Validate integer suffixes and exact column counts, including partial pairs.
        if ((tokens.size() != 10 && tokens.size() != 11
             && tokens.size() != 13)
            || parse_error || extra || hidden_n < 0 || hidden_m < 0) {
          io::messages.add("PERTRESTRAINTS: expected 10 legacy columns, 11 indexed columns, or 13 indexed columns with nonnegative integer N M.",
                           "In_Colvarres", io::message::error);
          continue;
        }
        specA.hidden_n = specB.hidden_n = hidden_n;
        specA.hidden_m = specB.hidden_m = hidden_m;

        specB.index = specA.index;
        specA.type = colvar_type_from_int(type);
        specB.type = specA.type;
        if (specA.type.empty()) {
          io::messages.add("PERTRESTRAINTS block: TYPE must be 0 DISTANCE, 1 ANGLE, 2 DIHEDRAL, 3 COORDNUM, or 4 DISTANCE_DIFFERENCE",
                           "In_Colvarres", io::message::error);
          continue;
        }
        if (specA.type == "DISTANCE_DIFFERENCE" && specA.index == 0) {
          io::messages.add(blockname + ": DISTANCE_DIFFERENCE requires a positive INDX.",
                           "In_Colvarres", io::message::error);
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
        validate_rah(rah, specA.type, blockname);
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
             << std::setw(8) << specA.hidden_n << std::setw(8) << specA.hidden_m
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
}

