/*
 * This file is part of GROMOS.
 * 
 * Copyright (c) 2011, 2012, 2016, 2018, 2021, 2023 Biomos b.v.
 * See <https://www.gromos.net> for details.
 * 
 * GROMOS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

/**
 * @file in_leus.cc
 * implements methods of In_Localelevspec
 */
#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../simulation/parameter.h"
#include "../../interaction/interaction_types.h"
#include "../../configuration/configuration.h"

#include "../../io/instream.h"
#include "../../io/blockinput.h"
#include "../../io/topology/in_leus.h"
#include "../../io/configuration/in_configuration.h"
#include "../../util/le_coordinate.h"
#include "../../util/le_colvar_coordinate.h"

#undef MODULE
#undef SUBMODULE
#define MODULE io
#define SUBMODULE topology

namespace {
bool read_le_sites(std::istream &input,size_t count,double dish,double disc,
                   const topology::Topology &topo,std::vector<util::Virtual_Atom> &sites) {
  for (size_t endpoint=0;endpoint<count;++endpoint) {
    int type=0;
    long long size=0;
    if (!(input>>type>>size) || size<=0
        || static_cast<unsigned long long>(size)>topo.num_atoms()) return false;
    long long expected=0;
    switch (type) {
      case 0: expected=1; break;
      case 1: expected=4; break;
      case 2: case 3: case 4: case 6: case 8: expected=3; break;
      case 5: case 7: expected=2; break;
      case -1: case -2: expected=size; break;
      default: return false;
    }
    if (size!=expected) return false;
    std::vector<int> atoms;
    std::set<int> unique;
    for (long long i=0;i<size;++i) {
      long long atom=0;
      if (!(input>>atom) || atom<=0 || static_cast<unsigned long long>(atom)>topo.num_atoms()
          || atom-1>std::numeric_limits<int>::max()
          || !unique.insert(static_cast<int>(atom-1)).second) return false;
      atoms.push_back(static_cast<int>(atom-1));
    }
    sites.emplace_back(util::virtual_type(type),atoms,dish,disc);
  }
  std::string extra;
  return !(input>>extra);
}
}

/**
 * @section localelevspec LOCALELEVSPEC block
 * LOCALELEVSPEC attaches geometric coordinates to umbrella potentials.
 * Existing types 1/2/6 are dihedral, distance and distance field. Type 3
 * defines an angle with I J K 0; type 4 defines d(I,J)-d(K,L). Type 5
 * uses RCUT N M NA NB followed by the two real-atom lists instead of I J K L.
 * See doc/colvar_le_types.md for units, grids and complete input examples.
 *
 * The block is read from the local elevation specification file
 * (\@led).
 *
 * @verbatim
LOCALELEVSPEC
#  NLEPID  TYPE  IPLE  JPLE  KPLE  LPLE
        1     1   1     2     3     4
        1     1   2     3     4     5
        2     1   1     2     3     4
END
@endverbatim
 */
void
io::In_Localelevspec::read(topology::Topology& topo,
		    simulation::Simulation & sim,
		    std::ostream & os){

  DEBUG(7, "reading in a local elevation definition file");

  std::vector<std::string> buffer;

  { // LOCALELEVSPEC

    buffer = m_block["LOCALELEVSPEC"];
    DEBUG(10, "LOCALELEVSPEC block : " << buffer.size());

    if (!buffer.size()){
      io::messages.add("no LOCALELEVSPEC block in local elevation definition file",
		       "In_Localelevspec", io::message::error);
      return;
    }

    std::vector<std::string>::const_iterator it = buffer.begin()+1,
      to = buffer.end()-1;

    DEBUG(10, "reading in LOCALELEVSPEC data");

    int id = 0;
    unsigned int n = 0, i = 0, j = 0, k = 0, l = 0, t = 0;
    for(n=0; it != to; ++n, ++it){
      DEBUG(11, "\tnr " << n);

      _lineStream.clear();
      _lineStream.str(*it);

      _lineStream >> id >> t;
      if (_lineStream.fail() || id <= 0) {
        io::messages.add("LOCALELEVSPEC: expected a positive umbrella ID and TYPE.",
                         "In_Localelevspec",io::message::error);
        return;
      }
      const auto geometry_start=_lineStream.tellg();
      std::string format;
      _lineStream >> format;
      if (format=="VA") {
        double dish=0,disc=0,rcut=0;
        int nn=0,mm=0;
        size_t count=t==2 ? 2 : (t==3 ? 3 : 4),first=0;
        if (!(_lineStream>>dish>>disc) || !std::isfinite(dish) || !std::isfinite(disc)
            || dish<=0 || disc<=0 || t<1 || t>5) {
          io::messages.add("LOCALELEVSPEC VA: types 1..5 require positive finite DISH DISC.",
                           "In_Localelevspec",io::message::error);
          return;
        }
        if (t==util::Umbrella::vt_coordnum) {
          long long na=0,nb=0;
          if (!(_lineStream>>rcut>>nn>>mm>>na>>nb) || !std::isfinite(rcut) || rcut<=0
              || nn<=0 || nn>std::numeric_limits<int>::max()/2 || mm<=nn
              || na<=0 || nb<=0 || static_cast<unsigned long long>(na)>topo.num_atoms()
              || static_cast<unsigned long long>(nb)>topo.num_atoms()) {
            io::messages.add("LOCALELEVSPEC VA CN: invalid switching parameters or group sizes.",
                             "In_Localelevspec",io::message::error);
            return;
          }
          first=static_cast<size_t>(na);
          count=first+static_cast<size_t>(nb);
        }
        std::vector<util::Virtual_Atom> sites;
        if (!read_le_sites(_lineStream,count,dish,disc,topo,sites)) {
          io::messages.add("LOCALELEVSPEC VA: expected VATYPE NATOMS ATOMS for each endpoint; check types, sizes, indices, duplicates and extra fields.",
                           "In_Localelevspec",io::message::error);
          return;
        }
        if (t==util::Umbrella::vt_coordnum) {
          for (int group=0;group<2;++group) {
            std::set<std::string> unique;
            for (size_t i=group ? first : 0;i<(group ? count : first);++i) {
              std::ostringstream key;
              key<<sites[i].type();
              for (int atom=0;atom<sites[i].size();++atom) key<<' '<<sites[i].atom(atom);
              if (!unique.insert(key.str()).second) {
                io::messages.add("LOCALELEVSPEC VA CN: duplicate endpoint within a group.",
                                 "In_Localelevspec",io::message::error);
                return;
              }
            }
          }
        }
        topo.le_coordinates().push_back(new util::LE_Colvar_Coordinate(id,t,sites,first,rcut,nn,mm));
        continue;
      }
      _lineStream.clear();
      _lineStream.seekg(geometry_start);
      if (t == util::Umbrella::vt_coordnum) {
        double rcut=0;
        int nn=0,mm=0;
        long long nfirst=0,nsecond=0;
        _lineStream >> rcut >> nn >> mm >> nfirst >> nsecond;
        if (_lineStream.fail() || !std::isfinite(rcut) || rcut <= 0 || nn <= 0
            || nn > std::numeric_limits<int>::max()/2 || mm <= nn
            || nfirst <= 0 || nsecond <= 0
            || static_cast<unsigned long long>(nfirst) > topo.num_atoms()
            || static_cast<unsigned long long>(nsecond) > topo.num_atoms()) {
          io::messages.add("LOCALELEVSPEC CN: expected RCUT > 0, 0 < N < M, positive group sizes, then atom lists.",
                           "In_Localelevspec",io::message::error);
          return;
        }
        std::vector<unsigned int> first,second;
        for (int group=0;group<2;++group) {
          std::set<unsigned int> unique;
          auto &atoms=group==0 ? first : second;
          const long long count=group==0 ? nfirst : nsecond;
          for (long long atom=0;atom<count;++atom) {
            long long index=0;
            _lineStream >> index;
            if (_lineStream.fail() || index <= 0
                || static_cast<unsigned long long>(index) > topo.num_atoms()
                || !unique.insert(static_cast<unsigned int>(index)).second) {
              io::messages.add("LOCALELEVSPEC CN: atom index out of range, missing, or duplicated within its group.",
                               "In_Localelevspec",io::message::error);
              return;
            }
            atoms.push_back(static_cast<unsigned int>(index-1));
          }
        }
        std::string extra;
        if (_lineStream >> extra) {
          io::messages.add("LOCALELEVSPEC CN: unexpected extra fields.","In_Localelevspec",io::message::error);
          return;
        }
        topo.le_coordinates().push_back(new util::LE_Colvar_Coordinate(id,first,second,rcut,nn,mm));
        continue;
      }
      _lineStream >> i >> j >> k >> l;
      if (_lineStream.fail()){
        io::messages.add("bad line in LOCALELEVSPEC block; expected format:\nLOCALELEVSPEC\n#  NLEPID  TYPE  IPLE  JPLE  KPLE  LPLE\nEND\n",
                         "In_Localelevspec", io::message::error);
        return;
      }

      if (t == util::Umbrella::vt_angle || t == util::Umbrella::vt_distance_difference) {
        std::string extra;
        if ((t == util::Umbrella::vt_angle && l != 0) || (_lineStream >> extra)) {
          io::messages.add("LOCALELEVSPEC angle/distance difference: expected exactly four atom fields; pad the angle's fourth field with 0.",
                           "In_Localelevspec",io::message::error);
          return;
        }
      }

      // convert to gromos
      --i; --j; --k; --l;

      switch(t){
        case util::Umbrella::vt_angle:
        case util::Umbrella::vt_distance_difference:
        {
          std::vector<unsigned int> atoms={i,j,k};
          if (t == util::Umbrella::vt_distance_difference) atoms.push_back(l);
          for (auto atom : atoms) {
            if (atom >= topo.num_atoms()) {
              io::messages.add("LOCALELEVSPEC angle/distance difference: atom index out of range.",
                               "In_Localelevspec",io::message::error);
              return;
            }
          }
          topo.le_coordinates().push_back(new util::LE_Colvar_Coordinate(id,t,atoms));
          break;
        }
	case util::Umbrella::vt_distance:
        {
          // i and j are unsigned int, if the user specified a negative value
          // this is interpreted as a very large value
          // no need to check for i<0 and j<0, as this will never happen
          // some compilers complain about this
          if (i >= topo.num_atoms() || j >= topo.num_atoms()){
            std::ostringstream msg;
            msg << "Distance (" << i+1 << "-" << j+1 
                << ") atom indices out of range.";
            io::messages.add(msg.str(), "In_Localelevspec", io::message::error);
            return;
          }
	  util::LE_Colvar_Coordinate *dist = new util::LE_Colvar_Coordinate(id,2,{i,j});
          topo.le_coordinates().push_back(dist);
          break;
        }
        case util::Umbrella::vt_dihedral:
        {
          // i, j, k and l are unsigned int, if the user specified a negative 
          // value this is interpreted as a very large value
          // no need to check for i<0 and j<0, as this will never happen
          // some compilers complain about this
          if (i >= topo.num_atoms() || j >= topo.num_atoms() ||
              k >= topo.num_atoms() || l >= topo.num_atoms()){
            std::ostringstream msg;
            msg << "Dihedral (" << i+1 << "-" << j+1 << "-" << k+1 << "-" << l+1
                 << ") atom indices out of range.";
            io::messages.add(msg.str(), "In_Localelevspec", io::message::error);
            return;
          }
          util::LE_Colvar_Coordinate *dih = new util::LE_Colvar_Coordinate(id,1,{i,j,k,l});
          topo.le_coordinates().push_back(dih);
          break;
        }
	case util::Umbrella::vt_distancefield:
	{
	  io::messages.add("Distancefield definition according to DISTANCEFIELD block in input file. Atom numbers ignored.",
			   "In_Localelevspec", io::message::notice);
	  if(sim.param().distancefield.distancefield!=1){
	    io::messages.add("When using local elevation with type 6 (distance field), the DISTANCEFIELD block should be switched on", 
			     "In_Localelevspec", io::message::error);
	    return;
	  }
	  if(topo.disfield_restraints().K!=0.0){
	    io::messages.add("Force constant in DISTANCEFIELD block != 0.0. This means that the distance field restraint is applied on top of the local elevation potential", 
			     "In_Localelevspec", io::message::warning);
	  }
	  
	  // the above message means that it needs to know about the parameters and about the topology
	  util::LE_DistanceField_Coordinate * df = new util::LE_DistanceField_Coordinate(id, topo, sim);
	  topo.le_coordinates().push_back(df);
	  break;
	}
        default:
        {
          std::ostringstream msg;
          msg << "LE Coordinate type " << t << " unknown";
          io::messages.add(msg.str(), "In_Localelevspec", io::message::error);
          return;
        }
      }
    }
  } // LOCALELEVSPEC
}

/**
 * @section leusbias LEUSBIAS block
 *
 * The block is read from the local elevation umbrella sampling database file
 * (\@lud) or the configuration. The reading from the input file is not implemented
 * in md++.
 *
 * @verbatim
LEUSBIAS
# NUMB:    number of umbrellas
# NLEPID:  ID of the umbrella
# NDIM:    number of dimensions of the umbrella
# CLES:    force constant of the umbrella
# VARTYPE: LE coordinate type (1 dihedral, 2 distance, 3 angle, 4 distance difference, 5 CN, 6 distance field)
# NTLEFU:  functional form of the potential in this dimension
# WLES:    width of the potential in grid units
# RLES:    cutoff in grid units
# NGRID:   number of grid points
# GRIDMIN: minumum value for the grid (origin)
# GRIDMAX: maximum value for the grid.
#    if GRIDMIN == GRIDMAX: periodic dihedral; bounded angles require distinct bounds
# NCONLE:  number of (different) configurations visited
# NVISLE:  number of times the confiuguration was visited
# ICONF(1..NDIM): the coordinates of the grid point (1..GRIDMAX)
#
# NUMB
2
# NLEPID NDIM CLES
         3         1    0.005000000
# VARTYPE(N) NTLEFU(N) WLES(N) RLES(N) NGRID(N) GRIDMIN(N) GRIDMAX(N)
         1         0    1.000000000    1.000000000        32    0.000000000    0.000000000
# NCONLE
         3
# NVISLE ICONF(1..NDIM)
         2         3
         5         4
         4         6
# NLEPID NDIM CLES
         7         2    0.004000000
# VARTYPE(N) NTLEFU(N) WLES(N) RLES(N) NGRID(N) GRIDMIN(N) GRIDMAX(N)
         1         1    1.000000000    1.000000000        16   15.000000000   15.000000000
         1         0    1.000000000    1.000000000        10    0.000000000    0.000000000
# NCONLE
         1
# NVISLE ICONF(1..NDIM)
         1         3         7
END
@endverbatim
 */
void
io::In_LEUSBias::read(topology::Topology& topo,
        configuration::Configuration & conf,
        simulation::Simulation & sim,
        std::ostream & os) {

  DEBUG(7, "reading in a local elevation definition file");

  std::vector<std::string> buffer;

  if (sim.param().localelev.read) { // LOCALELEVSPEC
    buffer = m_block["LEUSBIAS"];
    DEBUG(10, "LOCALELEVSPEC block : " << buffer.size());

    if (!buffer.size()){
      io::messages.add("no LEUSBIAS block in LEUS database file",
		       "In_Localelevspec", io::message::error);
      return;
    }

    // remove the title from the buffer. the title to the buffer
    std::vector<std::string> new_buffer;
    new_buffer.insert(new_buffer.begin(), buffer.begin()+1, buffer.end());

    io::In_Configuration::_read_leusbias(conf.special().umbrellas, new_buffer, false);

  } // LEUSBIAS
}



