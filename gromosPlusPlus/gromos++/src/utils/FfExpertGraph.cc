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
#include "FfExpertGraph.h"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "../gcore/AtomTopology.h"
#include "../gcore/BbSolute.h"
#include "../gcore/Bond.h"
#include "../gcore/BuildingBlock.h"
#include "../gcore/MoleculeTopology.h"
#include "../gio/Ginstream.h"
#include "../gromos/Exception.h"

namespace utils {
/**
 * write a vertex
 */
std::ostream &operator<<(std::ostream &os, const utils::Vertex &v) {
  os << "V(" << v.id << "," << v.name << ") ";
  return os;
}

/**
 * write the graph in the graphviz DOT format.
 * You can polt it with: neato -Tpng file.dot > file.png
 */
std::ostream &operator<<(std::ostream &os, const utils::FfExpertGraph &g) {
  g.print(os);
  return os;
}
} // namespace utils

void utils::FfExpertGraph::print(
    std::ostream &os, const utils::Vertex::equality_enum criterion) const {
  os << "graph G  {" << std::endl;
  for (size_t i = 0; i < vertices().size(); ++i) {
    // loop over neighbors
    for (size_t j = i; j < vertices().size(); ++j) {
      const double &w = weight(vertices()[i], vertices()[j]);
      if (w)
        os << i << "--" << j << "[len=" << w << "];" << std::endl;
    } // neighbors

    os << i << " [label=";
    switch (criterion) {
    case utils::Vertex::equal_type:
      os << vertices()[i].type;
      break;
    case utils::Vertex::equal_iac:
      os << vertices()[i].iac;
      break;
    case utils::Vertex::equal_charge:
      os << vertices()[i].charge;
      break;
    case utils::Vertex::equal_name:
    default: // name is default
      os << vertices()[i].name;
      break;
    }
    os << "];" << std::endl;
  }
  os << "}" << std::endl;
}

/**
 * gets the neighbors of a given vertex
 */
std::set<utils::Vertex>
utils::FfExpertGraph::get_neighbors(const utils::Vertex &v) const {
  std::set<utils::Vertex> neighbors;
  for (std::vector<utils::Vertex>::const_iterator it = vertices().begin(),
                                                  to = vertices().end();
       it != to; ++it) {
    if (connected(v, *it))
      neighbors.insert(*it);
  }
  return neighbors;
}
/**
 * cuts out a subgraph containing the centre and radius neighbors
 */
utils::FfExpertGraph
utils::FfExpertGraph::cut_subgraph(const utils::Vertex &centre,
                                   size_t radius) const {
  utils::FfExpertGraph g;

  // radius has to be a bit bigger - otherwise the result is not intuitive...
  radius++;

  std::vector<std::set<utils::Vertex>> level(radius + 1);
  level[0].insert(centre);
  for (size_t i = 0; i < radius; ++i) {
    for (std::set<utils::Vertex>::const_iterator it = level[i].begin(),
                                                 to = level[i].end();
         it != to; ++it) {

      // only insert if it is not there yet!
      if (std::find(g.vertices().begin(), g.vertices().end(), *it) ==
          g.vertices().end())
        g.vertices().push_back(*it);

      std::set<utils::Vertex> neighbors = get_neighbors(*it);
      for (std::set<utils::Vertex>::const_iterator
               neighbor_it = neighbors.begin(),
               neighbor_to = neighbors.end();
           neighbor_it != neighbor_to; ++neighbor_it) {

        // add it to the next level if it was not considered yet.
        if (std::find(g.vertices().begin(), g.vertices().end(), *neighbor_it) ==
            g.vertices().end())
          level[i + 1].insert(*neighbor_it);

        // connect it
        g.connect(*it, *neighbor_it, weight(*it, *neighbor_it));
      } // neighbors
    } // vertices in level
  } // level

  // return the sub graph
  return g;
}

/**
 * takes two nodes of two graphs and compares them and the neighbors
 */
bool utils::FfExpertGraph::match_neighbors(
    const utils::Vertex &v1, const utils::FfExpertGraph &g1,
    const utils::Vertex &v2, const utils::FfExpertGraph &g2,
    const utils::Vertex::equality_enum criterion, const unsigned int radius,
    const utils::Vertex &from_v1, const utils::Vertex &from_v2,
    const unsigned int my_radius) {
  typedef std::set<utils::Vertex> v_set;

  /* debug only
  std::cout << "Comparing " << v1 << " and " << v2 << std::endl;
   */

  // match the criterion
  if (!v1.equals(v2, criterion))
    return false;

  /* debug only
  if (my_radius != 0) {
    std::cout << "weight 1: " << g1.weight(v1, from_v1) << std::endl;
    std::cout << "weight 2: " << g1.weight(v2, from_v2) << std::endl;
  }*/

  if (my_radius != 0 && g1.weight(v1, from_v1) != g2.weight(v2, from_v2))
    return false;

  // don't compare neighbors at the end.
  if (my_radius == radius)
    return true;

  v_set neighbors_v1 = g1.get_neighbors(v1);
  v_set neighbors_v2 = g2.get_neighbors(v2);
  // remove precding vertex from neighbor list
  if (my_radius != 0) {
    neighbors_v1.erase(from_v1);
    neighbors_v2.erase(from_v2);
  }

  /* debug only
  std::cout << "Neighbors 1: ";
  std::copy(neighbors_v1.begin(), neighbors_v1.end(),
  std::ostream_iterator<utils::Vertex>(std::cout)); std::cout << std::endl;
  std::cout << "Neighbors 2: ";
  std::copy(neighbors_v2.begin(), neighbors_v2.end(),
  std::ostream_iterator<utils::Vertex>(std::cout)); std::cout << std::endl <<
  "-------" << std::endl;
   */

  if (neighbors_v1.size() != neighbors_v2.size())
    return false; // they can't match in this case!

  for (v_set::const_iterator it_v1 = neighbors_v1.begin(),
                             to_v1 = neighbors_v1.end();
       it_v1 != to_v1; ++it_v1) {

    for (v_set::const_iterator it_v2 = neighbors_v2.begin(),
                               to_v2 = neighbors_v2.end();
         it_v2 != to_v2; ++it_v2) {

      // check whether vertex was already used.
      // if (neighbors_v2_used.find(*it_v2) != neighbors_v2_used.)
      //  continue;

      if (utils::FfExpertGraph::match_neighbors(*it_v1, g1, *it_v2, g2,
                                                criterion, radius, v1, v2,
                                                my_radius + 1)) {
        neighbors_v2.erase(*it_v2);
        break;
      }
    }
  }
  bool result = neighbors_v2.empty();
  /* debug only
  std::cout << (result ? "match!" : "no match") << std::endl;
   */
  return result;
}

/**
 * find a subgraph in another graph
 */
std::vector<utils::Vertex> utils::FfExpertGraph::equal_subgraph(
    const utils::FfExpertGraph &g, const utils::Vertex::equality_enum criterion,
    size_t radius) const {
  // two graphs are equal when all nodes are physically equal and all
  // edges have the same weights

  std::vector<utils::Vertex> solutions;

  // first loop over the whole graph and try to find the first vertex
  for (std::vector<utils::Vertex>::const_iterator it = vertices().begin(),
                                                  to = vertices().end();
       it != to; ++it) {
    if (utils::FfExpertGraph::match_neighbors(g.vertices()[0], g, *it, *this,
                                              criterion, radius,
                                              g.vertices()[0], *it)) {
      solutions.push_back(*it);
    }
  }
  return solutions;
}

/**
 * create graph from a building block
 */
utils::FfExpertGraph::FfExpertGraph(const FfExpertGraphMapper &mapper,
                                    gcore::BbSolute const &bb) {
  for (int i = 0; i < bb.numAtoms(); ++i) {
    const gcore::AtomTopology &atom = bb.atom(i);
    vertices().push_back(utils::Vertex(i, mapper.map_mass(int(atom.mass())),
                                       atom.iac(), atom.charge(), atom.name(),
                                       bb.resName()));
  }

  gcore::BondIterator bi(bb);
  for (; bi; ++bi) {
    const int i = bi()[0];
    const int j = bi()[1];

    utils::Vertex *v_i = NULL;
    utils::Vertex *v_j = NULL;

    if (i < 0) {
      v_i = &vertices()[vertices().size() + i];
    } else if (i >= int(vertices().size())) {
      v_i = &vertices()[i - vertices().size()];
    } else {
      v_i = &vertices()[i];
    }

    if (j < 0) {
      v_j = &vertices()[vertices().size() + j];
    } else if (j >= int(vertices().size())) {
      v_j = &vertices()[j - vertices().size()];
    } else {
      v_j = &vertices()[j];
    }

    if (v_i != NULL && v_j != NULL)
      this->connect(*v_i, *v_j, mapper.map_bond(bi().type()));
  }
}

utils::FfExpertGraph::FfExpertGraph(const std::vector<std::string> &atom_name,
                                    const std::vector<std::string> &elements,
                                    const std::vector<gcore::Bond> &bond,
                                    const std::vector<int> &order) {
  assert(atom_name.size() == elements.size() &&
         atom_name.size() == order.size());
  // add vertices
  for (size_t i = 0; i < atom_name.size(); ++i) {
    vertices().push_back(
        utils::Vertex(i, utils::FfExpertGraphMapper::str2at(elements[order[i]]),
                      atom_name[order[i]]));
  }
  // add bonds
  for (const auto &b : bond) {
    // const unsigned int i = (*it)[0], j = (*it)[1];
    const unsigned int i = b[0];
    const unsigned int j = b[1];
    if (i >= vertices().size() || j >= vertices().size()) {
      throw gromos::Exception("FfExpertGraph", "invalid bond");
    }

    double bond_type;
    switch (b.type()) {
    case 1:
      bond_type = 1.0;
      break;
    case 2:
      bond_type = 2.0;
      break;
    case 4:
      bond_type = 1.5;
      break;
    case 3:
      bond_type = 3.0;
      break;
    default:
      throw gromos::Exception("FfExpertGraph",
                              "unkown bond type (1..4, 4 aromatic)");
    }
    connect(vertices()[i], vertices()[j], bond_type);
  }
}

utils::FfExpertGraphMapper::FfExpertGraphMapper(const std::string &filename) {
  gio::Ginstream file(filename);
  std::vector<std::string> buffer;
  file.getblock(buffer);
  if (buffer[0] != "FORCEFIELD")
    throw gromos::Exception(
        "MassElementMapper",
        "library file does not contain a FORCEFIELD block!");
  if (buffer[buffer.size() - 1].find("END") != 0)
    throw gromos::Exception("MassElementMapper", "No END in FORCEFIELD"
                                                 " block.");

  m_ff = buffer[1];

  file.getblock(buffer);
  if (buffer[0] != "ELEMENTMAPPING")
    throw gromos::Exception(
        "FfExpertGraphMapper",
        "library file does not contain a ELEMENTMAPPING block!");
  if (buffer[buffer.size() - 1].find("END") != 0)
    throw gromos::Exception("FfExpertGraphMapper", "No END in ELEMENTMAPPING"
                                                   " block.");
  // read in the lib
  for (size_t i = 1; i < buffer.size() - 1; i++) {
    int mass;
    std::string element;
    std::istringstream ss(buffer[i]);
    ss >> mass >> element;
    if (ss.fail())
      throw gromos::Exception("FfExpertGraphMapper",
                              "bad line in ELEMENTMAPPING"
                              " block.");

    m_mass_map[mass - 1] = str2at(element);
  }

  file.getblock(buffer);
  if (buffer[0] != "BONDMAPPING")
    throw gromos::Exception(
        "FfExpertGraphMapper",
        "library file does not contain a BONDMAPPING block!");
  if (buffer[buffer.size() - 1].find("END") != 0)
    throw gromos::Exception("FfExpertGraphMapper", "No END in BONDMAPPING"
                                                   " block.");
  // read in the lib
  for (size_t i = 1; i < buffer.size() - 1; i++) {
    int bond_type;
    double bond;
    std::istringstream ss(buffer[i]);
    ss >> bond_type >> bond;
    if (ss.fail())
      throw gromos::Exception("FfExpertGraphMapper",
                              "bad line in ELEMENTMAPPING"
                              " block.");
    m_bond_map[bond_type - 1] = bond;
  }
}

utils::Vertex::atom_type_t
utils::FfExpertGraphMapper::map_mass(const int &mass) const {
  std::map<int, utils::Vertex::atom_type_t>::const_iterator it =
      m_mass_map.find(mass);
  if (it == m_mass_map.end()) {
    std::ostringstream msg;
    msg.precision(8);
    msg << "Could not map mass: " << std::setw(15) << mass;
    throw gromos::Exception("FfExpertGraphMapper", msg.str());
  }

  return it->second;
}

double utils::FfExpertGraphMapper::map_bond(const int &bond) const {
  std::map<int, double>::const_iterator it = m_bond_map.find(bond);
  if (it == m_bond_map.end()) {
    std::ostringstream msg;
    msg.precision(8);
    msg << "Could not map bond: " << std::setw(15) << bond;
    throw gromos::Exception("MassElementMapper", msg.str());
  }

  return it->second;
}

utils::Vertex::atom_type_t
utils::FfExpertGraphMapper::str2at(const std::string &str) {
  if (str == "H")
    return utils::Vertex::at_H;
  if (str == "C")
    return utils::Vertex::at_C;
  if (str == "N")
    return utils::Vertex::at_N;
  if (str == "O")
    return utils::Vertex::at_O;
  if (str == "F")
    return utils::Vertex::at_F;
  if (str == "NA")
    return utils::Vertex::at_Na;
  if (str == "MG")
    return utils::Vertex::at_Mg;
  if (str == "SI")
    return utils::Vertex::at_Si;
  if (str == "P")
    return utils::Vertex::at_P;
  if (str == "S")
    return utils::Vertex::at_S;
  if (str == "CL")
    return utils::Vertex::at_Cl;
  if (str == "AR")
    return utils::Vertex::at_Ar;
  if (str == "CA")
    return utils::Vertex::at_Ca;
  if (str == "FE")
    return utils::Vertex::at_Fe;
  if (str == "CU")
    return utils::Vertex::at_Cu;
  if (str == "ZN")
    return utils::Vertex::at_Zn;
  if (str == "BR")
    return utils::Vertex::at_Br;
  throw gromos::Exception("MassElementMapper", "Don't know element: " + str);
}
