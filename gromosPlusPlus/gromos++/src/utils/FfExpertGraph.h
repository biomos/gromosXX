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
 * @file FfExpertGraph.h
 * holds a graph class for sub structure matching
 */

#ifndef INCLUDED_FFEXPERTGRAPH_H
#define	INCLUDED_FFEXPERTGRAPH_H

#include <string>
#include <map>

#include "../gcore/BbSolute.h"

namespace gcore
{
  class BuildingBlock;
  class Bond;
  class Angle;
  class Improper;
  class Dihedral;
}

namespace utils {
  class FfExpertGraphMapper;
  
  /**
   * @class Vertex
   * a vertex in the graph. Corresponds to a (solute) atom.
   */
  class Vertex {
  public:
    /**
     * @enum equality_enum
     * an enum for the equality criterion
     */
    enum equality_enum {
      equal_type, equal_iac, equal_charge, equal_name
    };

    /**
     * @enum atom_type_t
     * atom types (elements)
     */
    enum atom_type_t {
      at_H, at_C, at_N, at_O, at_F, at_Na, at_Mg, at_Si, at_P, at_S, at_Cl,
      at_Ar, at_Ca, at_Fe, at_Cu, at_Zn, at_Br
    };

    /**
     * constructor
     */
    Vertex() : id(0), type(at_H), iac(0), charge(0.0), name(""), residue("") {
    }
    
    /**
     * constructor with some arguments
     */
    Vertex(int i, atom_type_t t, const std::string & n) :
    id(i), type(t), iac(0), charge(0.0), name(n), residue("") {}

    /**
     * constructor with arguments
     */
    Vertex(int i, atom_type_t t, int iac, double charge, const std::string & n,
            const std::string & rn) : id(i), type(t), iac(iac), charge(charge),
    name(n), residue(rn) {
    };
    /**
     * the ID of the vertex
     */
    int id;
    /**
     * the element type of the vertex
     */
    atom_type_t type;
    /**
     * Integer atom code
     */
    int iac;
    /**
     * charge
     */
    double charge;
    /**
     * label of the vertex
     */
    std::string name;
    /**
     * residue to which the vertex belongs (for tracking)
     */
    std::string residue;

    /**
     * copy constructor. Does not copy neighbors!
     */
    Vertex(const Vertex & v) : id(v.id), type(v.type), iac(v.iac),
    charge(v.charge), name(v.name), residue(v.residue) {
    }

    /**
     * assignment operator
     */
    Vertex & operator=(const Vertex &v) {
      id = v.id;
      type = v.type;
      iac = v.iac, charge = v.charge;
      name = v.name;
      residue = v.residue;
      return *this;
    }

    /**
     * less than operator
     */
    bool operator<(const Vertex &v) const {
      return id < v.id;
    }

    /**
     * equality operator
     */
    bool operator==(const Vertex &v) const {
      return id == v.id;
    }
    
    /**
     * inequality operator
     */
    bool operator!=(const Vertex &v) const {
      return id != v.id;
    }

    /**
     * physical equality operator
     */
    bool equals(const Vertex &v, equality_enum crit) const {
      switch (crit) {
        case equal_type:
          return type == v.type;
        case equal_iac:
          return iac == v.iac;
        case equal_charge:
          return charge == v.charge;
        case equal_name:
          return name == v.name;
        default: // type is default
          return type == v.type;
      }
    }
    
    /**
     * write it for debug
     */
    friend std::ostream & operator<<(std::ostream & os, const Vertex & v);
  };

  /**
   * @class FfExpertGraph
   * a graph class to do substructure matching
   */
  class FfExpertGraph {
  public:
    
    /**
     * default constructor
     */
    FfExpertGraph() {}
    
    /**
     * create a graph from a GROMOS BuildingBlock
     */
    FfExpertGraph(const FfExpertGraphMapper & mapper, gcore::BbSolute const & bb);
    
    /**
     * create a graph from PDB file information
     */
    FfExpertGraph(const std::vector<std::string> & atom_name,
            const std::vector<std::string> & elements,
            const std::vector<gcore::Bond> & bond,
            const std::vector<int> & order);
    
    typedef std::pair<Vertex, Vertex> vertex_pair_t;
    typedef std::map<std::pair<Vertex, Vertex>, double> weight_matrix_t;

    /**
     * accessor to vertices
     */
    std::vector<Vertex> & vertices() {
      return m_vertices;
    }

    /**
     * const accessor to vertices
     */
    const std::vector<Vertex> & vertices() const {
      return m_vertices;
    }

    /**
     * get the index of a vertex pointer
     */
    unsigned int vertex_index(const Vertex & v) const {
        unsigned int a;
      for (unsigned int i = 0; i < m_vertices.size(); ++i)
        if (m_vertices[i] == v) {
            a=i;
        }
      return a;
    }

    /**
     * connect two vertices
     */
    void connect(const Vertex & v1, const Vertex & v2, const double & weight) {
      m_weights[vertex_pair_t(v1, v2)] = weight;
      m_weights[vertex_pair_t(v2, v1)] = weight;
    }

    /**
     * get the weight or 0 if not connected
     */
    double weight(const Vertex & v1, const Vertex & v2) const {
      vertex_pair_t p(v1, v2);
      weight_matrix_t::const_iterator it = m_weights.find(p);
      if (it == m_weights.end())
        return 0.0;
      else
        return it->second;
    }

    /**
     * check whether two vertices are connected
     */
    double connected(const Vertex & v1, const Vertex & v2) const {
      return weight(v1, v2) != 0;
    }
    /**
     * output in graphviz dot format
     */
    friend std::ostream & operator<<(std::ostream &os, const FfExpertGraph & g);
    /**
     * output in graphviz dot format
     */
    void print(std::ostream &os, const Vertex::equality_enum criterion = Vertex::equal_name) const;
    
    /**
     * get the neighbors of a vertex
     */
    std::set<Vertex> get_neighbors(const Vertex &v) const;
    /**
     * cut a subgraph from a vertex and a given radius
     */
    FfExpertGraph cut_subgraph(const Vertex &centre, size_t radius) const;
    /**
     * find a subgraph in the graph and return the matching vertices
     */
    std::vector<Vertex> equal_subgraph(const FfExpertGraph & g,
            const Vertex::equality_enum criterion, size_t radius) const;
  protected:
    /**
     * vector containing the vertices
     */
    std::vector<Vertex> m_vertices;
    /**
     * map containing the weights
     */
    weight_matrix_t m_weights;
    /**
     * match two vertices of two different graphs
     */
    static bool match_neighbors(const Vertex & v1, const FfExpertGraph & g1,
            const Vertex & v2, const FfExpertGraph & g2,
            const Vertex::equality_enum criterion, const unsigned int radius, 
            const utils::Vertex & from_v1, const utils::Vertex & from_v2,
            const unsigned int my_radius = 0
            );
  };
  
  /**
   * @class FfExpertGraphMapper
   * maps masses to elements, bonds types to bonds...
   */
  class FfExpertGraphMapper {
  public:
    static Vertex::atom_type_t str2at(const std::string & str);
    /**
     * constructor. Takes a file containing the mapping library
     */
    FfExpertGraphMapper(const std::string & filename);
    
    /**
     * map a mass
     */
    Vertex::atom_type_t map_mass(const int & mass) const;
    
    /**
     * map a bond
     */
    double map_bond(const int & bond) const;
    
    /** 
     * const accessor to force field
     */
    const std::string & ForceField() const {
      return m_ff;
    }
    /** 
     * accessor to force field
     */
    std::string & ForceField() {
      return m_ff;
    }
    
  protected:
    std::map<int, Vertex::atom_type_t> m_mass_map;
    std::map<int, double> m_bond_map;
    std::string m_ff;
    /**
     * get atom type from string
     */
  };
} // namespace

#endif	/* INDLUCDED_FFEXPERTGRAPH_H */

