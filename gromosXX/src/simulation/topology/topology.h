/**
 * @file topology.h
 * the topology class
 */

#ifndef INCLUDED_TOPOLOGY_H
#define INCLUDED_TOPOLOGY_H

namespace simulation
{
  /**
   * @class topology
   * holds the topological information of
   * the simulated system
   * @sa simulation::simulation
   * @sa simulation::system
   */
  class topology
  {
  public:
    /**
     * Constructor
     */
    explicit topology();

    /**
     * masses accessor
     */
    math::SArray &mass();

    /**
     * masses const accessor
     */
    math::SArray const & mass()const;

    /**
     * set the number of solute atoms
     */
    void num_solute_atoms(size_t atoms);
    
    /**
     * get the number of solute atoms
     */
    size_t num_solute_atoms()const;
    
  private:
    /**
     * the number of solute atoms.
     */
    size_t m_solute_atoms;
    
    /**
     * the atom masses.
     */
    math::SArray m_mass;
    
  }; // topology
  
} // simulation

// inline method definitions
#include "topology.tcc"

#endif
