/**
 * @file compound.h
 * base class for solute and solvent.
 *
 * NOTE: differences between solute and solvent:
 * 
 * solute 0..1                    --    solvent 0..many (ie H2O DMSO mixture)
 *
 * NPM = 1                        --    NSM 0..many
 * (no multiplying of topology)         loop topology
 *
 * multiple molecules             --    single (rigid?) molecule
 *
 */

#ifndef INCLUDED_COMPOUND_H
#define INCLUDED_COMPOUND_H

namespace simulation
{
  /**
   * @class compound
   * common features of solute and solvent.
   */
  class compound
  {
  public:
    /**
     * @struct distance_constraint_struct
     * hold distance constraints.
     */
    struct distance_constraint_struct
    {
      int i;
      int j;
      double b0;
    };
    
    /**
     * Constructor.
     */
    explicit compound();
    
    /**
     * accessor to the distance constraints.
     */
    std::vector<distance_constraint_struct> & distance_constraints();
    /**
     * accessor to a single distance constraint.
     */
    distance_constraint_struct &distance_constraint(size_t i);
    /**
     * add a distance constraint.
     */
    void add_distance_constraint(int i, int j, double b0);
    /**
     * number of atoms in the compound.
     */
    size_t num_atoms()const;

  protected:
    std::vector<distance_constraint_struct> m_distance_constraint;
    size_t m_num_atoms;
  };
  
} // simulation

// inline methods
#include "compound.tcc"

#endif
