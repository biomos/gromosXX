/**
 * @file solvent.h
 * solvent information.
 * allows multiple solvents!
 */

#ifndef INCLUDED_SOLVENT_H
#define INCLUDED_SOLVENT_H

namespace simulation
{
  /**
   * @class solvent
   * holds solvent information
   * (for one solute).
   */
  class solvent
  {
  public:
    /**
     * @struct solventatom_struct
     * the solvent atom information.
     */
    struct solventatom_struct
    {
      std::string name;
      int iac;
      double mass;
      double charge;
    };
    
    struct solventconstraint_struct
    {
      int i;
      int j;
      double b0;
    };
    
    /**
     * accessor to the solvent atom information.
     */
    std::vector<solventatom_struct> & atoms();
    /**
     * accessor to a single atom of solvent atom struct.
     */
    solventatom_struct & atom(size_t i);
    /**
     * add a solventatom.
     */
    void add_atom(std::string name, int iac, double mass, double charge);
    /**
     * accessor to the solvent constraint information.
     */
    std::vector<solventconstraint_struct> & constraints();
    /**
     * accessor to a single constraint.
     */
    solventconstraint_struct &constraint(size_t i);
    /**
     * add a constraint.
     */
    void add_constraint(int i, int j, double b0);

  private:
    std::vector<solventatom_struct> m_atom;
    std::vector<solventconstraint_struct> m_constraint;

  };
  
} // simulation

// inline methods
#include "solvent.tcc"

#endif

  
    
