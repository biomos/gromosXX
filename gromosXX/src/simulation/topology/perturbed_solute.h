/**
 * @file perturbed_solute.h
 * the perturbed part of the solute.
 */

#ifndef INCLUDED_PERTURBED_SOLUTE_H
#define INCLUDED_PERTURBED_SOLUTE_H

namespace simulation
{
  /**
   * @class Perturbed_Solute
   * holds the information about perturbation of
   * the solute.
   */
  class Perturbed_Solute
  {
  public:
    /**
     * @struct perturbed_distance_constraint_struct
     * hold perturbed distance constraint information.
     */
    struct perturbed_distance_constraint_struct
      : public compound::distance_constraint_struct
    {
      /**
       * bond length in state A.
       */
      double A_b0;
      /*
       * bond length in state B.
       */
      double B_b0;
    };
    
    /**
     * const perturbed bonds.
     */
    std::vector<Perturbed_Bond> const & bonds()const;

    /**
     * perturbed bonds.
     */
    std::vector<Perturbed_Bond> & bonds();

    /**
     * const perturbed angles.
     */
    std::vector<Perturbed_Angle> const & angles()const;
    
    /**
     * perturbed angles.
     */
    std::vector<Perturbed_Angle> & angles();
    
    /**
     * perturbed atoms accessor.
     */
    std::map<size_t, Perturbed_Atom> & atoms();

    /**
     * const perturbed atoms accessor.
     */
    std::map<size_t, Perturbed_Atom> const & atoms()const;

    /**
     * const perturbed atom accessoe
     */
    Perturbed_Atom & atom(const size_t i);
    
    /**
     * perturbed atompairs.
     */
    std::vector<Perturbed_Atompair> & atompairs();
    
    /**
     * const perturbed atompairs.
     */
    std::vector<Perturbed_Atompair> const & atompairs()const;

    /**
     * perturbed distance constraints accessor.
     */
    std::vector<perturbed_distance_constraint_struct> & distance_constraints();
    /**
     * perturbed distance constraints const accessor.
     */
    std::vector<perturbed_distance_constraint_struct> const & 
    distance_constraints()const;

    /**
     * add a distance constraint.
     */
    void add_distance_constraint(int const i, int const j, double const A_b0, 
				 double const B_b0);
    
    /**
     * set the distance constraints according to lambda.
     */
    void set_distance_constraints(double const lambda);
    

  private:
    /**
     * the perturbed bonds.
     */
    std::vector<Perturbed_Bond> m_bond;

    /**
     * the perturbed angles.
     */
    std::vector<Perturbed_Angle> m_angle;
    
    /**
     * the perturbed atoms.
     */
    std::map<size_t, Perturbed_Atom> m_atom;
    
    /**
     * the perturbed atompairs.
     */
    std::vector<Perturbed_Atompair> m_atompair;
    
    /**
     * the perturbed distance constraints.
     */
    std::vector<perturbed_distance_constraint_struct> m_distance_constraint;

  };
  
} // simulation

#include "perturbed_solute.tcc"

#endif
