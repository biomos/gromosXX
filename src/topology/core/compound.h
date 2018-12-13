/**
 * @file compound.h
 * base class for solute and solvent.
 *
 * NOTE: differences between solute and solvent:
 * 
 * solute 0 / 1                   --    solvent 0..many (ie H2O DMSO mixture)
 *
 * NPM = 1                        --    NSM 0..many
 * (no multiplying of topology)         loop topology
 *
 * multiple molecules             --    single (rigid?) molecule
 *
 */

#ifndef INCLUDED_COMPOUND_H
#define INCLUDED_COMPOUND_H

namespace topology
{
  /**
   * @class Compound
   * common features of solute and solvent.
   */
  class Compound
  {
  public:
    /**
     * @struct lincs_struct
     * lincs constraints information.
     * if a topology is multiplied (or copied)
     * lincs has to be reinitialized!
     */
    struct lincs_struct
    {
      std::vector<std::vector<unsigned int> > coupled_constr;
      std::vector<std::vector<double> > coef;
      std::vector<double> sdiag;
    };
    
    /**
     * Constructor.
     */
    explicit Compound() : m_num_atoms(0){}
    
    /**
     * accessor to the distance constraints.
     */
    std::vector<two_body_term_struct> & distance_constraints(){return m_distance_constraint;}
    /**
     * const accessor to the distance constraints.
     */
    std::vector<two_body_term_struct> const & distance_constraints()const{return m_distance_constraint;}
    
    /**
     * accessor to a single distance constraint.
     */
    two_body_term_struct &distance_constraint(unsigned int i){return m_distance_constraint[i];}
  
    /**
     * add a distance constraint.
     */
    void add_distance_constraint(two_body_term_struct const s) {m_distance_constraint.push_back(s);}

    /**
     * accessor to the coupled distance constraints.
     */
    lincs_struct & lincs() 
    { return m_lincs;
    }

    /**
     * const accessor to the lincs struct
     */
    lincs_struct const & lincs() const
    { return m_lincs; 
    }

    /**
     * number of atoms in the compound.
     */
    unsigned int num_atoms()const {return m_num_atoms;}

  protected:
    /**
     * the distance constraints.
     */
    std::vector<two_body_term_struct> m_distance_constraint;

    /**
     * the lincs information
     */
    lincs_struct m_lincs;
    
    /**
     * the number of atoms in the compound.
     */
    unsigned int m_num_atoms;
  };
  
} // simulation

#endif
