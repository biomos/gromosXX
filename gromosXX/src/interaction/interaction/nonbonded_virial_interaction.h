/**
 * @file nonbonded_virial_interaction.h
 * the nonbonded interactions with calculation of the
 * virial.
 */

#ifndef INCLUDED_NONBONDED_VIRIAL_INTERACTION_H
#define INCLUDED_NONBONDED_VIRIAL_INTERACTION_H

namespace interaction
{
  /**
   * @class Nonbonded_Virial_Interaction
   * calculate the nonbonded interactions
   * and the virial.
   */
  template<typename t_simulation, typename t_pairlist>
  class Nonbonded_Virial_Interaction : public Nonbonded_Interaction<t_simulation, t_pairlist>
  {
  public:
    /**
     * @enum virial_enum
     * type of virial to calculate.
     */
    enum virial_enum { atomic_virial = 0, molecular_virial = 1 };
    
    /**
     * Constructor.
     */
    Nonbonded_Virial_Interaction(virial_enum virial = molecular_virial);
    
    /**
     * Destructor.
     */
    virtual ~Nonbonded_Virial_Interaction();
    
    /**
     * calculate the interactions.
     */
    virtual void calculate_interactions(t_simulation &sim);
    
  protected:
    /**
     * calculate the normal interactions.
     * @param sim the simulation.
     * @param it start of the pairlist.
     * @param to end of the pairlist.
     * @param force the force array to store the forces in
     * (longrange for the longrange pairist, the normal system force
     * array for the shortrange pairlist).
     */
    virtual void do_interactions(t_simulation &sim,
				 typename t_pairlist::iterator it,
				 typename t_pairlist::iterator to,
				 typename Nonbonded_Interaction<t_simulation, t_pairlist>
				 ::nonbonded_type_enum range);
    
    /**
     * calculate the 1,4 interactions.
     */
    virtual void do_14_interactions(t_simulation &sim);
    
    /**
     * type of virial to calculate.
     */
    virial_enum m_virial_type;

    /**
     * the center of mass positions (for each atom).
     */
    math::VArray m_com_pos;

    /**
     * the total (molecular) kinetic energy
     */
    math::Vec m_tot_ekin;

    /**
     * the longrange virial.
     */
    math::Matrix m_longrange_virial;

  };
  
} // interaction

#include "nonbonded_virial_interaction.tcc"

#endif
