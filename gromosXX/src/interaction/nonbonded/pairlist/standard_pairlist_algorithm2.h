/**
 * @file standard_pairlist_algorithm2.h
 * create an atomic pairlist with a
 * chargegroup based cut-off criterion.
 */

#ifndef INCLUDED_STANDARD_PAIRLIST_ALGORITHM2_H
#define INCLUDED_STANDARD_PAIRLIST_ALGORITHM2_H

namespace interaction
{
  /**
   * @class Standard_Pairlist_Algorithm2
   * creates a pairlist.
   */
  template<typename t_interaction_spec, typename t_perturbation_spec>
  class Standard_Pairlist_Algorithm2 : 
    public Pairlist_Algorithm<t_interaction_spec, t_perturbation_spec>
  {
  public:
    typedef math::Periodicity<t_interaction_spec::boundary_type> Periodicity_type;

    /**
     * Constructor.
     */
    Standard_Pairlist_Algorithm2();
    /**
     * Destructor.
     */
    virtual ~Standard_Pairlist_Algorithm2(){}
    
    /**
     * prepare the pairlists
     */    
    virtual void prepare(topology::Topology & topo,
			 configuration::Configuration & conf,
			 simulation::Simulation &sim);

    /**
     * update the pairlist(s).
     */
    virtual void update(topology::Topology & topo,
			configuration::Configuration & conf,
			simulation::Simulation & sim,	
			Nonbonded_Set<t_interaction_spec, t_perturbation_spec> &nbs,
			size_t begin, size_t end, size_t stride);
        
  protected:

    /**
     * calculate the centres of geometry of the chargegroups
     * (if a chargegroup based cutoff scheme is applied)
     */
    void calc_cg_cog(topology::Topology & topo,
		configuration::Configuration & conf,
		simulation::Simulation & sim);

    /**
     * chargegroup center of geometry positions
     */
    math::VArray m_cog;

    /**
     * shortrange cutoff squared
     */
    double m_cutoff_short2;

    /**
     * longrange cutoff squared
     */
    double m_cutoff_long2;
    
  };
  
} // interaction


#include "standard_pairlist_algorithm2.cc"

#endif

  
