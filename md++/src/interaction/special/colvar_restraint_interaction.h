/**
 * @file colvar_restraint_interaction.h
 * collective variable restraining
 */

#ifndef INCLUDED_COLVAR_RESTRAINT_INTERACTION_H
#define INCLUDED_COLVAR_RESTRAINT_INTERACTION_H

#include "../../interaction/special/colvar/colvar.h"

namespace interaction
{
  /**
   * @class colvar_restraint_interaction
   * calculates the colvar restraining interaction
   */
  class Colvar_Restraint_Interaction : public Interaction
  {
  public:
    /**
     * Constructor.
     */
    Colvar_Restraint_Interaction() : Interaction("ColvarRestraint"), sum(false),
              Ctot0(0) {}
    
    /**
     * Destructor.
     */
    virtual ~Colvar_Restraint_Interaction();

    /**
     * init
     */
    virtual int init(topology::Topology &topo, 
		     configuration::Configuration &conf,
		     simulation::Simulation &sim,
		     std::ostream &os = std::cout,
		     bool quiet = false);
    /**
     * calculate the interactions.
     */
    virtual int calculate_interactions(topology::Topology & topo,
				       configuration::Configuration & conf,
				       simulation::Simulation & sim);
    
  private:
    /**
     * store pointers to all specified colvars
     */  
     std::vector<Colvar *> m_colvars;
     
    /**
     * sum of all target values
     */  
     double Ctot0;
     
    /** apply forces and calculate potential for given colvar
     * using the biasing function specified in the COLVARRES block
     */
     double apply_restraint(topology::Topology & topo,
      configuration::Configuration & conf, simulation::Simulation & sim, 
      std::vector< util::Virtual_Atom* > atoms, 
      math::VArray &derivatives, double &curr, double &target, double weight, int rah, double r_linear) ;
      
    /**
     * put the bias on the sum of all collective variables
     */  
     bool sum;
    
  };
  
} // interaction

#endif
