/**
 * @file symmetry_restraint_interaction.h
 * symmetry restraining
 */

#ifndef INCLUDED_SYMMETRY_RESTRAINT_INTERACTION_H
#define INCLUDED_SYMMETRY_RESTRAINT_INTERACTION_H

namespace interaction
{
  /**
   * @class symmetry_restraint_interaction
   * calculates the symmetry restraining interaction
   */
  class Symmetry_Restraint_Interaction : public Interaction
  {
  public:
    /**
     * Constructor.
     */
    Symmetry_Restraint_Interaction() : Interaction("Symmetry Restraint") {}
    
    /**
     * Destructor.
     */
    virtual ~Symmetry_Restraint_Interaction() {}

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

  protected:
    
  };
  
} // interaction

#endif
