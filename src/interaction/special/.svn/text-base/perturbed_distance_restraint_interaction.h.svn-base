/**
 * @file perturbed_distance_restraint_interaction.h
 * perturbed distance restraining
 */

#ifndef INCLUDED_PERTURBED_DISTANCE_RESTRAINT_INTERACTION_H
#define INCLUDED_PERTURBED_DISTANCE_RESTRAINT_INTERACTION_H

namespace interaction
{
  /**
   * @class Perturbed_Distance_Restraint_Interaction
   * calculates the perturbed distance restraining interaction
   */
  class Perturbed_Distance_Restraint_Interaction : public Interaction
  {
  public:
    /**
     * Constructor.
     */
    Perturbed_Distance_Restraint_Interaction() : Interaction("PerturbedDistanceRestraint"),
              exponential_term(0.0) {}
    
    /**
     * Destructor.
     */
    virtual ~Perturbed_Distance_Restraint_Interaction() {
      DEBUG(2, "Perturbed_Distance_Restraint_Interaction: destructor");
    }

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
    double exponential_term;
    /**
     * Determine the dimensions for which the distance restraint applies
     *
     * Distance restraints may be applied to selected dimensions only by specifying the appropriate
     * code for RAH in the DISTANCERESSPEC or PERTDISRESSPEC block: 
     * 
     * @verbatim
                                       value of RAH
                      --------------------------------------------
Dimensions to         Half harmonic   Full harmonic  Half harmonic
apply restraint       repulsive                      attractive

 x, y, z              -1                0             1
 x, y                  9               10            11
 x, z                 19               20            21
 y, z                 29               30            31
 x                    39               40            41
 y                    49               50            51
 z                    59               60            61
    * @endverbatim
    */
    std::map<int,math::Vec> rah_map;
    
  };
  
} // interaction

#endif
