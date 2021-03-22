/**
 * @file perturbed_distance_field_interaction.h
 * perturbed distance field restraining
 */

#ifndef INCLUDED_PERTURBED_DISTANCE_FIELD_INTERACTION_H
#define INCLUDED_PERTURBED_DISTANCE_FIELD_INTERACTION_H

namespace interaction
{
  /**
   * @class Perturbed_Distance_Field_Interaction
   * calculates the perturbed distance restraining interaction
   */
  class Perturbed_Distance_Field_Interaction : public Interaction
  {
  public:
    /**
     * Constructor.
     */
    Perturbed_Distance_Field_Interaction() : Interaction("PerturbedDistanceField"){};
    
    /**
     * Destructor.
     */
    virtual ~Perturbed_Distance_Field_Interaction() {
      DEBUG(2, "Perturbed_Distance_Field_Interaction: destructor");
    };

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
    util::Virtual_Atom va_i;
    util::Virtual_Atom va_j;

  };
  
} // interaction

#endif
