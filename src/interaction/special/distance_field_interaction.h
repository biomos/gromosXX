/**
 * @file distance_field_interaction.h
 * distance field restraining
 */

#ifndef INCLUDED_DISTANCE_FIELD_INTERACTION_H
#define INCLUDED_DISTANCE_FIELD_INTERACTION_H

namespace interaction
{
  /**
   * @class distance_field_interaction
   * calculates the distance field restraining interaction
   */
  class Distance_Field_Interaction : public Interaction
  {
  public:
    /**
     * Constructor.
     */
    Distance_Field_Interaction() : Interaction("DistanceField"){}
    
    /**
     * Destructor.
     */
    virtual ~Distance_Field_Interaction() {}

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
