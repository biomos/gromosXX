/**
 * @file electric_field_interaction.h
 * electric field
 */

#ifndef INCLUDED_ELECTRIC_FIELD_INTERACTION_H
#define INCLUDED_ELECTRIC_FIELD_INTERACTION_H

namespace interaction
{
  /**
   * @class electric_field_interaction
   * calculates the forces on an applied electric field
   */
  class Electric_Field_Interaction : public Interaction
  {
  public:
    /**
     * Constructor.
     */
    Electric_Field_Interaction() : Interaction("ElectricField") {}

    /**
     * Destructor.
     */
    virtual ~Electric_Field_Interaction() {}

    /**
     * init
     */
    virtual int init(topology::Topology &topo,
		     configuration::Configuration &conf,
		     simulation::Simulation &sim,
		     std::ostream &os = std::cout,
		     bool quiet = false)
    {
      os << "Electric field interaction\n";
      return 0;
    }
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

