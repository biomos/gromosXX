/**
 * @file order_parameter_restraint_interaction.h
 * order parameter restraining
 */
#ifndef ORDER_PARAMETER_RESTRAINT_INTERACTION_H
#define	ORDER_PARAMETER_RESTRAINT_INTERACTION_H

namespace interaction
{
  /**
   * @class Order_Parameter_Restraint_Interaction
   * calculates the order parameter restraining interaction
   */
  class Order_Parameter_Restraint_Interaction :
    public Interaction
  {
  public:
    /**
     * Constructor.
     */
    Order_Parameter_Restraint_Interaction() : Interaction("OrderParameterRestraint") {}

    /**
     * Destructor.
     */
    virtual ~Order_Parameter_Restraint_Interaction() {}

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

  };

} // interaction


#endif	/* ORDER_PARAMETER_RESTRAINT_INTERACTION_H */

