/**
 * @file src/interaction/interaction.h
 * the interaction interface.
 */

#ifndef INCLUDED_INTERACTION_H
#define INCLUDED_INTERACTION_H

namespace interaction
{
  /**
   * @class Interaction
   * @interface Interaction
   * declares the interaction interface.
   */
  class Interaction
  {
  public:
    /**
     * Constructor.
     */
    Interaction(std::string name) : name(name), m_timing(0.0) {};
    /**
     * Destructor.
     */
    virtual ~Interaction(){};
    /**
     * the name of the interaction.
     * can be used to identify a special class.
     */
    std::string name;
    /**
     * calculate the interactions.
     */
    virtual int calculate_interactions(topology::Topology & topo,
				       configuration::Configuration & conf,
				       simulation::Simulation & sim) = 0;

    /**
     * calculate the second derivative of the potential
     * energy term with respect to the position of the
     * specified atom.
     */
    virtual int calculate_hessian(topology::Topology & topo,
				  configuration::Configuration & conf,
				  simulation::Simulation & sim,
				  size_t const atom_i, size_t const atom_j,
				  math::Matrix & hessian){ return 0; }
    
    /**
     * timing information.
     */
    virtual void print_timing(std::ostream & os)
    {
      os << "        "
	 << std::setw(36) << std::left << name
	 << std::setw(20) << m_timing << "\n";
    }

  protected:
    /**
     * store time used in algorithm.
     */
    double m_timing;

  };  
  
} // interaction

#endif
