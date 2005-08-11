/**
 * @file pscale.h
 * periodic scaling
 */

#ifndef INCLUDED_PSCALE_H
#define INCLUDED_PSCALE_H

namespace interaction
{
  struct dihedral_type_struct;
  
  /**
   * @class Periodic_Scaling
   * implements periodic scaling (of topology parameters)
   */
  class Periodic_Scaling : public interaction::Interaction
  {
  public:
    /**
     * Constructor.
     */
    Periodic_Scaling(interaction::Forcefield & ff,
		     simulation::Parameter const & param);
    
    /**
     * Destructor.
     */
    virtual ~Periodic_Scaling(){}
    
    /**
     * Periodic Scaling algorithm
     */
    virtual int calculate_interactions(topology::Topology &topo, 
				       configuration::Configuration &conf,
				       simulation::Simulation &sim);

    /**
     * initialises data structures.
     * adds additional types for the dihedral potentials,
     * sets up a map between j-value restraint and dihedral angle.
     */
    virtual int init(topology::Topology &topo, 
		     configuration::Configuration &conf,
		     simulation::Simulation &sim,
		     std::ostream & os = std::cout,
		     bool quiet = false);
    
  private:
    /**
     * scale a force constant
     */
    double scale(double t, double T, double s);

    /**
     * Dihedral interaction
     */
    interaction::Dihedral_Interaction * m_DI;

  };
  
} // interaction

#endif

