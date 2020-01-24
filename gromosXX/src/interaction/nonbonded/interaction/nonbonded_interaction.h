/**
 * @file nonbonded_interaction.h
 * the non bonded interactions:
 * Lennard-Jones and Coulomb interactions.
 */

#ifndef INCLUDED_NONBONDED_INTERACTION_H
#define INCLUDED_NONBONDED_INTERACTION_H

namespace topology
{
  class Topology;
}
namespace configuration
{
  class Configuration;
}
namespace simulation
{
  class Simulation;
}

#include "interaction.h"
#include "nonbonded_parameter.h"

namespace interaction
{

  class Pairlist_Algorithm;
  class Nonbonded_Set_Interface;
  
  /**
   * @class Nonbonded_Interaction
   * calculates the nonbonded interactions.
   */
  class Nonbonded_Interaction : 
    public Interaction
  {
  public:    

    /**
     * Constructor.
     */
    Nonbonded_Interaction(Pairlist_Algorithm *pa);
    
    /**
     * Destructor.
     */
    virtual ~Nonbonded_Interaction();

    /**
     * calculate the interactions.
     */
    virtual int calculate_interactions(topology::Topology & topo,
				       configuration::Configuration & conf,
				       simulation::Simulation & sim);

    /**
     * calculate the interaction for a given atom pair.
     * SLOW! as it has to create the periodicity...
     */
    int calculate_interaction(topology::Topology & topo,
			      configuration::Configuration & conf,
			      simulation::Simulation & sim,
			      unsigned int atom_i, unsigned int atom_j,
			      math::Vec & force, 
			      double &e_lj, double &e_crf);
    

    /**
     * calculate the hessian for a given atom.
     */
    int calculate_hessian(topology::Topology & topo,
			  configuration::Configuration & conf,
			  simulation::Simulation & sim,
			  unsigned int atom_i, unsigned int atom_j,
			  math::Matrix & hessian);

    /**
     * Pairlist algorithm
     */
    Pairlist_Algorithm & pairlist_algorithm()
    {
      return *m_pairlist_algorithm;
    }

    /**
     * change the pairlist algorithm
     */
    void pairlist_algorithm(Pairlist_Algorithm *pa)
    {
      m_pairlist_algorithm = pa;
    }
      
    /**
     * size the arrays of storage.
     */
    virtual int init(topology::Topology & topo,
		     configuration::Configuration & conf,
		     simulation::Simulation & sim,
		     std::ostream & os = std::cout,
		     bool quiet = false);
  
    /**
     * parameter
     */
    Nonbonded_Parameter & parameter() { return m_parameter; }

  protected:
    /**
     * check whether "fast" loops
     * should be used
     */
    int check_special_loop
    (
     topology::Topology const & topo,
     configuration::Configuration const & conf,
     simulation::Simulation & sim,
     std::ostream & os = std::cout,
     bool quiet = false);

    /**
     * initialize the expanded configuration
     */
    void init_expand_configuration(topology::Topology const & topo,
     configuration::Configuration & conf,
     simulation::Simulation & sim,
     configuration::Configuration & exp_conf);

    /**
     * expand a configuration for
     * multiple unit cell
     * simulations
     */
    void expand_configuration
    (
     topology::Topology const & topo,
     configuration::Configuration & conf,
     simulation::Simulation & sim,
     configuration::Configuration & exp_conf
     );

    template<math::boundary_enum b>
    void _expand_configuration
    (
     topology::Topology const & topo,
     configuration::Configuration & conf,
     simulation::Simulation & sim,
     configuration::Configuration & exp_conf
     );

    /**
     * reduce a configuration for
     * multiple unit cell
     * simulations
     */
    void reduce_configuration
    (
     topology::Topology const & topo,
     configuration::Configuration & conf,
     simulation::Simulation & sim,
     configuration::Configuration & exp_conf
     );
    
    /**
     * store data from sets into the configuration
     */
    void store_set_data
    (
     topology::Topology const & topo,
     configuration::Configuration & conf,
     simulation::Simulation const & sim
     );
    
    /**
     * print the pairlist
     */
    int print_pairlist
    (
     topology::Topology & topo,
     configuration::Configuration & conf,
     simulation::Simulation & sim,
     std::ostream & os = std::cout
     );

    /**
     * the pairlist update algorithm.
     */
    Pairlist_Algorithm * m_pairlist_algorithm;

    /**
     * a vector of nonbonded sets
     */
    std::vector<Nonbonded_Set_Interface *> m_nonbonded_set;

    /**
     * nonbonded parameter
     */
    Nonbonded_Parameter m_parameter;

    /**
     * number of sets to create
     */
    int m_set_size;

    /**
     * the expanded configuration;
     */
    configuration::Configuration * m_exp_conf;

  };
  
} // interaction

#endif
