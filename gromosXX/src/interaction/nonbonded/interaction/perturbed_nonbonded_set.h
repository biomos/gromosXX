/**
 * @file perturbed_nonbonded_set.h
 * the (perturbed) non bonded interactions for a set of atoms:
 * Lennard-Jones and Coulomb interactions
 * perturbed and non-perturbed.
 */

#ifndef INCLUDED_PERTURBED_NONBONDED_SET_H
#define INCLUDED_PERTURBED_NONBONDED_SET_H

namespace interaction
{
  /**
   * @class Perturbed_Nonbonded_Set
   * calculates the (perturbed) nonbonded interactions.
   */
  class Perturbed_Nonbonded_Set : public Nonbonded_Set
  {
  public:    
    /**
     * Constructor.
     */
    Perturbed_Nonbonded_Set(Pairlist_Algorithm & pairlist_alg,
			    Nonbonded_Parameter & param);

    /**
     * Destructor
     */
    virtual ~Perturbed_Nonbonded_Set() {}
    
    /**
     * initialize some things
     */
    virtual int init(topology::Topology const & topo,
		     configuration::Configuration const & conf,
		     simulation::Simulation const & sim,
		     std::ostream & os = std::cout,
		     bool quiet = false);
    
    /**
     * calculate the interactions.
     */
    virtual int calculate_interactions(topology::Topology & topo,
				       configuration::Configuration & conf,
				       simulation::Simulation & sim,
				       int tid = 0, int num_threads = 1);

    /**
     * calculate the hessian for a given atom.
     */
    virtual int calculate_hessian(topology::Topology & topo,
				  configuration::Configuration & conf,
				  simulation::Simulation & sim,
				  unsigned int atom_i, unsigned int atom_j,
				  math::Matrix & hessian);
    
    Pairlist & perturbed_pairlist() { return m_perturbed_pairlist; }
    Pairlist const & perturbed_pairlist()const { return m_perturbed_pairlist; }

  protected:

    Pairlist m_perturbed_pairlist;

    Perturbed_Nonbonded_Outerloop m_perturbed_outerloop;

    Perturbed_Nonbonded_Pair m_perturbed_pair;

  };
  
} // interaction

#endif
