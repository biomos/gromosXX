/**
 * @file cuda_nonbonded.h
 * the non bonded interactions:
 * Lennard-Jones and Coulomb interactions.
 */
#include "nonbonded_interaction.h"

namespace interaction
{
  /**
   * @class CUDA_NonBonded
   * calculates the nonbonded interactions using CUDA
   */
  class CUDA_Nonbonded : 
    public Nonbonded_Interaction
  {
  public:    
    /**
     * Constructor.
     */
    CUDA_Nonbonded(Pairlist_Algorithm *pa);
    /**
     * Destructor.
     */
    virtual ~CUDA_Nonbonded();
    /**
     * calculate the interactions.
     */
    virtual int calculate_interactions(topology::Topology & topo,
				       configuration::Configuration & conf,
				       simulation::Simulation & sim);
    /**
     * size the arrays of storage.
     */
    virtual int init(topology::Topology & topo,
		     configuration::Configuration & conf,
		     simulation::Simulation & sim,
		     std::ostream & os = std::cout,
		     bool quiet = false);

  protected:

    /**
     * storage for stuff
     */
    Storage m_storage;

  };
  
} // interaction

