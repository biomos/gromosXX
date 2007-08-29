/**
 * @file standard_pairlist_algorithm.h
 * standard pairlist algorithm (reference implementation)
 */

#ifndef INCLUDED_VGRID_PAIRLIST_ALGORITHM_H
#define INCLUDED_VGRID_PAIRLIST_ALGORITHM_H

namespace interaction
{
  class Storage;
  class Pairlist;
 
  template<typename t_interaction_spec>
  class Nonbonded_Innerloop;
  
  /**
   * @class VGrid_Pairlist_Algorithm
   * create a triple - range (atom - range) pairlist
   * (suited for vectorization)
   * using a grid based algorithm
   * and calculate long-range forces
   */
  class VGrid_Pairlist_Algorithm : 
    public Pairlist_Algorithm
  {
  public:
    /**
     * Constructor.
     */
    VGrid_Pairlist_Algorithm();

    /**
     * Destructor.
     */
    virtual ~VGrid_Pairlist_Algorithm() {}

    /**
     * init
     */
    virtual int init(topology::Topology &topo, 
		     configuration::Configuration &conf,
		     simulation::Simulation &sim,
		     std::ostream &os = std::cout,
		     bool quiet = false);
    /**
     * prepare the pairlists
     */    
    virtual int prepare
    (
     topology::Topology & topo,
     configuration::Configuration & conf,
     simulation::Simulation &sim
     );
    
    /**
     * update the pairlist
     */
    virtual void update
    (
     topology::Topology & topo,
     configuration::Configuration & conf,
     simulation::Simulation & sim,
     interaction::Storage & storage,
     interaction::Pairlist & pairlist,
     unsigned int begin, unsigned int end,
     unsigned int stride
     );
    
    /**
     * update the pairlist, separating perturbed and nonperturbed interactions
     */
    virtual void update_perturbed
    (
     topology::Topology & topo,
     configuration::Configuration & conf,
     simulation::Simulation & sim,
     interaction::Storage & storage,
     interaction::Pairlist & pairlist,
     interaction::Pairlist & perturbed_pairlist,
     unsigned int begin, unsigned int end,
     unsigned int stride
     );

  protected:
    
    /**
     * timing information.
     */
    virtual void print_timing(std::ostream & os)
    {
      os << "            "
	 << std::setw(32) << std::left << "solv - solv"
	 << std::setw(20) << m_solvent_solvent_timing << "\n"
	 << "            "
	 << std::setw(32) << std::left << "spc - spc"
	 << std::setw(20) << m_spc_timing << "\n";
    }
      
  private:
    /**
     * solvent - solvent pairlist + longrange
     */
    double m_solvent_solvent_timing;
    /**
     * spc specialized solvent - solvent pairlist + longrange
     */
    double m_spc_timing;

  };
} // interaction

#endif
