/**
 * @file nonbonded_outerloop.h
 * the non bonded outerloops.
 */

#include "latticesum.h"


#ifndef INCLUDED_NONBONDED_OUTERLOOP_H
#define INCLUDED_NONBONDED_OUTERLOOP_H

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
namespace util {
  class Algorithm_Timer;
}

namespace interaction
{
  class Nonbonded_Parameter;
  class Storage;
  class Pairlist;
  class Lattice_Sum;
  struct KSpace_Element;
  
  /**
   * @class Nonbonded_Outerloop
   * loops the nonbonded interactions...
   */
  class Nonbonded_Outerloop
  {
  public:    
    /**
     * Constructor.
     */
    Nonbonded_Outerloop(Nonbonded_Parameter & nbp);
    
    /**
     * calculate the lj crf interactions.
     */
    void lj_crf_outerloop(topology::Topology & topo,
			  configuration::Configuration & conf,
			  simulation::Simulation & sim,
			  Pairlist const & pairlist_solute,
                          Pairlist const & pairlist_solvent,
			  Storage & storage,
                          bool longrange = false);

    void cg_exclusions_outerloop(topology::Topology & topo,
				 configuration::Configuration & conf,
				 simulation::Simulation & sim,
				 Storage & storage);

    /**
     * calculate the 1,4-interactions.
     */
    void one_four_outerloop(topology::Topology & topo,
			    configuration::Configuration & conf,
			    simulation::Simulation & sim,
			    Storage & storage);

    /**
     * calculate the RF contributions for excluded atoms.
     */
    void RF_excluded_outerloop(topology::Topology & topo,
			       configuration::Configuration & conf,
			       simulation::Simulation & sim,
			       Storage & storage);
    
    /**
     * calculate the self energy (polarization).
     */
    void self_energy_outerloop(topology::Topology & topo,
		            configuration::Configuration & conf,
		            simulation::Simulation & sim, 
			    Storage & storage);

    /**
     * calculate the electric field (polarization).
     */
    void electric_field_outerloop(topology::Topology & topo,
		            configuration::Configuration & conf,
		            simulation::Simulation & sim, 
		            PairlistContainer const & pairlist,
                            Storage & storage,
                            Storage & storage_lr,
                            int rank);
    
    /**
     * calculate the ls interactions. (lattice sum)
     * in real space
     */
    void ls_real_outerloop(topology::Topology & topo,
            configuration::Configuration & conf,
            simulation::Simulation & sim,
            Pairlist const & pairlist_solute,
            Pairlist const & pairlist_solvent,
            Storage & storage, int rank, int size);
    
    /**
     * calculate the ls interactions kspace interactions
     * using the Ewald method
     */
    void ls_ewald_kspace_outerloop(topology::Topology & topo,
            configuration::Configuration & conf,
            simulation::Simulation & sim,
            Storage & storage, int rank, int size);
    
    /**
     * calculate the ls interactions kspace interactions
     * using the P3M method
     */
    void ls_p3m_kspace_outerloop(topology::Topology & topo,
            configuration::Configuration & conf,
            simulation::Simulation & sim,
            Storage & storage, int rank, int size,
            util::Algorithm_Timer & timer);
    
    /**
     * calculate the ls self interactions
     */
    void ls_self_outerloop(topology::Topology & topo,
            configuration::Configuration & conf,
            simulation::Simulation & sim,
            Storage & storage, int rank, int size);

    /** 
     * calculate the ls surface interaction
     */
    void ls_surface_outerloop(topology::Topology & topo,
            configuration::Configuration & conf,
            simulation::Simulation & sim,
            Storage & storage, int rank, int size);
    
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
     * calculate the hessian for a given pair
     */
    int calculate_hessian(topology::Topology & topo,
			  configuration::Configuration & conf,
			  simulation::Simulation & sim,
			  unsigned int atom_i, unsigned int atom_j,
			  math::Matrix & hessian,
			  PairlistContainer const & pairlist);

  private:
    /**
     * the nonbonded parameter.
     */
    Nonbonded_Parameter & m_param;

    template<typename t_interaction_spec>
    void _lj_crf_outerloop(topology::Topology & topo,
			   configuration::Configuration & conf,
			   simulation::Simulation & sim,
			   Pairlist const & pairlist_solute,
                           Pairlist const & pairlist_solvent,
			   Storage & storage,
                           bool longrange);

    template<typename t_interaction_spec>
    void _one_four_outerloop(topology::Topology & topo,
			     configuration::Configuration & conf,
			     simulation::Simulation & sim,
			     Storage & storage);

    template<typename t_interaction_spec>
    void _cg_exclusions_outerloop(topology::Topology & topo,
				  configuration::Configuration & conf,
				  simulation::Simulation & sim,
				  Storage & storage);

    template<typename t_interaction_spec>
    void _RF_excluded_outerloop(topology::Topology & topo,
				configuration::Configuration & conf,
				simulation::Simulation & sim,
				Storage & storage);
    
    template<typename t_interaction_spec>
    void _self_energy_outerloop(topology::Topology & topo,
		            configuration::Configuration & conf,
		            simulation::Simulation & sim, 
			    Storage & storage);

    template<typename t_interaction_spec>
    void _electric_field_outerloop(topology::Topology & topo,
            configuration::Configuration & conf,
            simulation::Simulation & sim,
            PairlistContainer const & pairlist,
            Storage & storage,
            Storage & storage_lr,
            int rank);
    
    template<typename t_interaction_spec>
    void _ls_real_outerloop(topology::Topology & topo,
            configuration::Configuration & conf,
            simulation::Simulation & sim,
            Pairlist const & pairlist_solute,
            Pairlist const & pairlist_solvent,
            Storage & storage, int rank, int size);
    
    template<typename t_interaction_spec>
    void _ls_ewald_kspace_outerloop(topology::Topology & topo,
            configuration::Configuration & conf,
            simulation::Simulation & sim,
            Storage & storage, int rank, int size);

    template<typename t_interaction_spec>
    void _ls_p3m_kspace_outerloop(topology::Topology & topo,
            configuration::Configuration & conf,
            simulation::Simulation & sim,
            Storage & storage, int rank, int size,
            util::Algorithm_Timer & timer);

    template<typename t_interaction_spec>
    void _ls_surface_outerloop(topology::Topology & topo,
            configuration::Configuration & conf,
            simulation::Simulation & sim,
            Storage & storage, int rank, int size);
    
    template<typename t_interaction_spec>
    int _calculate_interaction(topology::Topology & topo,
			       configuration::Configuration & conf,
			       simulation::Simulation & sim,
			       unsigned int atom_i, unsigned int atom_j,
			       math::Vec & force,
			       double & e_lj, double & e_crf);

    template<typename t_interaction_spec>
    int _calculate_hessian(topology::Topology & topo,
			   configuration::Configuration & conf,
			   simulation::Simulation & sim,
			   unsigned int atom_i, unsigned int atom_j,
			   math::Matrix & hessian,
			   PairlistContainer const & pairlist);
    
  };
  
} // interaction

#endif
