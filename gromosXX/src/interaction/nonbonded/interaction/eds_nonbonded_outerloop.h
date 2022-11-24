/**
 * @file eds_nonbonded_outerloop.h
 * the eds-perturbed non bonded outer loop:
 * Lennard-Jones and Coulomb interactions.
 */

#ifndef INCLUDED_EDS_NONBONDED_OUTERLOOP_H
#define INCLUDED_EDS_NONBONDED_OUTERLOOP_H

namespace interaction
{
  /**
   * @class Eds_Nonbonded_Outerloop
   * outerloop for the eds-perturbed nonbonded interactions.
   */
  class Eds_Nonbonded_Outerloop
  {
  public:    

    /**
     * Constructor.
     * @param nbp the nonbonded parameters
     */
    Eds_Nonbonded_Outerloop(Nonbonded_Parameter & nbp);

    /**
     * calculate the eds-perturbed interactions.
     */
    void eds_lj_crf_outerloop(topology::Topology & topo,
				    configuration::Configuration & conf,
				    simulation::Simulation & sim,
				    Pairlist const & pairlist,
				    Storage & storage);  

    /**
     * calculate the eds-perturbed 1,4-interactions.
     */
    void eds_one_four_outerloop(topology::Topology & topo,
				      configuration::Configuration & conf,
				      simulation::Simulation & sim,
				      Storage & storage);

    /**
     * calculate the eds-perturbed RF contributions for excluded atoms.
     */
    void eds_RF_excluded_outerloop(topology::Topology & topo,
					 configuration::Configuration & conf,
					 simulation::Simulation & sim,
					 Storage & storage);

    /**
     * calculate the perturbed self energy (polarisation).
     */
    void perturbed_self_energy_outerloop(topology::Topology & topo,
		            configuration::Configuration & conf,
		            simulation::Simulation & sim, 
			    Storage & storage);

    /**
     * calculate the perturbed electric field (polarisation).
     */
    void perturbed_electric_field_outerloop(topology::Topology & topo,
		            configuration::Configuration & conf,
		            simulation::Simulation & sim, 
		            PairlistContainer const & pairlist,
                            PairlistContainer const & perturbed_pairlist,
                            Storage & storage,
                            Storage & storage_lr,
                            int rank);

  private:
    Nonbonded_Parameter & m_param;

    /**
     * calculate the eds-perturbed interactions.
     */
    template<typename t_interaction_spec, typename t_perturbation_spec>
    void _eds_lj_crf_outerloop(topology::Topology & topo,
				     configuration::Configuration & conf,
				     simulation::Simulation & sim,
				     Pairlist const & pairlist,
				     Storage & storage);
       
    /**
     * calculate the eds-perturbed 1,4-interactions.
     */
    template<typename t_interaction_spec, typename t_perturbation_spec>
    void _eds_one_four_outerloop(topology::Topology & topo,
				       configuration::Configuration & conf,
				       simulation::Simulation & sim,
				       Storage & storage);

    /**
     * calculate the eds-perturbed RF contributions for excluded atoms.
     */
    template<typename t_interaction_spec, typename t_perturbation_spec>
    void _eds_RF_excluded_outerloop(topology::Topology & topo,
					  configuration::Configuration & conf,
					  simulation::Simulation & sim,
					  Storage & storage);

    template<typename t_interaction_spec, typename t_perturbation_spec>
    void _perturbed_self_energy_outerloop(topology::Topology & topo,
		            configuration::Configuration & conf,
		            simulation::Simulation & sim, 
			    Storage & storage);

    template<typename t_interaction_spec, typename t_perturbation_spec>
    void _perturbed_electric_field_outerloop(topology::Topology & topo,
		            configuration::Configuration & conf,
		            simulation::Simulation & sim, 
		            PairlistContainer const & pairlist,
                            PairlistContainer const & perturbed_pairlist,
			    Storage & storage,
                            Storage & storage_lr,
                            int rank);
  };
  
} // interaction

#endif
