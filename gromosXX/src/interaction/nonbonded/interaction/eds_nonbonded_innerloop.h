/**
 * @file eds_nonbonded_innerloop.h
 * eds-perturbed inner loop class of the nonbonded routines.
 */

#ifndef INCLUDED_EDS_NONBONDED_INNERLOOP_H
#define INCLUDED_EDS_NONBONDED_INNERLOOP_H

namespace interaction
{
  /**
   * @class Eds_Nonbonded_Innerloop
   * eds-perturbed non bonded inner loop.
   */
  template<typename t_interaction_spec,
	   typename t_perturbation_details>
  class Eds_Nonbonded_Innerloop:
    public Eds_Nonbonded_Term
  {
  public:

    typedef math::Periodicity<t_interaction_spec::boundary_type> Periodicity_type;
    
    /**
     * Constructor
     */
    explicit Eds_Nonbonded_Innerloop(Nonbonded_Parameter &nbp): m_param(&nbp){}
    
    /**
     * eds-perturbed interaction
     */
    void eds_lj_crf_innerloop
    (
     topology::Topology & topo,
     configuration::Configuration & conf,
     unsigned int i, unsigned int j,
     Storage &storage,
     Periodicity_type const & periodicity
     );

    /**
     * eds-perturbed 1-4 interaction
     * (always shortrange)
     */
    void eds_one_four_interaction_innerloop
    (topology::Topology & topo, configuration::Configuration & conf,
     unsigned int i, unsigned int j,
     Periodicity_type const & periodicity);
    
    /**
     * eds-perturbed RF interaction (solute).
     * (always shortrange)
     */
    void eds_RF_excluded_interaction_innerloop
    ( topology::Topology & topo,
      configuration::Configuration & conf,
      std::map<unsigned int, topology::EDS_Perturbed_Atom>::const_iterator const & mit,
      Periodicity_type const & periodicity);
    
  protected:

    Nonbonded_Parameter * m_param;
    
  };
  
} // interaction

#include "eds_nonbonded_innerloop.cc"

#endif
