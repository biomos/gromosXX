/**
 * @file standard_pairlist_algorithm.h
 * standard pairlist algorithm (reference implementation)
 */

#ifndef INCLUDED_STANDARD_PAIRLIST_ALGORITHM_H
#define INCLUDED_STANDARD_PAIRLIST_ALGORITHM_H

namespace math
{
  template<math::boundary_enum>
  class Periodicity;
}

namespace interaction
{
  class Pairlist;
 
  
  template<typename t_interaction_spec>
  class Nonbonded_Innerloop;
  
  template<typename t_interaction_spec, typename t_perturbation_details>
  class Perturbed_Nonbonded_Innerloop; 
  
  /**
   * @class Standard_Pairlist_Algorithm
   * create an atomic pairlist with a
   * chargegroup based or atom based
   *  cut-off criterion.
   */
  class Standard_Pairlist_Algorithm : 
    public Pairlist_Algorithm
  {
  public:
    /**
     * Constructor.
     */
    Standard_Pairlist_Algorithm();

    /**
     * Destructor.
     */
    virtual ~Standard_Pairlist_Algorithm() {}

    /**
     * init
     */
    virtual int init(topology::Topology &topo, 
		     configuration::Configuration &conf,
		     simulation::Simulation &sim,
		     std::ostream &os = std::cout,
		     bool quiet = false) 
    {
      if (!quiet)
       os << "\tstandard pairlist algorithm\n";
      return 0;
    };

    /**
     * prepare the pairlists
     */    
    virtual int prepare(topology::Topology & topo,
			 configuration::Configuration & conf,
			 simulation::Simulation &sim);

    /**
     * update the pairlist
     */
    virtual void update(topology::Topology & topo,
			configuration::Configuration & conf,
			simulation::Simulation & sim,
			interaction::PairlistContainer & pairlist,
			unsigned int begin, unsigned int end,
			unsigned int stride);

    /**
     * update the pairlist, separating perturbed and nonperturbed interactions
     */
    virtual void update_perturbed
    (
     topology::Topology & topo,
     configuration::Configuration & conf,
     simulation::Simulation & sim,
     interaction::PairlistContainer & pairlist,
     interaction::PairlistContainer & perturbed_pairlist,
     unsigned int begin, unsigned int end,
     unsigned int stride
     );

    void update_perturbed_atomic
    (
     topology::Topology & topo,
     configuration::Configuration & conf,
     simulation::Simulation & sim,
     interaction::PairlistContainer & pairlist,
     interaction::PairlistContainer & perturbed_pairlist,
     unsigned int begin, unsigned int end,
     unsigned int stride
     );
        
  protected:

    void update_cg(topology::Topology & topo,
		   configuration::Configuration & conf,
		   simulation::Simulation & sim,
		   interaction::PairlistContainer & pairlist,
		   unsigned int begin, unsigned int end,
		   unsigned int stride);

    //template<typename t_interaction_spec>
    template<math::boundary_enum b>
    void _update_cg(topology::Topology & topo,
		    configuration::Configuration & conf,
		    simulation::Simulation & sim, 
		    interaction::PairlistContainer & pairlist,
		    unsigned int begin, unsigned int end,
		    unsigned int stride);

    template<math::boundary_enum b>
    void _update_pert_cg(topology::Topology & topo,
			 configuration::Configuration & conf,
			 simulation::Simulation & sim, 
			 interaction::PairlistContainer & pairlist,
			 interaction::PairlistContainer & perturbed_pairlist,
			 unsigned int begin, unsigned int end,
			 unsigned int stride);
    
    
    template<math::boundary_enum b>
    void _solvent_solvent
    (
     topology::Topology & topo,
     configuration::Configuration & conf,
     simulation::Simulation & sim,
     interaction::PairlistContainer & pairlist,
     int cg1, int stride,
     math::Periodicity<b> const & periodicity
     );

    void update_atomic(topology::Topology & topo,
		       configuration::Configuration & conf,
		       simulation::Simulation & sim,
		       interaction::PairlistContainer & pairlist,
		       unsigned int begin, unsigned int end,
		       unsigned int stride);

    template<math::boundary_enum b>
    void _update_atomic(topology::Topology & topo,
			configuration::Configuration & conf,
			simulation::Simulation & sim, 
			interaction::PairlistContainer & pairlist,
			unsigned int begin, unsigned int end,
			unsigned int stride);
    
    template<math::boundary_enum b>
    void _update_pert_atomic(topology::Topology & topo,
			     configuration::Configuration & conf,
			     simulation::Simulation & sim, 
			     interaction::PairlistContainer & pairlist,
			     interaction::PairlistContainer & perturbed_pairlist,
			     unsigned int begin, unsigned int end,
			     unsigned int stride);

    bool excluded_solute_pair(topology::Topology & topo,
			      unsigned int i, unsigned int j);

    void set_cutoff(double const cutoff_short, double const cutoff_long)
    {
      m_cutoff_long = cutoff_long;
      m_cutoff_short = cutoff_short;
      m_cutoff_short_2 = cutoff_short * cutoff_short;
      m_cutoff_long_2  = cutoff_long * cutoff_long;
    }
    
    bool insert_pair
    (
     topology::Topology & topo,
     interaction::Pairlist & pairlist,
     interaction::Pairlist & perturbed_pairlist,
     int a1, int a2,
     bool scaled_only
     );

    /**
     * timing information.
     */
    virtual void print_timing(std::ostream & os)
    {
      os << "            "
	 << std::setw(32) << std::left << "solv - solv pairlist"
	 << std::setw(20) << m_solvent_solvent_timing << "\n";
    }
      
  private:
    /**
     * the chargegroup center of geometries.
     */
    math::VArray m_cg_cog;
    /**
     * squared shortrange cutoff.
     */
    double m_cutoff_short_2;
    /**
     * squared longrange cutoff.
     */
    double m_cutoff_long_2;
    /**
     * longrange cutoff.
     */
    double m_cutoff_long;
    /**
     * shortrange cutoff.
     */
    double m_cutoff_short;
    /**
     * solvent - solvent pairlist 
     */
    double m_solvent_solvent_timing;
    
    math::Vec m_half_box;
    math::Vec m_box;

  };
} // interaction

#endif
