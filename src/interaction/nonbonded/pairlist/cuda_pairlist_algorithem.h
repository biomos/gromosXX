/* 
 * File:   CUDA_pairlist.h
 * Author: Axel && Mel
 *
 * Created on October 20, 2011, 11:34 AM
 */

#ifndef CUDA_PAIRLIST_ALGORITHEM_H
#define	CUDA_PAIRLIST_ALGORITHEM_H

namespace interaction
{

 class Pairlist;
 
  template<typename t_interaction_spec>
  class Nonbonded_Innerloop;

  template<typename t_interaction_spec, typename t_perturbation_details>
  class Perturbed_Nonbonded_Innerloop;

  /**
   * @class cuda_pairlist
   * create an atomic pairlist with a
   * chargegroup based or atom based
   *  cut-off criterion on a GPU.
   */


    class cuda_pairlist: public Pairlist_Algorithm
  {
  public:
      // constructor
      cuda_pairlist();

      // deconstructor
      virtual ~cuda_pairlist();

      // init
      virtual int init(topology::Topology &topo,
               configuration::Configuration &conf,
	       simulation::Simulation &sim,
               std::ostream &os = std::cout,
	       bool quiet = false);

      // prepare
      virtual int prepare(topology::Topology & topo,
                  configuration::Configuration & conf,
	          simulation::Simulation &sim); 
    
      // update
      virtual void update(topology::Topology & topo,
	          configuration::Configuration & conf,
		  simulation::Simulation & sim,
		  interaction::PairlistContainer & pairlist,
                  unsigned int begin,
                  unsigned int end,
	          unsigned int stride);

      virtual void update_atomic(topology::Topology & topo,
	          configuration::Configuration & conf,
		  simulation::Simulation & sim,
		  interaction::PairlistContainer & pairlist,
                  unsigned int begin,
                  unsigned int end,
	          unsigned int stride);

      virtual void update_cg(topology::Topology & topo,
	          configuration::Configuration & conf,
		  simulation::Simulation & sim,
		  interaction::PairlistContainer & pairlist,
                  unsigned int begin,
                  unsigned int end,
	          unsigned int stride);

    virtual void update_perturbed ( topology::Topology & topo,
                                    configuration::Configuration & conf,
                                    simulation::Simulation & sim,
                                    interaction::PairlistContainer & pairlist,
                                    interaction::PairlistContainer & perturbed_pairlist,
                                    unsigned int begin, unsigned int end,
                                    unsigned int stride);
    protected:
        void set_cutoff (double const cs, double const cl)
        {
            // set cutoff
                cutoff[0] = cs;
                cutoff[1] =  cs *  cs;
                cutoff[2] = cl;
                cutoff[3] =  cl *  cl;
        }
        bool excluded_solute_pair(topology::Topology & topo,
		       unsigned int i, unsigned int j);

#ifdef XXCUDA
        void prepare_exclusions(gcuda::TwoDArray<unsigned int> & earray, topology::Topology & topo);
#endif

    private:
        double m_solvent_solvent_timing; 
        double cutoff[4];
        bool atomic;

#ifdef XXCUDA
         // GPU memory map
         gcuda::memmap gpumemory;
#endif

    }; // cuda_pairlist


} // namespace interaction

#endif	/* CUDA_PAIRLIST_H */

