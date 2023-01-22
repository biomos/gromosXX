/**
 * @file qmmm_nonbonded_set.h
 * the parametric non bonded interactions between QM and MM atoms
 */

#ifndef INCLUDED_QMMM_NONBONDED_SET_H
#define INCLUDED_QMMM_NONBONDED_SET_H

namespace util
{
  class Algorithm_Timer;
}

namespace interaction
{
  class QMMM_Interaction;
  class QM_Zone;
  class Nonbonded_Parameter;
}

#include "../../../interaction/nonbonded/interaction/storage.h"
#include "../../../interaction/nonbonded/pairlist/pairlist.h"
#include "qmmm_nonbonded_outerloop.h"

namespace interaction
{
  /**
   * @class QMMM_Nonbonded_Set
   * calculates the nonbonded interactions within and around the QM zone
   */
  class QMMM_Nonbonded_Set
  {
  public:    
    /**
     * Constructor.
     */
    QMMM_Nonbonded_Set(QM_Zone& qm_zone
                     , util::Algorithm_Timer& timer
                     , Nonbonded_Parameter& param
		                 , int rank, int num_threads);
    
    /**
     * Destructor
     */
    virtual ~QMMM_Nonbonded_Set() {}
    
    /**
     * initialize some things
     */
    virtual int init(const topology::Topology& topo,
		                 const configuration::Configuration& conf,
		                 const simulation::Simulation& sim,
		                 std::ostream& os = std::cout,
		                 bool quiet = false);
    
    /**
     * calculate the interactions.
     */
    virtual int calculate_interactions(topology::Topology& topo,
				                               configuration::Configuration& conf,
				                               simulation::Simulation& sim);

    /**
     * store results in configuration
     */
    virtual int update_configuration(const topology::Topology& topo,
				                             configuration::Configuration& conf,
				                             const simulation::Simulation& sim);
    
    /**
     * calculate the hessian for a given atom.
     */
    virtual int calculate_hessian(topology::Topology & topo,
				  configuration::Configuration & conf,
				  simulation::Simulation & sim,
				  unsigned int atom_i, unsigned int atom_j,
				  math::Matrix & hessian);

    /**
     * pairlist accessor
     */
    PairlistContainer & pairlist() { return m_pairlist; }

    /**
     * pairlist mutator
     */
    PairlistContainer const & pairlist()const { return m_pairlist; }

    /**
     * shortrange storage accessor
     */
    Storage & storage()
    {
      return m_storage;
    }
    
    /**
     * timing function
     */
    void start_subtimer(std::string t) {
      if (m_rank == 0) m_timer.start_subtimer(t);
    }
    
    /**
     * timing function
     */
    void stop_subtimer(std::string t) {
      if (m_rank == 0) m_timer.stop_subtimer(t);
    }

  protected:
    /**
     * reference to QM zone
     */
    const QM_Zone& m_qm_zone;

    /**
     * reference to QMMM interaction timer
     */
    util::Algorithm_Timer& m_timer;

    /**
     * nonbonded outerloop
     */
    QMMM_Nonbonded_Outerloop m_outerloop;

    /**
     * set rank
     */
    const int m_rank;

    /**
     * number of threads
     */
    const int m_num_threads;

    /**
     * data storage
     */
    Storage m_storage;

    /**
     * pairlist
     */
    PairlistContainer m_pairlist;


  };
  
} // interaction

#endif
