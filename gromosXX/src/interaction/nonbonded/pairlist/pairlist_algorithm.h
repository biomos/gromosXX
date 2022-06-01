/**
 * @file pairlist_algorithm.h
 * base class for the pairlist(s) generating
 * algorithms.
 */

#ifndef INCLUDED_PAIRLIST_ALGORITHM_H
#define INCLUDED_PAIRLIST_ALGORITHM_H

#include <algorithm/algorithm.h>

namespace interaction
{
  class Pairlist;
  struct PairlistContainer;
  class Nonbonded_Parameter;
  
  /**
   * @class Pairlist_Algorithm
   * creates a pairlist.
   */
  class Pairlist_Algorithm : public algorithm::Algorithm
  {
  public:
    /**
     * Constructor.
     */
    Pairlist_Algorithm()
      : algorithm::Algorithm("PairlistAlgorithm"),
	m_param(NULL), p_timer(NULL)
    {}

    /**
     * destructor.
     */
    virtual ~Pairlist_Algorithm() {}

    void set_parameter(Nonbonded_Parameter * param) { m_param = param; }
    
    /**
     * prepare the pairlist(s).
     */
    virtual int prepare(topology::Topology & topo,
			 configuration::Configuration & conf,
			 simulation::Simulation &sim) = 0;
    /**
     * update the pairlist
     */
    virtual void update(topology::Topology & topo,
			configuration::Configuration & conf,
			simulation::Simulation &sim,
			interaction::PairlistContainer &pairlist,
			unsigned int begin, unsigned int end, 
			unsigned int stride) = 0;

    /**
     * update the pairlist, separating perturbed and non-perturbed interactions
     */
    virtual void update_perturbed(topology::Topology & topo,
				  configuration::Configuration & conf,
				  simulation::Simulation & sim,
                                  interaction::PairlistContainer & pairlist,
				  interaction::PairlistContainer & perturbed_pairlist,
				  unsigned int begin, unsigned int end, 
				  unsigned int stride) = 0;
    
    /**
     * print timing
     * no output... (done from forcefield)
     */
    virtual void print_timing(std::ostream & os) {}

    /**
     * apply the algorithm
     */
    virtual int apply(topology::Topology & topo,
		      configuration::Configuration & conf,
		      simulation::Simulation & sim) 
    {
      std::cerr << "don't call apply on a pairlist algorithm -- use update" << std::endl;
      assert(false);
      return 1;
    }
    
    /**
     * accessor to the timer pointer
     */
    void timer_pointer(util::Algorithm_Timer *t) {
      p_timer = t;
    }
    /**
     * accessor to the timer pointer
     */
    util::Algorithm_Timer * timer_pointer() {
      return p_timer;
    }
    
    /**
     * accessor to the timer
     */
    util::Algorithm_Timer & timer() {
      if (p_timer != NULL)
        return *p_timer;

      return m_timer;
    }
    /**
     * const accessor to the timer
     */
    const util::Algorithm_Timer & timer() const {
      if (p_timer != NULL)
        return *p_timer;

      return m_timer;
    }

    /**
     * Exclusion of QM atoms - for first atom of double-loop
     */
    static inline bool qm_excluded(const topology::Topology& topo
                                 , const simulation::qmmm_enum qmmm
                                 , const unsigned atom1) {
      /** 
       * With qmmm_mechanical, we exclude only QM-QM
       * With others we also exclude QM-MM
       */
      switch (qmmm) {
        case simulation::qmmm_off: break;
        case simulation::qmmm_mechanical: break;
        case simulation::qmmm_electrostatic:
        case simulation::qmmm_polarisable: {
          if (topo.is_qm(atom1))
            return true;
          break;
        }
        default:
          break;
      }
      return false;
    }

    /**
     * Exclusion of QM atoms - for second atom of double-loop or pair
     */
    static inline bool qm_excluded(const topology::Topology& topo
                                 , const simulation::qmmm_enum qmmm
                                 , const unsigned atom1
                                 , const unsigned atom2) {
      /** 
       * With qmmm_mechanical, we exclude only QM-QM
       * With others we also exclude QM-MM
       */
      switch (qmmm) {
        case simulation::qmmm_off: break;
        case simulation::qmmm_mechanical: {
          if ((topo.is_qm(atom1) &&
                (topo.is_qm(atom2) || topo.is_adaptive_qm_buffer(atom2)))
              ||
              (topo.is_qm(atom2) &&
                (topo.is_qm(atom1) || topo.is_adaptive_qm_buffer(atom1)))
          )
            return true;
          break;
        }
        case simulation::qmmm_electrostatic:
        case simulation::qmmm_polarisable: {
          if (topo.is_qm(atom1) || topo.is_qm(atom2))
            return true;
          break;
        }
        default: break;
      }
      return false;
    }

  protected:
    /**
     * nonbonded parameters (needed to construct the Innerloop).
     */
    Nonbonded_Parameter * m_param;
    /**
     * timer as pointer
     */
    util::Algorithm_Timer * p_timer;

  };

  class Failing_Pairlist_Algorithm : public Pairlist_Algorithm {
  public:
    Failing_Pairlist_Algorithm() : Pairlist_Algorithm(),
            failed(false), fallback_algorithm(NULL) {}
  protected:
    bool failed;
    Pairlist_Algorithm * fallback_algorithm;
  };
} // interaction

#endif
