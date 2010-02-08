/**
 * @file cuda_nonbonded_set.h
 * the non bonded interactions for a set of atoms:
 * Lennard-Jones and Coulomb interactions

 */

#ifndef INCLUDED_CUDA_NONBONDED_SET_H
#define INCLUDED_CUDA_NONBONDED_SET_H

#include "pairlist.h"
#include "storage.h"
#include "nonbonded_outerloop.h"
#include "nonbonded_set_interface.h"
#include "nonbonded_set.h"
#define gpu_status void

namespace interaction {

  /**
   * @class CUDA_Nonbonded_Set
   * calculates the nonbonded interactions.
   */
  class CUDA_Nonbonded_Set : public Nonbonded_Set {
  public:
    /**
     * Constructor.
     */
    CUDA_Nonbonded_Set(Pairlist_Algorithm & pairlist_alg, Nonbonded_Parameter & param,
            int rank, int num_threads);

    /**
     * Destructor
     */
    virtual ~CUDA_Nonbonded_Set();

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
            simulation::Simulation & sim);

    Storage & shortrange_storage() {
      return m_storage;
    }

    Storage & longrange_storage() {
      return m_longrange_storage;
    }

    /**
     * maximal number of neighbors per atom in shortrange pairlist
     */
    unsigned int estNeigh_short;
    /**
     * maximal number of neighbors per atom in longrange pairlist
     */
    unsigned int estNeigh_long;
    /**
     * Set the gpu_status pointer
     */
    void pre_init(gpu_status * gpu_s);

  private:
    gpu_status * gpu_stat;

  protected:
    Storage m_longrange_storage;
    Storage m_longrange_storage_cuda;

    Nonbonded_Outerloop m_outerloop;
  };

} // interaction

#endif
