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
#include <util/cycle_thread.h>
#define gpu_status void


namespace interaction {

  /**
   * @class CUDA_Nonbonded_Set
   * calculates the nonbonded interactions.
   */
  class CUDA_Nonbonded_Set : public Nonbonded_Set, public util::CycleThread {
  public:

    /**
     * Constructor.
     */
    CUDA_Nonbonded_Set(Pairlist_Algorithm & pairlist_alg, Nonbonded_Parameter & param,
            int rank, int num_threads, 
            unsigned int gpu_id);

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

  private:
    /**
     * Pointer to the gpu status
     */
    gpu_status * gpu_stat;
    /**
     * Pointer to the topology
     */
    topology::Topology * mytopo;
    /**
     * Pointer to the configuration
     */
    configuration::Configuration * myconf;
    /**
     * Pointer to the simulation
     */
    simulation::Simulation * mysim;
    /**
     * Pointer to the stream
     */
    std::ostream * mystream;
    /**
     * Defines, wheter the output is verbouse or not
     */
    bool amIquiet;
    /**
     * Which physical gpu device
     */
    unsigned int gpu_device_number;
    /**
     * Which identifier the prgramm has for this gpu
     */
    unsigned int mygpu_id;     

    Nonbonded_Parameter * m_parameter;

    /**
     * Calculates the contstants needed for
     * further calculation
     */
    virtual void init_run();
    /**
     * contains the calculations, executed at every step
     */
    virtual void cycle();
    /**
     * Clean up
     */
    virtual void end_run();
    

  protected:
    Storage m_longrange_storage;
  };

} // interaction

#endif
