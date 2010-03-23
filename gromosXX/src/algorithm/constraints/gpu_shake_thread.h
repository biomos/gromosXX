/**
 * @file gpu_shake_thread.h
 * Header file, which implements pthread for M_SHAKE on the GPU
 * Similar to cudaNonbondedSet
 */

#ifndef _GPU_SHAKE_THREAD_H
#define	_GPU_SHAKE_THREAD_H

#include <util/cycle_thread.h>
#ifdef HAVE_LIBCUKERNEL
#include <cudaKernel.h>
#else
#define gpu_status void
#endif
#include <iostream>

#include <math/gmath.h>

namespace interaction {
  struct bond_type_struct;
}

namespace algorithm {

  /**
   * @class GPU_Shake_Thread
   * Uses CycleThread to operate the GPUs
   */
  class GPU_Shake_Thread : public util::CycleThread {
  public:

    /**
     * Constructor
     */
    GPU_Shake_Thread(topology::Topology & topo,
            configuration::Configuration & conf,
            simulation::Simulation & sim,
            math::GenericMatrix<double> factor,
            math::Vec mass,
            math::Vec constr_length2,
            double tolerance,
            unsigned int num_gpus,
            unsigned int gpu_id,
            std::ostream & os = std::cout,
            bool quiet = false);
    /**
     * Destructor
     */
    ~GPU_Shake_Thread();
    /**
     * Apply
     */
    int apply();

  private:
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
     *  factor matrix
     */
    math::GenericMatrix<double> factor;
    /**
     * pointer to the mass
     */
    math::Vec mass;
    /**
     * constraint length
     */
    math::Vec constr_length2;
    /**
     * tolerance
     */
    double m_tolerance;
    /**
     * Number of GPUs
     */
    unsigned int num_gpus;
    /**
     * Which GPU?
     */
    unsigned int gpu_id;
    /**
     * Pointer to the stream
     */
    std::ostream * mystream;
    /**
     * Defines, wheter the output is verbouse or not
     */
    bool amIquiet;
    /**
     * The index of the molecule, which causes m_shake to fail
     */
    int shake_fail_mol;

  };
}


#endif	/* _GPU_SHAKE_THREAD_H */

