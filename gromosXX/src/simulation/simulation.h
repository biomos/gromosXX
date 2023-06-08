/**
 * @file simulation.h
 * the simulation class.
 */

#ifndef INCLUDED_SIMULATION_H
#define INCLUDED_SIMULATION_H

// necessary headers
#include "multibath.h"
#include "parameter.h"
#include "mpiControl.h"
#ifdef HAVE_HOOMD
#include <HOOMD_GROMOSXX_processor.h>
#endif

#include "cukernel/cuda_kernel.h"
#ifdef HAVE_LIBCUDART
#include "cukernel/cudaKernel.h"
#else
#define CUDA_KERNEL void
#endif



namespace simulation
{

  /**
   * @class Simulation
   * holds simulation data
   */
  class Simulation
  {
  public:

    /**
     * Constructor.
     */
    Simulation() : mpi(false), openmp(false), m_cuda_kernel(nullptr),
		   m_time_step_size(0),
		   m_steps(0), 
		   m_time(0) {}
    
    /**
     * the simulation parameter
     */
    simulation::Parameter & param() { return m_param;}
    /**
     * the simulation parameter as const
     */
    simulation::Parameter const & param() const { return m_param;}

    /**
     * multibath.
     */
    simulation::Multibath & multibath(){return m_param.multibath.multibath; }
    
    /**
     * multibath as const.
     */
    simulation::Multibath const & multibath()const{
      return m_param.multibath.multibath;
    }

     /**
     * MpiControl.
     */
    simulation::MpiControl & mpiControl(){
        return m_MpiControl;
    }
    /** 
     * MpiControl as const.
     */
    simulation::MpiControl const & mpiControl()const{
        return m_MpiControl; 
    }

     /**
     * CUDA_Kernel pointer mutator
     */
    void CUDA_Kernel(cudakernel::CUDA_Kernel * cuda_kernel) {
        m_cuda_kernel = cuda_kernel;
    }
    /** 
     * CUDA_kernel accessor
     */
    cudakernel::CUDA_Kernel * CUDA_Kernel() const{
        return m_cuda_kernel; 
    }

    
    /**
     * time step size
     */
    double & time_step_size() { return m_time_step_size; }

    /**
     * time step size
     */
    double time_step_size()const { return m_time_step_size; }

    /**
     * number of steps done.
     */
    unsigned int & steps() { return m_steps; }

    /**
     * number of steps done.
     */
    unsigned int steps()const { return m_steps; }

    /**
     * simulation time.
     */
    double & time() { return m_time; }

    /**
     * simulation time.
     */
    double time()const { return m_time; }

    /**
     * minimisation step size.
     */
    double & minimisation_step_size() { return m_minimisation_step_size; }

    /**
     * enable mpi?
     */
    bool mpi;

    /**
     * enable openmp?
     */
    bool openmp;

	/**
	 * Processor(s) choice to use for HOOMD code (CPU/GPUs)
	 */
#ifdef HAVE_HOOMD
	boost::shared_ptr<processor::Processor> proc;
#endif

  private:
    /**
     * the simulation parameters
     */
    Parameter m_param;

    /**
     * the MPI controller
     */
    MpiControl m_MpiControl;

    /**
     * the CUDA kernel
     */
    cudakernel::CUDA_Kernel * m_cuda_kernel;
    
    /**
     * the time step size
     */
    double m_time_step_size;

    /**
     * the number of steps done.
     */
    unsigned int m_steps;

    /**
     * the simulation time.
     */
    double m_time;

    /**
     * the minimisation step size.
     */
    double m_minimisation_step_size;

  };

} // simulation

#endif
