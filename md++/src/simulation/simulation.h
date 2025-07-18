/*
 * This file is part of GROMOS.
 * 
 * Copyright (c) 2011, 2012, 2016, 2018, 2021, 2023 Biomos b.v.
 * See <https://www.gromos.net> for details.
 * 
 * GROMOS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

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
#include "gpu/cuda/manager/cuda_manager.h"
#ifdef HAVE_HOOMD
#include <HOOMD_GROMOSXX_processor.h>
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
    Simulation() : 
		   m_time_step_size(0),
		   m_steps(0), 
		   m_time(0),
       m_minimisation_step_size(0),
       m_mpi(false), m_openmp(false), m_cuda(false) {
#ifdef XXMPI
      m_mpi = true;
#endif
#ifdef OMP
      m_openmp = true;
#endif
#ifdef USE_CUDA
      m_cuda = true;
#endif
    }
    
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
     * CudaManager mutator
     */
    gpu::CudaManager & cuda() {
        return m_cuda_manager;
    }

    /**
     * CudaManager accessor
     */
    const gpu::CudaManager & cuda() const {
        return m_cuda_manager;
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
     * MPI flag
     */
    bool mpi_enabled() const { return m_mpi; }
    void set_mpi_enabled(bool enabled) {
#ifdef XXMPI
      // allow change only if OMP available
      m_mpi = enabled;
#else
      (void)enabled; // avoids compiler warnings
#endif
    }

    /**
     * OpenMP flag
     */
    bool openmp_enabled() const { return m_openmp; }
    void set_openmp_enabled(bool enabled) {
#ifdef OMP
      // allow change only if OMP available
      m_openmp = enabled;
#else
      (void)enabled; // avoids compiler warnings
#endif
    }

    /**
     * CUDA flag
     */
    bool cuda_enabled() const { return m_cuda; }
    void set_cuda_enabled(bool enabled) {
#ifdef USE_CUDA
      // allow change only if CUDA available
      m_cuda = enabled;
#else
      (void)enabled; // avoids compiler warnings
#endif
    }

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
     * the CUDA Manager
     */
    gpu::CudaManager m_cuda_manager;
    
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

    /**
     * enable mpi?
     */
    bool m_mpi;

    /**
     * enable openmp?
     */
    bool m_openmp;

    /**
     * enable cuda?
     */
    bool m_cuda;

  };

} // simulation

#endif
