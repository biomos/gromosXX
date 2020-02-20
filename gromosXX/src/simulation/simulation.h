/**
 * @file simulation.h
 * the simulation class.
 */

#ifndef INCLUDED_SIMULATION_H
#define INCLUDED_SIMULATION_H

// necessary headers
#include "multibath.h"
#include "parameter.h"
#ifdef HAVE_HOOMD
#include <HOOMD_GROMOSXX_processor.h>
#endif
#ifdef XXMPI
    #include <mpi.h>
#endif

namespace simulation
{
  /**
   * @class Simulation
   * holds simulation data
   */
     #ifdef XXMPI
        struct mpi_control_struct
        {
            /**
             * Constructor
             * Default values:
             * - number_of_threads 0
             */
            mpi_control_struct() : simulationID(0), numberOfThreads(-1), masterID(0), threadID(-1), MPIColor(1)
            {
            }

            int simulationID; //local replica id of simulation
            int numberOfThreads;    //total_number_of_threads      
            int masterID; //local master of this 
            int threadID;
            int MPIColor;
            MPI_Comm COMM; 

            std::vector<unsigned int> simulationOwnedThreads; 

        } /** replica exchange parameters */;
    #endif

  class Simulation
  {
  public:
    
    /**
     * Constructor.
     */
    Simulation() : mpi(false), openmp(false),
		   m_time_step_size(0),
		   m_steps(0), 
		   m_time(0) {
        #ifdef XXMPI
            mpi_control = mpi_control_struct();  //mpi test
        #endif
    }

    
    /**
     *  Constructor determening mpi flag on construction needed in repex/ 
     */
    Simulation(bool mpi) : mpi(mpi), openmp(false),
               m_time_step_size(0),
               m_steps(0), 
               m_time(0) {
            #ifdef XXMPI
                mpi_control = mpi_control_struct();  //mpi test
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
     * @struct mpi_control_struct
     *  ToDO:
     */
    #ifdef XXMPI
        struct mpi_control_struct mpi_control;  //mpi_control_struct test
    #endif

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
