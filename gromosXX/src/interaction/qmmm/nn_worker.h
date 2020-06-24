/**
 * @file nn_worker.h
 * The worker class for the neural network interface
 */
#ifndef INCLUDED_NN_WORKER_H
#define	INCLUDED_NN_WORKER_H

#include <pybind11/embed.h>

namespace py = pybind11;

namespace interaction {
  /**
   * @class NN_Worker
   * a worker class which calls the NN interface
   */
  class __attribute__((visibility("default"))) NN_Worker : public QM_Worker {
  public:
    /**
     * Constructor
     */
    NN_Worker();
    /**
     * Destructor
     */
    virtual ~NN_Worker();
    /**
     * initialise the NN worker
     * @return 0 if successful, non-zero on failure
     */
    virtual int init(simulation::Simulation& sim);

  private:
    /**
     * Pointer to simulation parameters
     */
    simulation::Parameter::qmmm_struct::nn_param_struct* param;

    /**
     * Interface to neural network
     */
    //interaction::NN_Interface * nn_interface;

    /**
     * PyBind interpreter scope guard
     */
    py::scoped_interpreter guard;

    /**
     * Python modules
     */
    std::unordered_map<std::string, py::module> py_modules;

    /**
     * ASE calculator
     */
    py::object ml_calculator;

    /**
     * ASE calculator for validation
     */
    py::object val_calculator;

    /**
     * run the NN worker
     */
    int run_QM(topology::Topology& topo
                     , configuration::Configuration& conf
                     , simulation::Simulation& sim, interaction::QM_Zone & qm_zone);
    
    // int write_input(const topology::Topology& topo
    //                       , const configuration::Configuration& conf
    //                       , const simulation::Simulation& sim
    //                       , const interaction::QM_Zone & qm_zone) {return 0;};
    // /**
    //  * Open input file for QM
    //  */
    // int open_input(std::ofstream & inputfile_stream, const std::string & input_file) {return 0;};
    
    // /**
    //  * Call external QM program
    //  */
    // int system_call() {return 0;};
    // /**
    //  * read QM output files
    //  */
    // int read_output(topology::Topology& topo
    //                       , configuration::Configuration& conf
    //                       , simulation::Simulation& sim
    //                       , interaction::QM_Zone & qm_zone) {return 0;};
    // /**
    //  * Open QM output file
    //  */
    // int open_output(std::ifstream & outputfile_stream, const std::string & output_file) {return 0;};
  };
}

#endif	/* NN_WORKER_H */

