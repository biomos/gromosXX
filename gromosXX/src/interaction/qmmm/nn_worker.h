/**
 * @file nn_worker.h
 * The worker class for the neural network interface (SchNetPack)
 * https://github.com/atomistic-machine-learning/schnetpack
 */
#ifndef INCLUDED_NN_WORKER_H
#define	INCLUDED_NN_WORKER_H

#ifdef HAVE_PYBIND11
  #include <pybind11/embed.h>
  namespace py = pybind11;
  #define PY_SCOPED_INTERPRETER   py::scoped_interpreter
  #define PY_MODULE_MAP           std::unordered_map<std::string, py::module>
  #define PY_OBJECT               py::object
#else
  #define DUMMY_TYPE              char
  #define PY_SCOPED_INTERPRETER   DUMMY_TYPE
  #define PY_MODULE_MAP           DUMMY_TYPE
  #define PY_OBJECT               DUMMY_TYPE
#endif

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
    virtual int init(const topology::Topology& topo
                   , const configuration::Configuration& conf
                   , simulation::Simulation& sim
                   , const interaction::QM_Zone& qm_zone) override;

  private:
    /**
     * Pointer to simulation parameters
     */
    simulation::Parameter::qmmm_struct::nn_param_struct* param;

    /**
     * PyBind interpreter scope guard
     */
    PY_SCOPED_INTERPRETER guard;

    /**
     * Python modules
     */
    PY_MODULE_MAP py_modules;

    /**
     * ASE calculator
     */
    PY_OBJECT ml_calculator;

    /**
     * ASE calculator for validation
     */
    PY_OBJECT val_calculator;

    /**
     * ASE calculator for charges
     */
    PY_OBJECT charge_calculator;

    /**
     * run the NN worker
     */
    int run_QM(topology::Topology& topo
                     , configuration::Configuration& conf
                     , simulation::Simulation& sim, interaction::QM_Zone & qm_zone) override;
  };
}

#endif	/* NN_WORKER_H */

