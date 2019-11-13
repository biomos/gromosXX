/**
 * @file dftb_worker.h
 * The worker class for the DFTB+ software
 */ 

#ifndef INCLUDED_DFTB_WORKER_H
#define INCLUDED_DFTB_WORKER_H

#include "qm_worker.h"


namespace interaction {
  class QM_Worker;
  /**
   * @class dftb_Worker
   * a worker class which calls the dftb software
   */
  class DFTB_Worker : public QM_Worker {
  public:
    /**
     * Constructor
     */
    DFTB_Worker() : QM_Worker("DFTB Worker") {
        //this->get_new_qmID();
    }
    /**
     * Destructor
     */
    virtual ~DFTB_Worker();
    /**
     * initialise the QM worker
     * @return 0 if successful, non-zero on failure
     */
    virtual int init(topology::Topology & topo,
            configuration::Configuration & conf,
            simulation::Simulation & sim);
    /**
     * run a QM job in DFTB
     * @param qm_pos a vector containing the QM atom positions
     * @param mm_atoms the MM atoms to include
     * @param storage the energies, forces, charges obtained
     * @return 0 if successful, non-zero if not.
     */
    virtual int run_QM(const topology::Topology & topo,
                       const configuration::Configuration & conf,
                       const simulation::Simulation & sim,
                       interaction::QM_Zone & qm_zone);
  private:
    /**
     * file name for DFTB input file
     */
    std::string input_file;
    /**
     * file name for DFTB output file
     */
    std::string output_file;
    /**
     * file name for DFTB gradient output file
     */
    std::string output_gradient_file;
    /**
     * file name for DFTB charg file
     */
    std::string output_charg_file;
    /**
     * file name for DFTB geom file
     */
    std::string geom_file;
  };
}

#endif /* DFTB_WORKER_H */


