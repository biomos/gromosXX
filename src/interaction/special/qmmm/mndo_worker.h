/**
 * @file mndo_worker.h
 * The worker class for the MNDO QM software
 */
#ifndef MNDO_WORKER_H
#define	MNDO_WORKER_H

namespace interaction {
  class QM_Worker;
  /**
   * @class MNDO_Worker
   * a worker class which calls the MNDO software
   */
  class MNDO_Worker : public QM_Worker {
  public:
    /**
     * Constructor
     */
    MNDO_Worker() : QM_Worker("MNDO Worker") {}
    /**
     * Destructor
     */
    virtual ~MNDO_Worker();
    /**
     * initialise the QM worker
     * @return 0 if successful, non-zero on failure
     */
    virtual int init(topology::Topology & topo,
            configuration::Configuration & conf,
            simulation::Simulation & sim);
    /**
     * run a QM job in MNDO
     * @param qm_pos a vector containing the QM atom positions
     * @param mm_atoms the MM atoms to include
     * @param storage the energies, forces, charges obtained
     * @return 0 if successful, non-zero if not.
     */
    virtual int run_QM(topology::Topology & topo,
            configuration::Configuration & conf,
            simulation::Simulation & sim,
            const math::VArray & qm_pos,
            const std::vector<MM_Atom> & mm_atoms,
            interaction::QM_Storage & storage);
  private:
    /**
     * file name for MNDO input file
     */
    std::string input_file;
    /**
     * file name for MNDO output file
     */
    std::string output_file;
    /**
     * file name for MNDO gradient output file
     */
    std::string output_gradient_file;
  };
}

#endif	/* MNDO_WORKER_H */

