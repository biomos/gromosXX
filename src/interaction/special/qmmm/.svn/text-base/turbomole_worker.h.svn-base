/**
 * @file turbomole_worker.h
 * The worker class for the Turbomole QM software
 */
#ifndef TURBOMOLE_WORKER_H
#define	TURBOMOLE_WORKER_H

namespace interaction {
  class QM_Worker;
  /**
   * @class Turbomole_Worker
   * a worker class which calls the Turbomole software
   */
  class Turbomole_Worker : public QM_Worker {
  public:
    /**
     * Constructor
     */
    Turbomole_Worker() : QM_Worker("Turbomole Worker") {}
    /**
     * Destructor
     */
    virtual ~Turbomole_Worker();
    /**
     * initialise the QM worker
     * @return 0 if successful, non-zero on failure
     */
    virtual int init(topology::Topology & topo,
            configuration::Configuration & conf,
            simulation::Simulation & sim);
    /**
     * run a QM job in Turbomole
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
  };
}

#endif	/* TURBOMOLE_WORKER_H */

