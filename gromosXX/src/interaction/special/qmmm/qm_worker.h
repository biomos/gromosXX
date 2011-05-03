/**
 * @file qm_worker.h
 * the worker class which calls the QM program and collectes the data.
 */
#ifndef QM_WORKER_H
#define	QM_WORKER_H

namespace util {
  class Algorithm_Timer;
}

namespace interaction {
  class QM_Storage;
  class MM_Atom;
  /**
   * @class QM_Worker
   * interface for a class implementing a call to an external QM software
   */
  class QM_Worker {
  public:
    /**
     * Constructor
     */
    QM_Worker(std::string name) : m_name(name), m_timer(name) {}
    /**
     * Destructor
     */
    virtual ~QM_Worker() {}

    /**
     * initialise the QM worker
     * @return 0 if successful, non-zero on failure
     */
    virtual int init(topology::Topology & topo,
            configuration::Configuration & conf,
            simulation::Simulation & sim) {
      return 0;
    }

    /**
     * run the QM worker
     * @param[in] qm_pos a vector containing the QM atom positions
     * @param[in] mm_atoms a vector containing the MM atoms to include
     * @param[out] storage the storage for energies, forces etc...
     * @return 0 if successful, non-zero on failure
     */
    virtual int run_QM(topology::Topology & topo,
            configuration::Configuration & conf,
            simulation::Simulation & sim,
            const math::VArray & qm_pos,
            const std::vector<MM_Atom> & mm_atoms,
            interaction::QM_Storage & storage) = 0;

    /**
     * Get an instance of a QM worker for the provided parameters
     * @return the instance or NULL on failure
     */
    static QM_Worker * get_instance(const simulation::Simulation & sim);

    /**
     * accessor the the name of the QM worker
     */
    const std::string & name() const { return m_name; }
    /**
     * accessor to the timer
     */
    util::Algorithm_Timer & timer() { return m_timer; }
    /**
     * const accessor to the timer
     */
    const util::Algorithm_Timer & timer() const { return m_timer; }
  private:
    /**
     * name of the QM worker
     */
    std::string m_name;

    /**
     * timer of this worker
     */
    util::Algorithm_Timer m_timer;
  };
}

#endif	/* QM_WORKER_H */

