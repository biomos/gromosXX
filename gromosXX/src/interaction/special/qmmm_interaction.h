/**
 * @file qmmm_interaction.h
 * QM/MM interaction
 */

#ifndef QMMM_INTERACTION_H
#define	QMMM_INTERACTION_H

namespace interaction {
  class QM_Worker;
  class QM_Storage;

  /**
   * @class QMMM_Interaction
   * calculates the QM/MM interaction
   */ class QMMM_Interaction : public Interaction {
  public:

    /**
     * Constructor.
     */
    QMMM_Interaction() : Interaction("QM/MM Interaction"), worker(NULL) {}
    /**
     * Destructor.
     */
    virtual ~QMMM_Interaction();

    /**
     * init
     */
    virtual int init(topology::Topology &topo,
            configuration::Configuration &conf,
            simulation::Simulation &sim,
            std::ostream &os = std::cout,
            bool quiet = false);
    /**
     * calculate the interactions.
     */
    virtual int calculate_interactions(topology::Topology & topo,
            configuration::Configuration & conf,
            simulation::Simulation & sim);
  protected:
    QM_Worker * worker;
    QM_Storage storage;
  };
} // interaction

#endif	/* QMMM_INTERACTION_H */

