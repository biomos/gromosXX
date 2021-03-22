/**
 * @file qmmm_interaction.h
 * QM/MM interaction
 */

#ifndef QMMM_INTERACTION_H
#define	QMMM_INTERACTION_H

namespace interaction {
  class QM_Worker;
  class QM_Storage;
  class MM_Atom;

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
    
    /**
     * Add the electric field contributions of the QM part
     * to the MM atoms.
     * 
     * @param electric_field the electric field
     * @return zero on success, non-zero on failure.
     */
    virtual int add_electric_field_contribution(topology::Topology & topo,
            configuration::Configuration & conf,
            simulation::Simulation & sim,
            math::VArray & electric_field);
    
    /**
     * prepare for a force/energy or electric field calculation 
     * @return zero on success
     */
    virtual int prepare(topology::Topology& topo,
        configuration::Configuration& conf,
        simulation::Simulation& sim);
  protected:
    QM_Worker * worker;
    QM_Storage storage;
    math::VArray qm_pos;
    std::vector<interaction::MM_Atom> mm_atoms;
  };
} // interaction

#endif	/* QMMM_INTERACTION_H */

