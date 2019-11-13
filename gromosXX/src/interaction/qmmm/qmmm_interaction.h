/**
 * @file qmmm_interaction.h
 * QM/MM interaction
 */

#ifndef INCLUDED_QMMM_INTERACTION_H
#define	INCLUDED_QMMM_INTERACTION_H




namespace io {
  class IFP;
}

#include <interaction/nonbonded/interaction/nonbonded_parameter.h>

namespace interaction {
  class QM_Worker;
  class QM_Zone;
  class QMMM_Nonbonded_Set;

  /**
   * @class QMMM_Interaction
   * calculates the QM/MM interaction
   * This will be a hybrid of new Interaction and 
   * Nonbonded_Interaction
   */
  class QMMM_Interaction : public Interaction {
  public:

    /**
     * Constructor.
     */
    QMMM_Interaction();
    /**
     * Destructor.
     */
    virtual ~QMMM_Interaction();

    /**
     * init
     */
    virtual int init(topology::Topology &topo
                   , configuration::Configuration &conf
                   , simulation::Simulation &sim
                   , std::ostream &os = std::cout
                   , bool quiet = false) override;
    /**
     * calculate the interactions.
     */
    virtual int calculate_interactions(topology::Topology & topo,
            configuration::Configuration & conf,
            simulation::Simulation & sim) override;
    
    /**
     * Add the electric field contributions of the QM part
     * to the MM atoms.
     * 
     * @param electric_field the electric field
     * @return zero on success, non-zero on failure.
     */
    /*virtual int add_electric_field_contribution(topology::Topology & topo,
            configuration::Configuration & conf,
            simulation::Simulation & sim,
            math::VArray & electric_field);*/

    int update_electric_field(topology::Topology& topo,
                              configuration::Configuration& conf,
                              simulation::Simulation& sim) {return 0;};

    /**
     * QM Zone accessor
     */
    const QM_Zone & qm_zone() const {return *m_qm_zone;}

    /**
     * parameter accessor
     */
    const Nonbonded_Parameter & parameter() const {return m_parameter;}

    /**
     * parameter mutator
     */
    Nonbonded_Parameter & parameter() {return m_parameter;}

    /**
     * timer accessor
     */
    util::Algorithm_Timer & timer() {return m_timer;}

    /**
     * timing information.
     */
    virtual void print_timing(std::ostream & os);

  protected:
    int init_nonbonded(topology::Topology& topo
                      , configuration::Configuration& conf
                      , simulation::Simulation& sim
                      , std::ostream &os
                      , bool quiet);

    int calculate_nonbonded(topology::Topology& topo
                          , configuration::Configuration& conf
                          , simulation::Simulation& sim);
    
    void remove_bonded_terms(topology::Topology& topo
                           , std::ostream& os
                           , bool quiet);

    void remove_exclusions(topology::Topology& topo
                         , const simulation::Simulation& sim
                         , std::ostream &os
                         , bool quiet);

    /**
     * store data from sets into the configuration
     */
    void store_set_data(const topology::Topology& topo
                      , configuration::Configuration& conf
                      , const simulation::Simulation& sim);
    
    /**
     * print the pairlist
     */
    int print_pairlist(const topology::Topology& topo
                     , std::ostream& os = std::cout);

    /**
     * a vector of nonbonded QMMM sets
     */
    std::vector<QMMM_Nonbonded_Set *> m_qmmm_nonbonded_set;

    /**
     * nonbonded parameter
     */
    Nonbonded_Parameter m_parameter;

    /**
     * number of sets to create (OpenMP size)
     */
    unsigned m_set_size;

    /**
     * MPI rank
     */
    unsigned m_rank;

    /**
     * MPI size
     */
    unsigned m_size;

    /**
     * QM worker
     */
    QM_Worker * m_worker;

    /**
     * QM zone
     */
    QM_Zone * m_qm_zone;
  };
  
} // interaction

#endif	/* QMMM_INTERACTION_H */

