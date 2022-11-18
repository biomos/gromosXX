/**
 * @file ghost_worker.h
 * The worker class for the Ghost software (does nothing)
 */
#ifndef INCLUDED_GHOST_WORKER_H
#define	INCLUDED_GHOST_WORKER_H

namespace interaction {
  class QM_Worker;
  class QM_Zone;
  /**
   * @class Ghost_Worker
   * The worker class for the Ghost software (does nothing)
   */
  class Ghost_Worker : public QM_Worker {
    
  public:
    /**
     * Constructor
     */
    Ghost_Worker();

    /**
     * Destructor
     */
    virtual ~Ghost_Worker();

    /**
     * initialise the QM worker
     * @param topo Topology
     * @param conf Configuration
     * @param sim Simulation
     * @param qm_zone QM Zone
     * @return 0 if successful, non-zero on failure
     */
    virtual int init(const topology::Topology& topo
                   , const configuration::Configuration& conf
                   , simulation::Simulation& sim
                   , const interaction::QM_Zone& qm_zone) override; 

  private:
    
    /**
     * Write input files for QM program
     * @param topo Topology
     * @param conf Configuration
     * @param sim Simulation
     * @param qm_zone QM Zone
     */
    int process_input(const topology::Topology& topo
                  , const configuration::Configuration& conf
                  , const simulation::Simulation& sim
                  , const interaction::QM_Zone& qm_zone) override;              

    /**
     * Do nothing
     */
    int run_calculation() override;
    
    /**
     * Read output file from the QM program
     * @param topo Topology
     * @param conf Configuration
     * @param sim Simulation
     * @param qm_zone QM Zone
     */
    int process_output(topology::Topology& topo
                  , configuration::Configuration& conf
                  , simulation::Simulation& sim
                  , interaction::QM_Zone& qm_zone) override;


    /**
     * Helper function to write the header in coordinate file trajectory
     */
    void write_coordinate_header(std::ofstream& ifs
                               , const QM_Zone& qm_zone) const override;

    /**
     * Helper function to write the footer in coordinate file trajectory
     */
    void write_coordinate_footer(std::ofstream& inputfile_stream) const override;
    
    /**
     * Pointer to simulation parameters
     */
    simulation::Parameter::qmmm_struct::ghost_param_struct* param;

  };
}

#endif	/* GHOST_WORKER_H */