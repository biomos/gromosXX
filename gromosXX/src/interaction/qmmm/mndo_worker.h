/**
 * @file mndo_worker.h
 * The worker class for the MNDO QM software
 */
#ifndef INCLUDED_MNDO_WORKER_H
#define	INCLUDED_MNDO_WORKER_H

namespace simulation {
  class Parameter;
  class Simulation;
}

namespace interaction {
  class QM_Worker;
  class QM_Zone;
  /**
   * @class MNDO_Worker
   * a worker class which calls the MNDO software
   */
  class MNDO_Worker : public QM_Worker {
  public:
    /**
     * Constructor
     */
    MNDO_Worker();
    /**
     * Destructor
     */
    virtual ~MNDO_Worker() = default;
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
     * Pointer to simulation parameters
     */
    simulation::Parameter::qmmm_struct::mndo_param_struct* param;

    /**
     * Write input file for QM
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
     * System call - MNDO
     */
    int run_calculation() override;

    /**
     * Read outputs
     */
    int process_output(topology::Topology& topo
                  , configuration::Configuration& conf
                  , simulation::Simulation& sim
                  , interaction::QM_Zone& qm_zone) override;

    /**
     * Write QM atom
     */
    void write_qm_atom(std::ofstream& inputfile_stream
                  , const int atomic_number
                  , const math::Vec& pos
                  , const int opt_flag = 0) const;

    /**
     * Write MM atom
     */
    void write_mm_atom(std::ofstream& inputfile_stream
                      , const math::Vec& pos
                      , const double charge) const;

    /**
     * Parse charges
     */
    int parse_charges(std::ifstream& ofs, interaction::QM_Zone& qm_zone) const;

    /**
     * Parse coordinates
     */
    int parse_coordinates(std::ifstream& ofs, interaction::QM_Zone& qm_zone) const;

    /**
     * Parse energy
     */
    int parse_energy(std::ifstream& ofs, interaction::QM_Zone& qm_zone) const;

    /**
     * Parse gradients
     */
    int parse_gradients(const simulation::Simulation& sim
                      , std::ifstream& ofs
                      , interaction::QM_Zone& qm_zone) const;

    /**
     * Parse gradient line
     */
    int parse_gradient(std::ifstream& ofs
                     , math::Vec& force) const;
  };
}

#endif	/* MNDO_WORKER_H */

