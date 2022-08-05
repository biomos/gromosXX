/**
 * @file gaussian_worker.h
 * The worker class for the Gaussian QM software
 */
#ifndef INCLUDED_GAUSSIAN_WORKER_H
#define	INCLUDED_GAUSSIAN_WORKER_H

namespace simulation {
  class Parameter;
  class Simulation;
}

namespace interaction {
  class QM_Worker;
  class QM_Zone;
  /**
   * @class Gaussian_Worker
   * a worker class which calls the Gaussian software
   */
  class Gaussian_Worker : public QM_Worker {
  public:
    /**
     * Constructor
     */
    Gaussian_Worker();
    /**
     * Destructor
     */
    virtual ~Gaussian_Worker() = default;
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
                   , const interaction::QM_Zone& qm_zone); 

  private:
    /**
     * Pointer to simulation parameters
     */
    simulation::Parameter::qmmm_struct::gaussian_param_struct* param;

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
     * System call
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
                  , const math::Vec& pos);

    /**
     * Write MM atom
     */
    void write_mm_atom(std::ofstream& inputfile_stream
                      , const math::Vec& pos
                      , const double charge);

    /**
     * Write MM position
     */
    void write_mm_pos(std::ofstream& inputfile_stream
                      , const math::Vec& pos);

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
     * Parse forces
     */
    int parse_forces(const simulation::Simulation& sim
                  , std::ifstream& ofs
                  , interaction::QM_Zone& qm_zone) const;
    
    /**
     * Parse force line
     */
    int parse_force(std::ifstream& ofs
                  , math::Vec& force) const;
  };
}

#endif	/* GAUSSIAN_WORKER_H */

