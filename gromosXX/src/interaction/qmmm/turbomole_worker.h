/**
 * @file turbomole_worker.h
 * The worker class for the Turbomole QM software
 */
#ifndef INCLUDED_TURBOMOLE_WORKER_H
#define	INCLUDED_TURBOMOLE_WORKER_H

namespace simulation {
  class Parameter;
  class Simulation;
}

namespace interaction {
  class QM_Worker;
  class QM_Zone;
  /**
   * @class Turbomole_Worker
   * a worker class which calls the Turbomole software
   */
  class Turbomole_Worker : public QM_Worker {
  public:
    /**
     * Constructor
     */
    Turbomole_Worker();
    /**
     * Destructor
     */
    virtual ~Turbomole_Worker() = default;
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
     * Current working directory (where GROMOS is called from)
     */
    std::string cwd;

    /**
     * Pointer to simulation parameters
     */
    simulation::Parameter::qmmm_struct::turbomole_param_struct* param;

    /**
     * Write input file for the QM program
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
     * Write QM atom line
     * @param inputfile_stream ofstream to input file
     * @param atomic_number atomic number of the atom
     * @param pos position of the atom
     */
    void write_qm_atom(std::ofstream& inputfile_stream
                     , const int atomic_number
                     , const math::Vec& pos) const;

    /**
     * Write MM atom line
     * @param inputfile_stream ofstream to the input file
     * @param pos position of the atom
     * @param charge charge of the atom
     */
    void write_mm_atom(std::ofstream& inputfile_stream
                     , const math::Vec& pos
                     , const double charge) const;

    /**
     * Call external QM program - Turbomole
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
     * Parse charges
     * @param ofs ifstream from the output file
     * @param qm_zone QM Zone
     */
    int parse_charges(std::ifstream& ofs, interaction::QM_Zone& qm_zone) const;

    /**
     * Parse energy
     * @param ofs ifstream from the output file
     * @param qm_zone QM Zone
     */
    int parse_energy(std::ifstream& ofs, interaction::QM_Zone& qm_zone) const;

    /**
     * Parse QM gradients
     * @param ofs ifstream from the output file
     * @param qm_zone QM Zone
     */
    int parse_qm_gradients(std::ifstream& ofs, interaction::QM_Zone& qm_zone) const;

    /**
     * Parse MM gradients
     * @param ofs ifstream from the output file
     * @param qm_zone QM Zone
     */
    int parse_mm_gradients(std::ifstream& ofs, interaction::QM_Zone& qm_zone) const;

    /**
     * Parse gradient line of QM or MM atom
     * @param ofs ifstream from the output file
     * @param force reference for writing the force
     */
    int parse_gradient(std::ifstream& ofs, math::Vec& force) const;
  };
}
#endif	/* TURBOMOLE_WORKER_H */

