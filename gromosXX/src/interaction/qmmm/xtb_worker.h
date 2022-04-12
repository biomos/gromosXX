/**
 * @file xtb_worker.h
 * The worker class for the XTB software
 */
#ifndef INCLUDED_XTB_WORKER_H
#define	INCLUDED_XTB_WORKER_H

namespace interaction {
  class QM_Worker;
  class QM_Zone;
  /**
   * @class XTB_Worker
   * a worker class which calls the XTB software
   */
  class XTB_Worker : public QM_Worker {
  public:
    /**
     * Constructor
     */
    XTB_Worker();

    /**
     * Destructor
     */
    virtual ~XTB_Worker() = default;

    /**
     * initialise the QM worker
     * @return 0 if successful, non-zero on failure
     */
    virtual int init(simulation::Simulation& sim);

  private:
    /**
     * Pointer to simulation parameters
     */
    simulation::Parameter::qmmm_struct::xtb_param_struct* param;

    /**
     * Tests if a temporary file should be generated and performs some sanity checks
     * @param file_name Reference to the param structure file
     * @param new_name Potential new name for the file (for temporary files)
     * @return int error code
     */
    int initialize_file(std::string& file_name, const std::string& new_name);

    /**
     * Write input files for QM program
     * @param topo Topology
     * @param conf Configuration
     * @param sim Simulation
     * @param qm_zone QM Zone
     */
    int write_input(const topology::Topology& topo
                  , const configuration::Configuration& conf
                  , const simulation::Simulation& sim
                  , const interaction::QM_Zone& qm_zone);

    /**
     * Write input parameter file
     * @param inputfile_stream ofstream to input file
     * @param topo Topology
     * @param conf Configuration
     * @param sim Simulation
     * @return int error code
     */
    int write_input_parameters(std::ofstream& inputfile_stream
                             , const topology::Topology& topo
                             , const configuration::Configuration& conf
                             , const simulation::Simulation& sim
                             , const interaction::QM_Zone& qm_zone);

    /**
     * Write input coordinates file (xyz file format)
     * @param inputfile_stream ofstream to input file
     * @param topo Topology
     * @param conf Configuration
     * @param sim Simulation
     * @return int error code
     */
    int write_input_coordinates(std::ofstream& inputfile_stream
                              , const topology::Topology& topo
                              , const configuration::Configuration& conf
                              , const simulation::Simulation& sim
                              , const interaction::QM_Zone& qm_zone);

    /**
     * Write input file for point charges
     * @param inputfile_stream ofstream to input file
     * @param topo Topology
     * @param conf Configuration
     * @param sim Simulation
     * @return int error code
     */
    int write_input_pointcharges(std::ofstream& inputfile_stream
                               , const topology::Topology& topo
                               , const configuration::Configuration& conf
                               , const simulation::Simulation& sim
                               , const interaction::QM_Zone& qm_zone);

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
     * Call external QM program - XTB
     */
    int system_call();

    /**
     * Read output file from the QM program
     * @param topo Topology
     * @param conf Configuration
     * @param sim Simulation
     * @param qm_zone QM Zone
     */
    int read_output(topology::Topology& topo
                  , configuration::Configuration& conf
                  , simulation::Simulation& sim
                  , interaction::QM_Zone& qm_zone);
    
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
  };
}

#endif	/* XTB_WORKER_H */

