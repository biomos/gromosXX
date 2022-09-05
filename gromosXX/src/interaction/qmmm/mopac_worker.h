/**
 * @file mopac_worker.h
 * The worker class for the MOPAC QM software
 */
#ifndef INCLUDED_MOPAC_WORKER_H
#define	INCLUDED_MOPAC_WORKER_H

namespace interaction {
  class QM_Worker;
  class QM_Zone;
  /**
   * @class MOPAC_Worker
   * a worker class which calls the MOPAC software
   */
  class MOPAC_Worker : public QM_Worker {

  public:
    /**
     * Constructor
     */
    MOPAC_Worker();

    /**
     * Destructor
     */
    virtual ~MOPAC_Worker() = default;

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
    simulation::Parameter::qmmm_struct::mopac_param_struct* param;

    /**
     * Write input file for QM program
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
     * System call - MOPAC
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
                  , const int var_flag = 1) const;

    /**
     * Write potential from MM atoms
     */
    void write_mm_potential(std::ofstream& inputfile_stream
                          , const int atomic_number
                          , const math::Vec& pos
                          , double potential) const;
    
    /**
     * Calculate total potential from MM atoms on QM atom
     */
    double total_potential(const topology::Topology& topo
                         , const simulation::Simulation& sim
                         , const QM_Zone& qm_zone
                         , const QM_Atom& qm_atom) const;
                         
    /**
     * Get a set of excluded MM atoms for link atom
     */
    void get_excluded_mm(const topology::Topology& topo
                       , const simulation::Simulation& sim
                       , const QM_Zone& qm_zone
                       , const QM_Atom& qm_atom
                       , std::set<unsigned> excluded) const;
    /**
     * Calculate potential from MM atom at the given position
     */
    double pair_potential(const math::Vec& pos
                        , const MM_Atom& mm_atom) const;

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
     * Parse gradients of QM atoms
     */
    int parse_qm_gradients(const simulation::Simulation& sim
                         , std::ifstream& ofs
                         , interaction::QM_Zone& qm_zone) const;

    /**
     * calculate forces between QM and MM atoms
     */
    void calculate_mm_forces(const topology::Topology& topo
                           , const simulation::Simulation& sim
                           , interaction::QM_Zone& qm_zone) const;

    /**
     * calculate force of QM/MM pair
     */
    inline void calculate_pair_force(const math::Vec& qm_pos
                                   , const math::Vec& mm_pos
                                   , math::Vec& qm_force
                                   , math::Vec& mm_force
                                   , const double qmq_mmq_four_pi_eps_i) const;
  };
}

#endif	/* MOPAC_WORKER_H */

