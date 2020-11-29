/**
 * @file mopac_worker.h
 * The worker class for the MOPAC QM software
 */
#ifndef INCLUDED_MOPAC_WORKER_H
#define	INCLUDED_MOPAC_WORKER_H

//#include "../../../simulation/simulation.h"

//#include <interaction/qmmm/qm_worker.h>
//#include <interaction/qmmm/qm_zone.h>

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
     * @return 0 if successful, non-zero on failure
     */
    virtual int init(simulation::Simulation& sim);

  private:
    /**
     * Pointer to simulation parameters
     */
    simulation::Parameter::qmmm_struct::mopac_param_struct* param;

    /**
     * Write input file for QM
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
     * System call
     */
    int system_call();

    /**
     * Read outputs
     */
    int read_output(topology::Topology& topo
                  , configuration::Configuration& conf
                  , simulation::Simulation& sim
                  , interaction::QM_Zone& qm_zone);

    /**
     * Write QM atom
     */
    void write_qm_atom(std::ofstream& inputfile_stream
                  , const int atomic_number
                  , const math::Vec& pos
                  , const int var_flag = 1);

    /**
     * Write potential from MM atoms
     */
    void write_mm_potential(std::ofstream& inputfile_stream
                          , const int atomic_number
                          , const math::Vec& pos
                          , double potential) const;
    
    /**
     * Calculate total potential from MM atoms at QM atom
     */
    double total_potential(const QM_Atom& qm_atom
                         , const std::set<MM_Atom>& mm_atoms) const;
    
    /**
     * Calculate total potential from MM atoms at link atom with possible exclusions
     */
    double link_total_potential(const topology::Topology& topo
                              , const simulation::Simulation& sim
                              , const QM_Link& link
                              , const std::set<MM_Atom>& mm_atoms) const;
    
    /**
     * Calculate potential from MM atom at the given position
     */
    double pair_potential(const math::Vec& pos
                        , const MM_Atom& mm_atom) const;

    /**
     * Parse charges
     */
    int parse_charges(std::ifstream& ofs, interaction::QM_Zone& qm_zone);

    /**
     * Parse coordinates
     */
    int parse_coordinates(std::ifstream& ofs, interaction::QM_Zone& qm_zone);

    /**
     * Parse energy
     */
    int parse_energy(std::ifstream& ofs, interaction::QM_Zone& qm_zone);

    /**
     * Parse gradients
     */
    int parse_gradients(const simulation::Simulation& sim
                      , std::ifstream& ofs
                      , interaction::QM_Zone& qm_zone);
  };
}

#endif	/* MOPAC_WORKER_H */

