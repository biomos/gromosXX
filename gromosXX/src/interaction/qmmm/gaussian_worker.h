/**
 * @file gaussian_worker.h
 * The worker class for the Gaussian QM software
 */
#ifndef INCLUDED_GAUSSIAN_WORKER_H
#define	INCLUDED_GAUSSIAN_WORKER_H

//#include "qm_worker.h"

//#include "../../../simulation/simulation.h"

namespace simulation {
  struct gaussian_param_struct;
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
    virtual ~Gaussian_Worker();
    /**
     * initialise the QM worker
     * @return 0 if successful, non-zero on failure
     */
    virtual int init(simulation::Simulation& sim);
    /**
     * run a QM job in Gaussian
     * @param qm_pos a vector containing the QM atom positions
     * @param mm_atoms the MM atoms to include
     * @param storage the energies, forces, charges obtained
     * @return 0 if successful, non-zero if not.
     */

  private:
    /**
     * Pointer to simulation parameters
     */
    simulation::Parameter::qmmm_struct::gaussian_param_struct* param;

    /**
     * Write input file for QM
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
     * Parse gradients wrapper
     */
    int parse_gradients(const simulation::Simulation& sim
                      , std::ifstream& ofs
                      , interaction::QM_Zone& qm_zone);

    /**
     * Parse gradients
     */
    template<class AtomType>
    int _parse_gradients(std::ifstream& ofs, std::set<AtomType>& atom_set);

    /**
     * Parse gradient line
     */
    int parse_gradient(std::ifstream& ofs
                     , const int index
                     , math::Vec& force
                     , const double unit_factor);
  };

  /**
   * Parse gradients of polarisable MM atoms
   */
  template<>
  int Gaussian_Worker::_parse_gradients<interaction::MM_Atom>
        (std::ifstream& ofs, std::set<interaction::MM_Atom>& atom_set);
}

#endif	/* GAUSSIAN_WORKER_H */

