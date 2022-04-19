/**
 * @file xtb_worker.h
 * The worker class for the XTB software
 */
#ifndef INCLUDED_XTB_WORKER_H
#define	INCLUDED_XTB_WORKER_H

#include "xtb.h"

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
    virtual ~XTB_Worker();

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
     * Write input files for QM program
     * @param topo Topology
     * @param conf Configuration
     * @param sim Simulation
     * @param qm_zone QM Zone
     */
    int process_input(const topology::Topology& topo
                  , const configuration::Configuration& conf
                  , const simulation::Simulation& sim
                  , const interaction::QM_Zone& qm_zone);              

    /**
     * Call external QM program - XTB
     */
    int run_calculation();

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

    /**
     * Pointer to simulation parameters
     */
    simulation::Parameter::qmmm_struct::xtb_param_struct* param;

    /**
     * Environment object to the XTB C API
     */
    xtb_TEnvironment env;

    /**
     * Calculator object to the XTB C API
     */
    xtb_TCalculator calc;

    /**
     * Molecule object to the XTB C API
     */
    xtb_TMolecule mol;
  
    /**
     * Results object to the XTB C API
     */
    xtb_TResults res;

    /**
     * Number of unpaired electrons
     */
    int uhf;

    /**
     * Charge of the QM zone
     */
    double charge;

    /**
     * Number of atoms in the QM zone 
     */
    int natoms;

    /**
     * Coordinates of the QM zone as one-dimensional vector, necessary to access C style array for xtb
     */
    std::vector<double> coord;

    /**
     * Atom types of the QM zone
     */
    std::vector<int> attyp;
  };
}

#endif	/* XTB_WORKER_H */

