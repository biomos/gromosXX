/*
 * This file is part of GROMOS.
 * 
 * Copyright (c) 2011, 2012, 2016, 2018, 2021, 2023 Biomos b.v.
 * See <https://www.gromos.net> for details.
 * 
 * GROMOS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

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
     * Initialize the QM atoms and QM capping atoms
     * 
     * @param qm_zone The QM zone
     */
    void init_input_atom_types(const interaction::QM_Zone& qm_zone);
    
    /**
     * Extracts the coordinates of the QM zone, transfers them into a one-dimensional
     * array, and scales them according to the given conversion factor
     * 
     * @param topo Topology
     * @param conf Configuration
     * @param sim Simulation
     * @param qm_zone The QM zone
     */
    int process_input_coordinates(const topology::Topology& topo
                                 , const configuration::Configuration& conf
                                 , const simulation::Simulation& sim
                                 , const interaction::QM_Zone& qm_zone);

    /**
     * Extracts the coordinates of the MM zone, transfers them into a one-dimensional
     * array, and scales them according to the given conversion factor
     * 
     * @param topo Topology
     * @param conf Configuration
     * @param sim Simulation
     * @param qm_zone The QM zone
     */
    int process_input_pointcharges(const topology::Topology& topo
                                  , const configuration::Configuration& conf
                                  , const simulation::Simulation& sim
                                  , const interaction::QM_Zone& qm_zone);

    /**
     * Call external QM program - XTB
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
     * @param qm_zone QM Zone
     */
    int parse_charges(interaction::QM_Zone& qm_zone) const;

    /**
     * Parse energy
     * @param qm_zone QM Zone
     */
    int parse_energy(interaction::QM_Zone& qm_zone) const;

    /**
     * Parse QM gradients
     * @param qm_zone QM Zone
     */
    int parse_qm_gradients(interaction::QM_Zone& qm_zone) const;

    /**
     * Parse MM gradients
     * @param qm_zone QM Zone
     */
    int parse_mm_gradients(interaction::QM_Zone& qm_zone) const;

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
     * Write QM atom line
     * @param inputfile_stream ofstream to input file
     * @param atomic_number atomic number of the atom
     * @param pos position of the atom
     */
    void write_qm_atom(std::ofstream& inputfile_stream
                     , const int atomic_number
                     , const math::Vec& pos) const override;

    /**
     * Write MM atom line
     * @param inputfile_stream ofstream to the input file
     * @param pos position of the atom
     * @param charge charge of the atom
     */
    void write_mm_atom(std::ofstream& inputfile_stream
                     , const int atomic_number
                     , const math::Vec& pos
                     , const double charge) const override;
    
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

    /**
     * Number of point charges
     */
    int ncharges;

    /**
     * These correspond to element types - required for XTB
     * calculations to match internally hard-coded hardness
     * parameters: https://xtb-docs.readthedocs.io/en/latest/pcem.html and
     * https://github.com/grimme-lab/xtb/pull/584#discussion_r814312190 
     */
    std::vector<int> numbers;

    /**
     * Charges of the point charges
     */
    std::vector<double> charges;

    /**
     * Cartesian coordinates of the point charges
     */
    std::vector<double> point_charges;
  };
}

#endif	/* XTB_WORKER_H */