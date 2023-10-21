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
 * @file orca_worker.h
 * The worker class for the Orca QM software
 */
#ifndef INCLUDED_ORCA_WORKER_H
#define	INCLUDED_ORCA_WORKER_H

namespace interaction {
  class QM_Worker;
  class QM_Zone;
  /**
   * @class Orca_Worker
   * a worker class which calls the Orca software
   */
  class Orca_Worker : public QM_Worker {
  public:
    /**
     * Constructor
     */
    Orca_Worker();

    /**
     * Destructor
     */
    virtual ~Orca_Worker() = default;

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
    simulation::Parameter::qmmm_struct::orca_param_struct* param;
    
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
    int process_input(const topology::Topology& topo
                  , const configuration::Configuration& conf
                  , const simulation::Simulation& sim
                  , const interaction::QM_Zone& qm_zone) override;

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
     * @param qm_zone The QM zone
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
     * @param qm_zone The QM zone
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
     * Call external QM program - Orca
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
     * Parse gradient line of QM or MM atom (x,y,z components in one line)
     * @param ofs ifstream from the output file
     * @param force reference for writing the force
     */
    int parse_gradient1D(std::ifstream& ofs, math::Vec& force) const;

    /**
     * Parse gradient line of QM or MM atom (x,y,z components in three lines)
     * @param ofs ifstream from the output file
     * @param force reference for writing the force
     */
    int parse_gradient3D(std::ifstream& ofs, math::Vec& force) const;

    /**
     * Helper function to write the header in coordinate file trajectory
     */
    void write_coordinate_header(std::ofstream& ifs
                               , const QM_Zone& qm_zone) const override;
  };
}

#endif	/* ORCA_WORKER_H */

