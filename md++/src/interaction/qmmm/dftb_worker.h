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
 * @file dftb_worker.h
 * The worker class for the DFTB+ QM software
 */
#ifndef INCLUDED_DFTB_WORKER_H
#define	INCLUDED_DFTB_WORKER_H

namespace simulation {
  class Parameter;
  class Simulation;
}

namespace interaction {
  class QM_Worker;
  class QM_Zone;
  /**
   * @class DFTB_Worker
   * a worker class which calls the DFTB software
   */
  class DFTB_Worker : public QM_Worker {
  public:
    /**
     * Constructor
     */
    DFTB_Worker();
    /**
     * Destructor
     */
    virtual ~DFTB_Worker() = default;
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
    simulation::Parameter::qmmm_struct::dftb_param_struct* param;
    
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
     * Write list of atom types and create a map
     * @param input_file_stream Stream pointing to coordinate input file
     * @param type_indices map of atomic numbers to atom type indices
     * @param qm_zone QM Zone
     */
    void write_atom_types_list(std::ofstream& input_file_stream
                             , std::map<unsigned, unsigned>& type_indices
                             , const interaction::QM_Zone& qm_zone) const;
    
    /**
     * Write QM atom line
     * @param inputfile_stream ofstream to input file
     * @param atomic_number atomic number of the atom
     * @param pos position of the atom
     */
    void write_qm_atom(std::ofstream& inputfile_stream
                     , const int id
                     , const int atomic_type_id
                     , const math::Vec& pos) const;

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
     * Call external QM program - DFTB
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
     * Parse QM forces
     * @param ofs ifstream from the output file
     * @param qm_zone QM Zone
     */
    int parse_qm_forces(std::ifstream& ofs, interaction::QM_Zone& qm_zone) const;

    /**
     * Parse MM forces
     * @param ofs ifstream from the output file
     * @param qm_zone QM Zone
     */
    int parse_mm_forces(std::ifstream& ofs, interaction::QM_Zone& qm_zone) const;

    /**
     * Parse force line of QM or MM atom
     * @param ofs ifstream from the output file
     * @param force reference for writing the force on atom
     */
    int parse_force(std::ifstream& ofs, math::Vec& force) const;
  };
}
#endif	/* DFTB_WORKER_H */

