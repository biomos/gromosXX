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
 * @file qm_worker.h
 * the worker class which calls the QM program and collectes the data.
 */
#ifndef INCLUDED_QM_WORKER_H
#define	INCLUDED_QM_WORKER_H

#define MAXPATH 4096

namespace util {
  class Algorithm_Timer;
}

namespace interaction {
  class QM_Zone;
  struct MM_Atom;
  /**
   * @class QM_Worker
   * interface for a class implementing a call to an external QM software
   */
  class QM_Worker {
  public:
    /**
     * Constructor
     */
    QM_Worker(std::string name);

    virtual ~QM_Worker();

    /**
     * Get an instance of a QM worker for the provided parameters
     * @return the instance or nullptr on failure
     */
    static QM_Worker * get_instance(const simulation::Simulation & sim);

    /**
     * initialise the QM worker
     */
    virtual int init(const topology::Topology& topo
                   , const configuration::Configuration& conf
                   , simulation::Simulation& sim
                   , const interaction::QM_Zone& qm_zone);

    /**
     * run the QM worker
     */
    virtual int run_QM(topology::Topology& topo
                     , configuration::Configuration& conf
                     , simulation::Simulation& sim
                     , interaction::QM_Zone & qm_zone);

    /**
     * accessor to the name of the QM worker
     */
    const std::string & name() const { return m_name; }

    /**
     * accessor to the timer
     */
    util::Algorithm_Timer & timer() { return m_timer; }

    /**
     * const accessor to the timer
     */
    const util::Algorithm_Timer & timer() const { return m_timer; }

    /**
     * Print units conversion factors
     */
    virtual void print_unit_factors(std::ostream & os) const {
      os << this->param->unit_factor_length << ", "
         << this->param->unit_factor_energy << ", "
         << this->param->unit_factor_force << ", "
         << this->param->unit_factor_charge;
      };
  
  protected:
    /**
     * the timer
     */
    util::Algorithm_Timer m_timer;
    
    /**
     * List of used temporary files to cleanup gracefully
     */
    std::set<std::string> tmp_files;

    /**
     * name of the QM worker
     */
    std::string m_name;

    /**
     * the pointer to QM parameters
     */
    simulation::Parameter::qmmm_struct::qm_param_struct* param;

    /**
     * Flag to enable transfer of coordinates from external QM program to GROMOS
     */
    bool minimisation;

    /**
     * Handle to the input coordinate trajectory (QM) 
     */
    mutable std::ofstream input_coordinate_stream;

    /**
     * Handle to the input point charge trajectory (QM) 
     */
    mutable std::ofstream input_point_charge_stream;

    /**
     * Handle to the output gradient trajectory (QM) 
     */
    mutable std::ofstream output_gradient_stream;

    /**
     * Handle to the output point charge gradient trajectory (QM) 
     */
    mutable std::ofstream output_point_charge_gradient_stream;

    /**
     * Handle to the output charges trajectory (QM) 
     */
    mutable std::ofstream output_charges_stream;
    
    /**
     * Write input file for QM
     */
    virtual int process_input(const topology::Topology& topo
                            , const configuration::Configuration& conf
                            , const simulation::Simulation& sim
                            , const interaction::QM_Zone & qm_zone);  
    
    /**
     * Open input file for QM
     */
    virtual int open_input(std::ofstream & inputfile_stream, const std::string & input_file) const;
    
    /**
     * Call external QM program
     */
    virtual int run_calculation();

    /**
     * Helper function to write information on step size in QM trajectory
     */
    virtual void write_step_size(std::ofstream& inputfile_stream
                               , const unsigned int step) const;

    /**
     * Helper function to write the header in coordinate file trajectory
     */
    virtual void write_coordinate_header(std::ofstream& inputfile_stream
                                       , const QM_Zone& qm_zone) const;

    /**
     * Helper function to write the footer in coordinate file trajectory
     */
    virtual void write_coordinate_footer(std::ofstream& inputfile_stream) const;

    /**
     * Writes a gradient (default implementation)
     */
    virtual void write_gradient(const math::Vec& gradient, 
                                std::ofstream& inputfile_stream) const;

    /**
     * Writes a QM atom (default implementation)
     */
    virtual void write_qm_atom(std::ofstream& inputfile_stream
                             , const int atomic_number
                             , const math::Vec& pos) const;

    /**
     * Writes a MM atom (default implementation)
     */
    virtual void write_mm_atom(std::ofstream& inputfile_stream
                             , const int atomic_number
                             , const math::Vec& pos
                             , const double charge) const;

    /**
     * Writes a charge
     */
    virtual void write_charge(std::ofstream& inputfile_stream
                            , const int atomic_number
                            , const double charge) const;
    
    /**
     * Writes the QM trajectory to special file (coordinates, point charges)
     */
    virtual void save_input(const unsigned int step
                          , const simulation::Simulation& sim
                          , const interaction::QM_Zone & qm_zone) const;

    /**
     * Writes the QM trajectory to special file (gradients, point charge gradients)
     */
    virtual void save_output(const unsigned int step
                           , const simulation::Simulation& sim
                           , const interaction::QM_Zone & qm_zone) const;

    /**
     * Writes QM input coordinates in QM units
     */
    virtual void save_input_coord(std::ofstream& inputfile_stream
                                , const unsigned int step
                                , const interaction::QM_Zone & qm_zone) const;

    /**
     * Writes QM input point charge coordinates in QM units
     */
    virtual void save_input_point_charges(std::ofstream& inputfile_stream
                                        , const unsigned int step
                                        , const unsigned int ncharges
                                        , const interaction::QM_Zone & qm_zone) const;

    /**
     * Writes QM gradients in QM units

     */
    virtual void save_output_gradients(std::ofstream& inputfile_stream
                                     , const unsigned int step
                                     , const interaction::QM_Zone & qm_zone) const;

    /**
     * Writes QM point charge gradients in QM units

     */
    virtual void save_output_pc_gradients(std::ofstream& inputfile_stream
                                        , const unsigned int step
                                        , const interaction::QM_Zone & qm_zone) const;

    /**
     * Writes QM charges in QM units

     */
    virtual void save_output_charges(std::ofstream& inputfile_stream
                                   , const unsigned int step
                                   , const interaction::QM_Zone & qm_zone) const;

    /**
     * read QM output files
     */
    virtual int process_output(topology::Topology& topo
                          , configuration::Configuration& conf
                          , simulation::Simulation& sim
                          , interaction::QM_Zone & qm_zone);
    
    /**
     * Open QM output file
     */
    virtual int open_output(std::ifstream & outputfile_stream, const std::string & output_file) const;
    
    /**
     * Replace exponent character to parse correctly
     */
    inline void defortranize(std::string& str) const {
      std::replace(str.begin(), str.end(), 'D', 'E');
    }

    /**
     * Calculate number of charges
     */
    virtual int get_num_charges(const simulation::Simulation& sim
                              , const interaction::QM_Zone & qm_zone) const;

    /**
     * Get current working directory
     */
    inline std::string getcwd() {
  #ifdef HAVE_GETCWD
      char buff[MAXPATH];
      if (::getcwd(buff, MAXPATH) == NULL) {
        io::messages.add("Cannot get current working directory. "
            "Path of the current working directory is too long.", 
            this->name(), io::message::error);
        return "";
      }
      std::string cwd(buff);
      return cwd;
  #else
      io::messages.add("getcwd function is not available on this platform.", 
          this->name(), io::message::error);
      return "";
  #endif
    }

    /**
     * Change directory
     * @param path directory to go to
     */
    inline int chdir(std::string path){
  #ifdef HAVE_CHDIR
      if (::chdir(path.c_str()) != 0) {
        io::messages.add("Cannot change into QM program working directory",
                this->name(), io::message::error);
        return -1;
      }
      return 0;
  #else
      io::messages.add("chdir function is not available on this platform.", 
            this->name(), io::message::error);
      return -1;
  #endif
    }
  };
}

#endif	/* QM_WORKER_H */

