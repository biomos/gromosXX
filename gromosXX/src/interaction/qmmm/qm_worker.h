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
    virtual int init(simulation::Simulation & sim) = 0;

    /**
     * run the QM worker
     */
    virtual int run_QM(topology::Topology& topo
                     , configuration::Configuration& conf
                     , simulation::Simulation& sim, interaction::QM_Zone & qm_zone);

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
     * Write input file for QM
     */
    virtual int write_input(const topology::Topology& topo
                          , const configuration::Configuration& conf
                          , const simulation::Simulation& sim
                          , const interaction::QM_Zone & qm_zone);

    /**
     * Open input file for QM
     */
    virtual int open_input(std::ofstream & inputfile_stream, const std::string & input_file);
    
    /**
     * Call external QM program
     */
    virtual int system_call();

    /**
     * read QM output files
     */
    virtual int read_output(topology::Topology& topo
                          , configuration::Configuration& conf
                          , simulation::Simulation& sim
                          , interaction::QM_Zone & qm_zone);
    
    /**
     * Open QM output file
     */
    virtual int open_output(std::ifstream & outputfile_stream, const std::string & output_file);
    
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

