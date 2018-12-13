/**
 * @file argument.h
 * argument parsing.
 */

#ifndef INCLUDED_ARGUMENT_H
#define INCLUDED_ARGUMENT_H

namespace util
{
  class Known;
}

namespace io{

  /**
   * topology argument
   */
  extern std::string argname_topo;
  /**
   * perturbation topology argument
   */
  extern std::string argname_pttopo;
  /**
   * configuration argument
   */
  extern std::string argname_conf;
  /**
   * input argument
   */
  extern std::string argname_input;
  /**
   * trajectory argument
   */
  extern std::string argname_trj;
  /**
   * final configuration argument
   */
  extern std::string argname_fin;
  /**
   * energy trajectory argument
   */
  extern std::string argname_tre;
  /**
   * velocity trajectory argument
   */
  extern std::string argname_trv;
  /**
   * force trajectory argument
   */
  extern std::string argname_trf;  
  /**
   * special trajectory argument
   */
  extern std::string argname_trs;
  /**
   * replica exchange trajectory argument
   */
  extern std::string argname_re;
  /**
   * free energy trajectory argument
   */
  extern std::string argname_trg;
  /**
   * block averaged energy trajectory argument
   */
  extern std::string argname_bae;
  /**
   * block averaged free energy trajectory argument
   */
  extern std::string argname_bag;
  /**
   * gzip compression argument
   */
  extern std::string argname_gzip;
  /**
   * Class Argument
   * Purpose: Parse arguments from the command line or an input file.
   * 
   * Description:
   * This class is used to parse arguments from the command line or an input file .
   * 
   * @class Argument
   * @version $Date: Mon Jul 15 14:17:41 MEST 2002
   * @author  R. Buergi
   */

  class Argument: public std::multimap<std::string,std::string>{
    
  public:
    /**
     * Argument constructor.
     * Details.
     */
    Argument();
    
    // not implemented
    //Argument(const Argument &);
    //Argument &operator=(const Argument &);
    
    /**
     * parse the command line arguments.
     */
    int parse(int argc, char **argv, util::Known &known, bool empty_ok = false); 
    /**
     * Argument destructor.
     * Details.
     */
    ~Argument();

    /** 
     * Checks for whether the argument string had num_arg arguments
     * @param &str Takes a std::string as argument.
     * @param num_args Integer of the number of arguments.
     * @return check Integer to check for failure.
     */
    int check(const std::string &str, int num_args=0)const;

    /**
     * Returns the number of arguments that follow string
     * @param &str Takes a std::string as argument.
     * @return the number of arguments found for this argument. Returns -1 
     * if string was not found at all in the argument list.
     */
    int count(const std::string &str)const;
  
    // This has to be in to fix a bug in gcc Solaris 2.6 (?)
    //   const const_iterator begin()const
    //   {return this->multimap<string,string>::begin();}

    /** 
     * Member operator >>.
     * Details.
     */
    int parse_line(std::istream &is);
  
    /** 
     * Member operator [] used to access the members.
     * Details.
     */
    const std::string &operator[](const std::string &str)const;
    
    /**
     * add some elements
     * @param key
     * @param values
     */
    void put(const std::string & key, const std::vector<std::string> & values);

  private:
    std::string d_usage;
    std::string d_prog;
    std::set<std::string> d_known;

    std::string const empty;

  };

}

#endif

