/**
 * @file argument.h
 * argument parsing.
 */

#ifndef INCLUDED_ARGUMENT_H
#define INCLUDED_ARGUMENT_H

class Known;

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
  // not implemented
  Argument(const Argument &);
  Argument &operator=(const Argument &);
public:
  /**
   * Argument constructor.
   * Details.
   */
  Argument();
  /**
   * parse the command line arguments.
   */
  int parse(int argc, char **argv, Known &known, bool empty_ok = false); 
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

private:
  std::string d_usage;
  std::string d_prog;
  std::set<std::string> d_known;

  std::string const empty;

};

#endif

