/**
 * @file in_parameter.h
 * read in a G96 input file.
 */

#ifndef INCLUDED_IN_PARAMETER_H
#define INCLUDED_IN_PARAMETER_H

#include "instream.h"

class Parameter;

/**
 * @class In_Parameter
 * reads in an input file and parses it.
 */
class In_Parameter : public GInStream {

public:
  /**
   * Default Constructor.
   */
  In_Parameter() {};
  /**
   * Constructor.
   * read in the complete file at once.
   */
  In_Parameter(std::istream& is) : GInStream(is) { readStream(); };

  /**
   * Store standard parameters in the simulation.
   */
  void read(Parameter &param, std::ostream & os = std::cout);
    
  /**
   * read minimise block.
   */
  void read_MINIMISE(Parameter &param, std::ostream & os = std::cout);

  /**
   * read step block.
   */
  void read_STEP(Parameter &param, std::ostream & os = std::cout);

  /**
   * read the SYSTEM block.
   */
  void read_SYSTEM(Parameter &param, std::ostream & os = std::cout);
    
  /**
   * read the SHAKE block.
   */
  void read_SHAKE(Parameter &param, std::ostream & os = std::cout);

  /**
   * read the CONSTRAINTS block.
   */
  void read_CONSTRAINTS(Parameter &param, std::ostream & os = std::cout);

  /**
   * read the FLEXCON block.
   */
  // void read_FLEXCON(Parameter &param, std::ostream & os = std::cout);

  /**
   * read FORCE block.
   */
  void read_FORCE(Parameter &param, std::ostream & os = std::cout);

  /**
   * read ENERGYGROUP block. (right now from the FORCE block)
   */
  void read_ENERGYGROUP(Parameter &param, std::ostream & os = std::cout);

  /**
   * read FORCEFIELD block.
   */
  void read_FORCEFIELD(Parameter &param, std::ostream & os = std::cout);

  /**
   * read PRINT block.
   */
  void read_PRINT(Parameter &param, std::ostream & os = std::cout);

  /**
   * read WRITE block.
   */
  void read_WRITE(Parameter &param, std::ostream & os = std::cout);

  /**
   * read the PCOUPLE03 block.
   */
  void read_PCOUPLE(Parameter &param, std::ostream & os = std::cout);
    
  /**
   * read the BOUNDARY block.
   */
  void read_BOUNDARY(Parameter &param, std::ostream & os = std::cout);

  /**
   * read the PERTURB block.
   */
  void read_PERTURB(Parameter &param, std::ostream & os = std::cout);

  /**
   * read START block.
   */
  void read_START(Parameter &param, std::ostream & os = std::cout);

  /**
   * read CENTREOFMASS block.
   */
  void read_CENTREOFMASS(Parameter &param, std::ostream & os = std::cout);

  /**
   * read PLIST block.
   */
  void read_PLIST(Parameter &param, std::ostream & os = std::cout);

  /**
   * read POSRES block.
   */
  void read_POSREST(Parameter &param, std::ostream & os = std::cout);

  /**
   * read DISTREST block.
   */
  void read_DISTREST(Parameter &param, std::ostream & os = std::cout);

  /**
   * read LONGRANGE block.
   */
  void read_LONGRANGE(Parameter &param, std::ostream & os = std::cout);

  /**
   * read SUBMOLECULES block.
   */
  void read_SUBMOLECULES(Parameter &param, std::ostream & os = std::cout);

  /**
   * read MULTIBATH block.
   */
  void read_MULTIBATH(Parameter &param, std::ostream & os = std::cout);

  /**
   * read J-VAL block.
   */
  void read_JVALUE(Parameter &param, std::ostream & os = std::cout);

  /**
   * read PSCALE block.
   */
  void read_PSCALE(Parameter &param, std::ostream & os = std::cout);

  /**
   * read ROTTRANS block.
   */
  void read_ROTTRANS(Parameter &param, std::ostream & os = std::cout);

  /**
   * read REPLICA block.
   */
  void read_REPLICA(Parameter &param, std::ostream & os = std::cout);

  /**
   * read REPLICA03 block.
   */
  void read_REPLICA03(Parameter &param, std::ostream & os = std::cout);

  /**
   * read INNERLOOP block.
   */
  void read_INNERLOOP(Parameter &param, std::ostream & os = std::cout);

  /**
   * read MULTICELL block.
   */
  void read_MULTICELL(Parameter & param, std::ostream & os = std::cout);

  /**
   * read ANALYZE block.
   */
  void read_ANALYZE(Parameter & param, std::ostream & os = std::cout);
    
};
  
#endif
