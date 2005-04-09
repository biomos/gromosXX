/**
 * @file in_parameter.h
 * read in a G96 input file.
 */

#ifndef INCLUDED_IN_PARAMETER_H
#define INCLUDED_IN_PARAMETER_H

#include <gromosXX/io/instream.h>

namespace io {

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
    void read(simulation::Parameter &param, std::ostream & os = std::cout);
    
    /**
     * read minimise block.
     */
    void read_MINIMISE(simulation::Parameter &param, std::ostream & os = std::cout);

    /**
     * read step block.
     */
    void read_STEP(simulation::Parameter &param, std::ostream & os = std::cout);

    /**
     * read the SYSTEM block.
     */
    void read_SYSTEM(simulation::Parameter &param, std::ostream & os = std::cout);
    
    /**
     * read the SHAKE block.
     */
    void read_SHAKE(simulation::Parameter &param, std::ostream & os = std::cout);

    /**
     * read the CONSTRAINTS block.
     */
    void read_CONSTRAINTS(simulation::Parameter &param, std::ostream & os = std::cout);

    /**
     * read the FLEXCON block.
     */
    // void read_FLEXCON(simulation::Parameter &param, std::ostream & os = std::cout);

    /**
     * read FORCE block.
     */
    void read_FORCE(simulation::Parameter &param, std::ostream & os = std::cout);

    /**
     * read ENERGYGROUP block. (right now from the FORCE block)
     */
    void read_ENERGYGROUP(simulation::Parameter &param, std::ostream & os = std::cout);

    /**
     * read FORCEFIELD block.
     */
    void read_FORCEFIELD(simulation::Parameter &param, std::ostream & os = std::cout);

    /**
     * read PRINT block.
     */
    void read_PRINT(simulation::Parameter &param, std::ostream & os = std::cout);

    /**
     * read WRITE block.
     */
    void read_WRITE(simulation::Parameter &param, std::ostream & os = std::cout);

    /**
     * read the PCOUPLE03 block.
     */
    void read_PCOUPLE(simulation::Parameter &param, std::ostream & os = std::cout);
    
    /**
     * read the BOUNDARY block.
     */
    void read_BOUNDARY(simulation::Parameter &param, std::ostream & os = std::cout);

    /**
     * read the PERTURB block.
     */
    void read_PERTURB(simulation::Parameter &param, std::ostream & os = std::cout);

    /**
     * read START block.
     */
    void read_START(simulation::Parameter &param, std::ostream & os = std::cout);

    /**
     * read CENTREOFMASS block.
     */
    void read_CENTREOFMASS(simulation::Parameter &param, std::ostream & os = std::cout);

    /**
     * read PLIST block.
     */
    void read_PLIST(simulation::Parameter &param, std::ostream & os = std::cout);

    /**
     * read POSRES block.
     */
    void read_POSREST(simulation::Parameter &param, std::ostream & os = std::cout);

    /**
     * read DISTREST block.
     */
    void read_DISTREST(simulation::Parameter &param, std::ostream & os = std::cout);

    /**
     * read LONGRANGE block.
     */
    void read_LONGRANGE(simulation::Parameter &param, std::ostream & os = std::cout);

    /**
     * read SUBMOLECULES block.
     */
    void read_SUBMOLECULES(simulation::Parameter &param, std::ostream & os = std::cout);

    /**
     * read MULTIBATH block.
     */
    void read_MULTIBATH(simulation::Parameter &param, std::ostream & os = std::cout);

    /**
     * read J-VAL block.
     */
    void read_JVALUE(simulation::Parameter &param, std::ostream & os = std::cout);

    /**
     * read PSCALE block.
     */
    void read_PSCALE(simulation::Parameter &param, std::ostream & os = std::cout);

    /**
     * read ROTTRANS block.
     */
    void read_ROTTRANS(simulation::Parameter &param, std::ostream & os = std::cout);

    /**
     * read REPLICA block.
     */
    void read_REPLICA(simulation::Parameter &param, std::ostream & os = std::cout);

    /**
     * read REPLICA03 block.
     */
    void read_REPLICA03(simulation::Parameter &param, std::ostream & os = std::cout);

    /**
     * read INNERLOOP block.
     */
    void read_INNERLOOP(simulation::Parameter &param, std::ostream & os = std::cout);

    /**
     * read MULTICELL block.
     */
    void read_MULTICELL(simulation::Parameter & param, std::ostream & os = std::cout);
    
  };
  
} // io

#endif
