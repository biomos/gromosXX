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
    void read(simulation::Parameter &param);
    
    /**
     * read minimise block.
     */
    void read_MINIMISE(simulation::Parameter &param);

    /**
     * read step block.
     */
    void read_STEP(simulation::Parameter &param);

    /**
     * read the SYSTEM block.
     */
    void read_SYSTEM(simulation::Parameter &param);
    
    /**
     * read the SHAKE block.
     */
    void read_SHAKE(simulation::Parameter &param);

    /**
     * read the CONSTRAINTS block.
     */
    void read_CONSTRAINTS(simulation::Parameter &param);

    /**
     * read the FLEXCON block.
     */
    // void read_FLEXCON(simulation::Parameter &param);

    /**
     * read FORCE block.
     */
    void read_FORCE(simulation::Parameter &param);

    /**
     * read ENERGYGROUP block. (right now from the FORCE block)
     */
    void read_ENERGYGROUP(simulation::Parameter &param);

    /**
     * read FORCEFIELD block.
     */
    void read_FORCEFIELD(simulation::Parameter &param);

    /**
     * read PRINT block.
     */
    void read_PRINT(simulation::Parameter &param);

    /**
     * read WRITE block.
     */
    void read_WRITE(simulation::Parameter &param);

    /**
     * read the PCOUPLE03 block.
     */
    void read_PCOUPLE(simulation::Parameter &param);
    
    /**
     * read the BOUNDARY block.
     */
    void read_BOUNDARY(simulation::Parameter &param);

    /**
     * read the PERTURB block.
     */
    void read_PERTURB(simulation::Parameter &param);

    /**
     * read START block.
     */
    void read_START(simulation::Parameter &param);

    /**
     * read CENTREOFMASS block.
     */
    void read_CENTREOFMASS(simulation::Parameter &param);

    /**
     * read PLIST block.
     */
    void read_PLIST(simulation::Parameter &param);

    /**
     * read POSRES block.
     */
    void read_POSREST(simulation::Parameter &param);

    /**
     * read LONGRANGE block.
     */
    void read_LONGRANGE(simulation::Parameter &param);

    /**
     * read SUBMOLECULES block.
     */
    void read_SUBMOLECULES(simulation::Parameter &param);

    /**
     * read MULTIBATH block.
     */
    void read_MULTIBATH(simulation::Parameter &param);

    /**
     * read J-VAL block.
     */
    void read_JVALUE(simulation::Parameter &param);

    /**
     * read PSCALE block.
     */
    void read_PSCALE(simulation::Parameter &param);

    /**
     * read ROTTRANS block.
     */
    void read_ROTTRANS(simulation::Parameter &param);
    
  };
  
} // io

#endif
