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
     * read the CONSTRAINT block.
     */
    void read_CONSTRAINT(simulation::Parameter &param, std::ostream & os = std::cout);

    /**
     * read the FLEXCON block.
     */
    // void read_FLEXCON(simulation::Parameter &param, std::ostream & os = std::cout);

    /**
     * read FORCE block.
     */
    void read_FORCE(simulation::Parameter &param, std::ostream & os = std::cout);
    
    /**
     * read COVALENTFORM block.
     */
    void read_COVALENTFORM(simulation::Parameter &param, std::ostream & os = std::cout);

    /**
     * read ENERGYGROUP block. (right now from the FORCE block)
     */
    void read_ENERGYGROUP(simulation::Parameter &param, std::ostream & os = std::cout);

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
     * read INITIALISE block.
     */
    void read_INITIALISE(simulation::Parameter &param, std::ostream & os = std::cout);

    /**
     * read CENTREOFMASS block.
     */
    void read_CENTREOFMASS(simulation::Parameter &param, std::ostream & os = std::cout);

    /**
     * read PAIRLIST block.
     */
    void read_PAIRLIST(simulation::Parameter &param, std::ostream & os = std::cout);

    /**
     * read POSRES block.
     */
    void read_POSREST(simulation::Parameter &param, std::ostream & os = std::cout);

    /**
     * read DISTREST block.
     */
    void read_DISTREST(simulation::Parameter &param, std::ostream & os = std::cout);

    /**
     * read DIHEDRALRES block.
     */
    void read_DIHREST(simulation::Parameter &param, std::ostream & os = std::cout);

    /**
     * read LONGRANGE block.
     */
    void read_LONGRANGE(simulation::Parameter &param, std::ostream & os = std::cout);

    /**
     * read MULTIBATH block.
     */
    void read_MULTIBATH(simulation::Parameter &param, std::ostream & os = std::cout);

    /**
     * read J-VAL block.
     */
    void read_JVALUE(simulation::Parameter &param, std::ostream & os = std::cout);

    /**
     * read PERSCALE block.
     */
    void read_PERSCALE(simulation::Parameter &param, std::ostream & os = std::cout);

    /**
     * read ROTTRANS block.
     */
    void read_ROTTRANS(simulation::Parameter &param, std::ostream & os = std::cout);

    /**
     * read REPLICA block.
     */
    void read_REPLICA(simulation::Parameter &param, std::ostream & os = std::cout);

    /**
     * read INNERLOOP block.
     */
    void read_INNERLOOP(simulation::Parameter &param, std::ostream & os = std::cout);

    /**
     * read CGRAIN block.
     */
    void read_CGRAIN(simulation::Parameter &param, std::ostream & os = std::cout);

    /**
     * read MULTICELL block.
     */
    void read_MULTICELL(simulation::Parameter & param, std::ostream & os = std::cout);

    /**
     * read READTRAJ block.
     */
    void read_READTRAJ(simulation::Parameter & param, std::ostream & os = std::cout);

    /**
     * read INTEGRATE block.
     */
    void read_INTEGRATE(simulation::Parameter & param, std::ostream & os = std::cout);
    
    /**
     * read STOCHASTIC block.
     */
    void read_STOCHASTIC(simulation::Parameter & param, std::ostream & os = std::cout);

    /**
     * read EWARN block.
     */
    void read_EWARN(simulation::Parameter & param, std::ostream & os = std::cout);

    /**
     * read MULTISTEP block.
     */
    void read_MULTISTEP(simulation::Parameter & param, std::ostream & os = std::cout);
    
    /**
     * read MONTECARLO block.
     */
    void read_MONTECARLO(simulation::Parameter & param, std::ostream & os = std::cout);
    
    /**
     * read RAND block.
     */
    void read_RAMD(simulation::Parameter & param, std::ostream & os = std::cout);
    
    /**
     * read CONSISTENCYCHECK block.
     */
    void read_CONSISTENCYCHECK(simulation::Parameter & param, std::ostream & os = std::cout);

    /**
     * read THERMOSTAT block.
     */
    void read_THERMOSTAT(simulation::Parameter & param, std::ostream & os = std::cout);
    
    /**
     * read BAROSTAT block.
     */
    void read_BAROSTAT(simulation::Parameter & param, std::ostream & os = std::cout);
    
    /**
     * read VIRIAL block.
     */
    void read_VIRIAL(simulation::Parameter & param, std::ostream & os = std::cout);
    
    /**
     * read GROMOS96COMPAT block.
     */
    void read_GROMOS96COMPAT(simulation::Parameter & param, std::ostream & os = std::cout);
    
    /**
     * read PATHINT block.
     */
    void read_PATHINT(simulation::Parameter & param, std::ostream & os = std::cout);
    
    /**
     * read NEIGHBOURLIST block.
     */
    void read_NEIGHBOURLIST(simulation::Parameter & param, std::ostream & os = std::cout);
    
    /**
     * read NONBONDED block.
     */
    void read_NONBONDED(simulation::Parameter & param, std::ostream & os = std::cout);
    
    /**
     * read LOCALELEVATION block.
     */
    void read_LOCALELEVATION(simulation::Parameter & param, std::ostream & os = std::cout);
    
    /**
     * read UMBRELLA block.
     */
    void read_UMBRELLA(simulation::Parameter & param, std::ostream & os = std::cout);
    
    /**
     * read FORCEFIELD block.
     */
    void read_FORCEFIELD(simulation::Parameter &param, std::ostream & os = std::cout);
    /**
     * read POLARIZE block.
     */
    void read_POLARIZE(simulation::Parameter &param, std::ostream & os = std::cout);
    /**
     * read RANDOMNUMBERS block.
     */
    void read_RANDOMNUMBERS(simulation::Parameter &param, std::ostream & os = std::cout);

  };
  
} // io

#endif
