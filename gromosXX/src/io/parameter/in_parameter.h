/**
 * @file in_parameter.h
 * read in a md++ input file.
 */
/**
 * @page input input file format
 *
 * @date 28-10-2008
 *
 * - @ref  title
 * - @ref  system
 * - @ref  initialise
 * - @ref  step
 * - @ref  stochdyn
 * - @ref  energymin
 * - @ref  boundcond
 * - @ref  multicell
 * - @ref  multibath
 * - @ref  pressurescale
 * - @ref  comtransrot
 * - @ref  printout
 * - @ref  writetraj
 * - @ref  constraint
 * - @ref  force
 * - @ref  covalentform
 * - @ref  hoomd
 * - @ref  pairlist
 * - @ref  nonbonded
 * - @ref  positionres
 * - @ref  posres
 * - @ref  perturbation
 * - @ref  precalclam
 * - @ref  pttopo
 * - @ref  jval
 * - @ref  jvalue
 * - @ref  perscale
 * - @ref  replica
 * - @ref  readtraj
 * - @ref  integrate
 * - @ref  cgrain
 * - @ref  rottrans
 * - @ref  distanceres
 * - @ref  distancefield
 * - @ref  disres
 * - @ref  angleres
 * - @ref  dihedralres
 * - @ref  dihrest
 * - @ref  multistep
 * - @ref  ewarn
 * - @ref  polarise
 * - @ref  randomnumbers
 * - @ref  spc_loops
 * - @ref  EDS
 * - @ref  REEDS
 * - @ref  AEDS
 * - @ref  LAMBDAS
 * - @ref  localelev
 * - @ref  bsleusparam
 * - @ref  electric
 * - @ref  sasa
 * - @ref  symres
 * - @ref  nemd
 * - @ref  multigradient
 * - @ref  addecouple
 * - @ref  orderparamres
 * - @ref  rdcres
 * - @ref  qmmm
 * - @ref  xrayres
 * - @ref  amber
 * - @ref  dfunct
 */



#ifndef INCLUDED_IN_PARAMETER_H
#define INCLUDED_IN_PARAMETER_H

#include "../instream.h"

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
     * read energymin block.
     */
    void read_ENERGYMIN(simulation::Parameter &param, std::ostream & os = std::cout);

    /**
     * read step block.
     */
    void read_STEP(simulation::Parameter &param, std::ostream & os = std::cout);

    /**
     * read the SYSTEM block.
     */
    void read_SYSTEM(simulation::Parameter &param, std::ostream & os = std::cout);

    /**
     * read the CONSTRAINT block.
     */
    void read_CONSTRAINT(simulation::Parameter &param, std::ostream & os = std::cout);

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
     * read PRINTOUT block.
     */
    void read_PRINTOUT(simulation::Parameter &param, std::ostream & os = std::cout);

    /**
     * read WRITETRAJ block.
     */
    void read_WRITETRAJ(simulation::Parameter &param, std::ostream & os = std::cout);

    /**
     * read the PRESSURESCALE block.
     */
    void read_PRESSURESCALE(simulation::Parameter &param, std::ostream & os = std::cout);

    /**
     * read the BOUNDCOND block.
     */
    void read_BOUNDCOND(simulation::Parameter &param, std::ostream & os = std::cout);

    /**
     * read the PERTURBATION block.
     */
    void read_PERTURBATION(simulation::Parameter &param, std::ostream & os = std::cout);

    /**
     * read INITIALISE block.
     */
    void read_INITIALISE(simulation::Parameter &param, std::ostream & os = std::cout);

    /**
     * read COMTRANSROT block.
     */
    void read_COMTRANSROT(simulation::Parameter &param, std::ostream & os = std::cout);

    /**
     * read HOOMD block.
     */
    void read_HOOMD(simulation::Parameter &param, std::ostream & os = std::cout);

	/**
     * read PAIRLIST block.
     */
    void read_PAIRLIST(simulation::Parameter &param, std::ostream & os = std::cout);

    /**
     * read POSRES block.
     */
    void read_POSITIONRES(simulation::Parameter &param, std::ostream & os = std::cout);

    /**
     * read XRAYRES block.
     */
    void read_XRAYRES(simulation::Parameter &param, std::ostream & os = std::cout);

    /**
     * read DISTANCERES block.
     */
    void read_DISTANCERES(simulation::Parameter &param, std::ostream & os = std::cout);

    /**
     * read DISTANCEFIELD block.
     */
    void read_DISTANCEFIELD(simulation::Parameter &param, std::ostream & os = std::cout);

    /**
     * read ANGLERES block.
     */
    void read_ANGLERES(simulation::Parameter &param, std::ostream & os = std::cout);

    /**
     * read DIHEDRALRES block.
     */
    void read_DIHEDRALRES(simulation::Parameter &param, std::ostream & os = std::cout);

    /**
     * read MULTIBATH block.
     */
    void read_MULTIBATH(simulation::Parameter &param, std::ostream & os = std::cout);

    /**
     * read JVALUERES block.
     */
    void read_JVALUERES(simulation::Parameter &param, std::ostream & os = std::cout);
    /**
     * read ORDERPARAMRES block.
     */
    void read_ORDERPARAMRES(simulation::Parameter &param, std::ostream & os = std::cout);

    /**
     * read RDCRES block.
     */
    void read_RDCRES(simulation::Parameter &param, std::ostream & os = std::cout);

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
    * read REPLICA_EDS block.
    */
    void read_REPLICA_EDS(simulation::Parameter &param, std::ostream & os = std::cout);

     
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
     * read STOCHDYN block.
     */
    void read_STOCHDYN(simulation::Parameter & param, std::ostream & os = std::cout);

    /**
     * read EWARN block.
     */
    void read_EWARN(simulation::Parameter & param, std::ostream & os = std::cout);

    /**
     * read MULTISTEP block.
     */
    void read_MULTISTEP(simulation::Parameter & param, std::ostream & os = std::cout);

    /**
     * read CHEMICALMONTECARLO block.
     */
    void read_CHEMICALMONTECARLO(simulation::Parameter & param, std::ostream & os = std::cout);

    /**
     * read RAND block.
     */
    void read_RAMD(simulation::Parameter & param, std::ostream & os = std::cout);

    /**
     * read POLARISE block.
     */
    void read_POLARISE(simulation::Parameter &param, std::ostream & os = std::cout);

    /**
     * read RANDOMNUMBERS block.
     */
    void read_RANDOMNUMBERS(simulation::Parameter &param, std::ostream & os = std::cout);

    /**
     * read LAMBDAS block.
     */
    void read_LAMBDAS(simulation::Parameter &param, std::ostream & os = std::cout);
    
    /** 
     * read PRECALCLAM block.
     */
    void read_PRECALCLAM(simulation::Parameter &param, std::ostream & os = std::cout);

    /**
     * read blocks that are either g96 or promd specific and tell user what to do
     */
    void read_known_unsupported_blocks();
    /**
     * read EDS block.
     */
    void read_EDS(simulation::Parameter &param, std::ostream & os = std::cout);
    /**
    * read AEDS block.
    */
    void read_AEDS(simulation::Parameter &param, std::ostream & os = std::cout);
    /**
     * read NONBONDED block.
     */
    void read_NONBONDED(simulation::Parameter &param, std::ostream & os = std::cout);
    /**
     * read LOCALELEV block.
     */
    void read_LOCALELEV(simulation::Parameter &param, std::ostream & os = std::cout);
    /**
     * read the BSLEUS block
     */
    void read_BSLEUS(simulation::Parameter &param, std::ostream &os = std::cout);
    /**
     * read ELECTRIC block.
     */
    void read_ELECTRIC(simulation::Parameter &param, std::ostream & os = std::cout);
    /**
     * read SASA block.
     */
    void read_SASA(simulation::Parameter &param, std::ostream & os = std::cout);
    /**
     * read NEMD block.
     */
    void read_NEMD(simulation::Parameter &param, std::ostream & os = std::cout);
    /**
     * read MULTIGRADIENT block.
     */
    void read_MULTIGRADIENT(simulation::Parameter &param, std::ostream & os = std::cout);
     /**
     * read ADDECOUPLE block.
     */
    void read_ADDECOUPLE(simulation::Parameter &param, std::ostream & os = std::cout);
    /**
     * read QMMM block
     */
    void read_QMMM(simulation::Parameter &param, std::ostream & os = std::cout);
    /**
     * read SYMRES block
     */
    void read_SYMRES(simulation::Parameter &param, std::ostream & os = std::cout);
    /**
     * read AMBER block
     */
    void read_AMBER(simulation::Parameter &param, std::ostream & os = std::cout);
    /**
     * read DFUNCT block 
     */
    void read_DFUNCT(simulation::Parameter &param, std::ostream & os = std::cout);
  };

} // io

#endif
