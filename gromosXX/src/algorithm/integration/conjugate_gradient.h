/**
 * @file conjugate_gradient.h
 * conjugate gradient energy minimisation
 */

#ifndef INCLUDED_CONJUGATE_GRADIENT_H
#define INCLUDED_CONJUGATE_GRADIENT_H

#include "../../algorithm/constraints/shake.h"

namespace algorithm
{
  /**
   * @class Conjugate_Gradient
   * implements conjugate gradient energy minimisation
   */
  class Conjugate_Gradient : public Algorithm
  {
  public:
    /**
     * Constructor.
     */
    explicit Conjugate_Gradient(
      algorithm::Algorithm_Sequence &md_seq
      ) : Algorithm("ConjugateGradient"), cgrad_seq(md_seq){}

    /**
     * Destructor.
     */
    virtual ~Conjugate_Gradient(){}
    
    /**
     * conjugate gradient step.
     */
    virtual int apply(topology::Topology &topo, 
		      configuration::Configuration &conf,
		      simulation::Simulation &sim);

    /**
     * init an algorithm
     * print out input parameter, what it does...
     */
    virtual int init(topology::Topology & topo,
		     configuration::Configuration & conf,
		     simulation::Simulation & sim,
		     std::ostream & os = std::cout,
		     bool quiet = false);



  private:

    /**
     * MD sequence
     */
    algorithm::Algorithm_Sequence & cgrad_seq;

    /**
     * Positional constraints 
     */
    algorithm::Position_Constraints * cgrad_posres;

    /**
     * SHAKE algorithm
     */
    algorithm::Shake * cgrad_shake;

    /**
     * Calculation of interactions
     */
    interaction::Forcefield * cgrad_ff;

    /**
     * Do SHAKE within Conjugate Gradient step
     */
    bool do_shake;

    /**
     * Apply positional restraints within Conjugate Gradient step
     */
    bool do_posres;

    /**
     * Minimum from cubic interpolation along the search direction
     */
    double Smin;
    
    /**
    * Squared magnitude of search direction
    * used in the step size and with the SHAKE algorithm
    */
    double p_squared;

    /**
    * Total number of performed energy calculations
    */
    unsigned total_iterations;

    /**
    * Configuration to use with SHAKE to not mess up the main configuration
    */
    configuration::Configuration conf_sh;

    /**
     * Initialise separate configuration for SHAKE
     */
    inline void conf_sh_init(const configuration::Configuration & conf) {
      conf_sh = conf;
      conf_sh.old().constraint_force = conf_sh.current().constraint_force = conf.old().constraint_force;
    }
    
    /**
     * Shake old search directions
     */
    inline void shake_cgrads(const topology::Topology & topo, configuration::Configuration & conf) {
      for(unsigned int i=0; i<topo.num_atoms(); ++i) {
        conf.old().cgrad(i) = (conf.current().pos(i) - conf.old().pos(i)) / Smin;
      }
    }

    /**
     * Shake current forces
     */
    int shake_forces(topology::Topology & topo,
         configuration::Configuration & conf,
         simulation::Simulation & sim,
         double shake_step);

    /**
     * Calculate the search direction coefficient
     */
    double calculate_beta(const topology::Topology & topo,
         const configuration::Configuration & conf,
         const simulation::Simulation & sim);
    
    /**
     * Update search directions and return the gradient
     */
    double calculate_cgrad(const topology::Topology & topo,
         configuration::Configuration & conf,
         const double & beta);
    
    /**
    * Calculate interactions and energies of conformation
    * Optionally also apply constraints
    */
    int evaluate_configuration(topology::Topology & topo,
         configuration::Configuration & conf,
         simulation::Simulation & sim,
         double & ene,
         bool do_pos_shake);
  
  };
} // algorithm

#endif

