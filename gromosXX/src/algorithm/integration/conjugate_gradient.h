/**
 * @file conjugate_gradient.h
 * conjugate gradient energy minimisation
 */

#ifndef INCLUDED_CONJUGATE_GRADIENT_H
#define INCLUDED_CONJUGATE_GRADIENT_H

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
    Conjugate_Gradient(interaction::Forcefield *ff) : Algorithm("ConjugateGradient"), cg_ff(*ff) {}

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
     * Forcefield
     */
    interaction::Forcefield & cg_ff;

    /**
     * Calculates an old search direction coefficient
     */
    double calculate_beta(const topology::Topology & topo,
         const configuration::Configuration & conf,
         const simulation::Simulation & sim);
    
    /**
     * Updates search directions and returns sum of their squared sizes
     */
    double calculate_cgrad(const topology::Topology & topo,
         configuration::Configuration & conf,
         const double & beta);

  };
  
} // algorithm

#endif

