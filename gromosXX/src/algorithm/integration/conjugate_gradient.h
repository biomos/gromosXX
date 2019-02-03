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
    Conjugate_Gradient(
      algorithm::Algorithm_Sequence &md_seq
      ) : Algorithm("ConjugateGradient"),
          cgrad_seq(md_seq){}

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

    bool do_shake;
    bool do_posres;

    /**
     * Calculate the search direction coefficient
     */
    double calculate_beta(const topology::Topology & topo,
         const configuration::Configuration & conf,
         const simulation::Simulation & sim);
    
    /**
     * Update search directions and return sum of their squared sizes
     */
    double calculate_cgrad(const topology::Topology & topo,
         configuration::Configuration & conf,
         const double & beta);
    
    /**
    * Calculate interactions and energies of conformation
    * Optionally also apply constraints
    */
    int evaluate_conf(topology::Topology & topo,
         configuration::Configuration & conf,
         simulation::Simulation & sim);
  
  };
} // algorithm

#endif

