/**
 * @file perturbation_md.h
 * md algorithm including perturbation.
 * should maybe be changed to a mix in class
 * using multiple inheritance???
 */

#ifndef INCLUDED_PERTURBATION_MD_H
#define INCLUDED_PERTURBATION_MD_H

namespace algorithm
{
  /**
   * @class Perturbation_MD
   * MD algorithm with perturbation.
   */
  template<typename t_md_spec = perturbed_MD_spec<>,
	   typename t_interaction_spec = Interaction_spec<
    typename t_md_spec::simulation_type,
    true>
  >
  class Perturbation_MD : public MD<t_md_spec, t_interaction_spec>
  {
  public:
    /**
     * Constructor.
     */
    Perturbation_MD();
    
  protected:
    /**
     * initialize the input.
     */
    virtual void init_input(io::Argument &args, io::InTopology &topo,
			    io::InTrajectory &sys, io::InInput &input);
    
    /**
     * read the input and setup a standard simulation.
     */
    virtual void read_input(io::Argument &args, io::InTopology &topo,
			    io::InTrajectory &sys, io::InInput &input);
    
    virtual void G96Forcefield(io::InTopology &topo,
			       io::InInput &input,
			       io::Argument &args);
    
    /**
     * initialize the perturbation parameters.
     */
    int init_perturbation(io::Argument &args, io::InInput &input);
    
    /**
     * print pairlists.
     */
    virtual void print_pairlists(); 
    
    /**
     * energy calculation.
     */
    void do_perturbed_energies();
    
    /**
     * post step
     */
    virtual void post_step();
    
    /**
     * pre md.
     */
    virtual int pre_md(io::InInput &input);
    
    /**
     * post md
     */
    virtual void post_md();
    
    /**
     * store the lambda change per step.
     */
    double m_d_lambda;
    
  };
  
} // algorithm

#include "perturbation_md.tcc"

#endif
