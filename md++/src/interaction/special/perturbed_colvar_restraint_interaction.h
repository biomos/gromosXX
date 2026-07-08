/**
 * @file perturbed_colvar_restraint_interaction.h
 * perturbed collective variable restraining
 */

#ifndef INCLUDED_PERTURBED_COLVAR_RESTRAINT_INTERACTION_H
#define INCLUDED_PERTURBED_COLVAR_RESTRAINT_INTERACTION_H

#include "../../interaction/special/colvar/colvar.h"
#include "../../interaction/special/colvar/colvar_bias.h"
#include "../../simulation/parameter.h"

namespace interaction
{
  /**
   * @class perturbed_colvar_restraint_interaction
   * calculates the perturbed colvar restraining interaction
   */
  class Perturbed_Colvar_Restraint_Interaction : public Interaction
  {
  public:
    /**
     * Constructor.
     */
    Perturbed_Colvar_Restraint_Interaction() : Interaction("PerturbedColvarRestraint"), sum(false),
              Ctot0(0) {}
    
    /**
     * Destructor.
     */
    virtual ~Perturbed_Colvar_Restraint_Interaction();

    /**
     * init
     */
    virtual int init(topology::Topology &topo, 
                     configuration::Configuration &conf,
                     simulation::Simulation &sim,
                     std::ostream &os = std::cout,
                     bool quiet = false);
    /**
     * calculate the interactions.
     */
    virtual int calculate_interactions(topology::Topology & topo,
                                       configuration::Configuration & conf,
                                       simulation::Simulation & sim);
    
  private:
    struct Term {
      Term() : cv(NULL), bias(NULL) {}
      Colvar *cv;
      Colvar_Bias *bias;
      simulation::Parameter::colvar_bias_spec A;
      simulation::Parameter::colvar_bias_spec B;
      std::string type;
    };

    /**
     * store all specified perturbed colvars
     */  
     std::vector<Term> m_terms;

    /**
     * put the bias on the sum of all collective variables
     */  
     bool sum;
     
    /**
     * sum of all target values
     */  
     double Ctot0;
     
    /** apply forces and calculate potential for given colvar
     * using the biasing function specified in the COLVARRES block
     */
     double apply_restraint(topology::Topology & topo,
      configuration::Configuration & conf, simulation::Simulation & sim,
      Term &term, double lambda, double lambda_derivative,
      unsigned int energy_group);

     Colvar_Bias::Settings settings_from_spec(
      const simulation::Parameter::colvar_bias_spec &spec,
      simulation::Simulation &sim) const;
    
  };
  
} // interaction

#endif
