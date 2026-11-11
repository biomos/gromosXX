/**
 * @file colvar_restraint_interaction.h
 * collective variable restraining
 */

#ifndef INCLUDED_COLVAR_RESTRAINT_INTERACTION_H
#define INCLUDED_COLVAR_RESTRAINT_INTERACTION_H

#include <vector>

#include "../../interaction/interaction.h"
#include "../../simulation/parameter.h"
#include "../../interaction/special/colvar/colvar.h"
#include "../../interaction/special/colvar/colvar_bias.h"

namespace interaction
{
  /**
   * @class colvar_restraint_interaction
   * calculates the colvar restraining interaction
   */
  class Colvar_Restraint_Interaction : public Interaction
  {
  public:
    /** Constructor. */
    Colvar_Restraint_Interaction() : Interaction("ColvarRestraint"), sum(false),
              Ctot0(0), m_average_step(-1) {}

    /** Destructor. */
    virtual ~Colvar_Restraint_Interaction();

    /** init */
    virtual int init(topology::Topology &topo,
             configuration::Configuration &conf,
             simulation::Simulation &sim,
             std::ostream &os = std::cout,
             bool quiet = false);

    /** calculate the interactions. */
    virtual int calculate_interactions(topology::Topology & topo,
                       configuration::Configuration & conf,
                       simulation::Simulation & sim);

  private:
    /** store pointers to all specified colvars */
    std::vector<Colvar *> m_colvars;

    /** one reusable restraint/bias object per colvar */
    std::vector<Colvar_Bias *> m_biases;

    /** type name matching each colvar */
    std::vector<std::string> m_colvar_types;
    std::vector<unsigned int> m_spec_ordinals;

    /** put the bias on the sum of all collective variables */
    bool sum;

    /** sum of all target values */
    double Ctot0;
    long long m_average_step;

    /** apply forces and calculate potential for given colvar */
    double apply_restraint(topology::Topology & topo,
      configuration::Configuration & conf,
      simulation::Simulation & sim,
      const Colvar &cv,
      const std::string &type,
      Colvar_Bias &bias, bool advance_average);

    /** construct reusable bias settings from one COLVARRES line */
    Colvar_Bias::Settings settings_from_spec(
      const simulation::Parameter::colvar_bias_spec &spec,
      simulation::Simulation &sim) const;
  };

} // interaction

#endif
