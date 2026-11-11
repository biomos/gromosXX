/** State-resolved CV biases accumulated before EDS mixing. */
#ifndef INCLUDED_EDS_COLVAR_RESTRAINT_INTERACTION_H
#define INCLUDED_EDS_COLVAR_RESTRAINT_INTERACTION_H
#include "../../interaction/interaction.h"
namespace interaction {
class Eds_Colvar_Restraint_Interaction : public Interaction {
public:
  Eds_Colvar_Restraint_Interaction() : Interaction("EdsColvarRestraint") {}
  ~Eds_Colvar_Restraint_Interaction();
  int init(topology::Topology &, configuration::Configuration &, simulation::Simulation &,
           std::ostream &os = std::cout, bool quiet = false);
  int calculate_interactions(topology::Topology &, configuration::Configuration &,
                             simulation::Simulation &);
private:
  struct Term;
  std::vector<Term*> m_terms;
};
}
#endif
