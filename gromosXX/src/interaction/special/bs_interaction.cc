#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"
#include "../../interaction/interaction.h"

#include "../../math/periodicity.h"

// special interactions
#include "../../util/bs_umbrella.h"
#include "bs_interaction.h"

#include "../../util/template_split.h"
#include "../../util/debug.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE special

int interaction::BS_Interaction::init(topology::Topology& topo, 
                                      configuration::Configuration& conf, 
                                      simulation::Simulation& sim, 
                                      std::ostream& os, 
                                      bool quiet)
{
  DEBUG(4, "Initialze the B&S-LEUS interaction.");
  if(!quiet){
    os << "BSLEUS\n";
    os << conf.special().bs_umbrella.str();
    os << "END\n";
  }
  return 0;
}

int interaction::BS_Interaction::calculate_interactions(
                                        topology::Topology& topo, 
                                        configuration::Configuration& conf, 
                                        simulation::Simulation& sim)
{
  m_timer.start();
  conf.special().bs_umbrella.apply(conf, sim);
  m_timer.stop();
  return 0;
}