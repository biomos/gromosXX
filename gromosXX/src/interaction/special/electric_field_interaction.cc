/**
 * @file electric_field_interaction.cc
 * template external Electric_Field_Interaction
 */


#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"
#include "../../interaction/interaction.h"

#include "../../math/periodicity.h"

// special interactions
#include "../../interaction/interaction_types.h"

#include "../../interaction/special/electric_field_interaction.h"

#include "../../util/template_split.h"
#include "../../util/debug.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE special

/**
 * calculate electric field interactions
 */

int interaction::Electric_Field_Interaction::
calculate_interactions(topology::Topology& topo,
                       configuration::Configuration& conf,
                       simulation::Simulation& sim)
{
  // loop over the atoms
  // unit E = e . nm^-2
  math::Vec E(sim.param().electric.Ef_x, sim.param().electric.Ef_y, sim.param().electric.Ef_z);
  switch (sim.param().force.interaction_function) {
    case simulation::lj_crf_func:
    case simulation::pol_lj_crf_func:
    case simulation::pol_off_lj_crf_func:
    {
      // electric field in a spherical cavity is higher!
      E *= 3 * sim.param().nonbonded.rf_epsilon / (2 * sim.param().nonbonded.rf_epsilon
              + sim.param().nonbonded.epsilon);
      break;
    }
    case simulation::cggromos_func:
    {
      E *= 3 * sim.param().nonbonded.rf_epsilon / (2 * sim.param().nonbonded.rf_epsilon
              + sim.param().cgrain.EPS);
      break;
    }
    default:
      io::messages.add("Electric_Field_Interaction",
              "interaction function not implemented",
              io::message::critical);
  }

  if (simulation::pol_off_lj_crf_func != 0)
    for (unsigned int i = 0; i < topo.num_atoms(); ++i){
    // math::four_pi_eps_i contains already epsilon of cutoff-sphere (param().nonbonded.epsilon)
     math::Vec force = (math::four_pi_eps_i*sim.param().nonbonded.epsilon) * topo.charge(i) * E;
     conf.current().force(i) +=(1-topo.gamma(i))*force;

     if(topo.gamma(i)!=0.0){
         conf.current().force(topo.gamma_j(i)) +=topo.gamma(i)/2*force;
         conf.current().force(topo.gamma_k(i)) +=topo.gamma(i)/2*force;
     }
  } 
  else
    for (unsigned int i = 0; i < topo.num_atoms(); ++i){
    // math::four_pi_eps_i contains already epsilon of cutoff-sphere (param().nonbonded.epsilon)
      conf.current().force(i) += (math::four_pi_eps_i*sim.param().nonbonded.epsilon) * topo.charge(i) * E;
    }

  //The electric energies were not added to the total energy
  /* F = qE
   * U = - integral(Fdr) => U = -qEr
   * The energy depends on the (absolute) position of the particle
   * Maybe it doesn't make sense to use this energy within
   * periodic boundary conditions
   */

  

  return 0;
}

