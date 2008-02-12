/**
 * @file vgrid_pairlist_algorithm.cc
 * vector grid pairlist algorithm
 */

#include <stdheader.h>

#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <simulation/simulation.h>
#include <configuration/configuration.h>

#include <interaction/interaction_types.h>

#include <interaction/nonbonded/pairlist/pairlist.h>
#include <interaction/nonbonded/interaction/storage.h>
#include <interaction/nonbonded/interaction/nonbonded_parameter.h>

#include <interaction/nonbonded/pairlist/pairlist_algorithm.h>
#include <interaction/nonbonded/pairlist/vgrid_pairlist_algorithm.h>

#include <interaction/nonbonded/vector/vgrid.h>

#include <util/debug.h>

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE pairlist

interaction::VGrid_Pairlist_Algorithm::VGrid_Pairlist_Algorithm()
  : interaction::Pairlist_Algorithm(),
    m_solvent_solvent_timing(0.0),
    m_spc_timing(0.0)
{
}
/**
 * initialize
 */
int interaction::VGrid_Pairlist_Algorithm::
init(topology::Topology &topo,
        configuration::Configuration &conf,
        simulation::Simulation &sim,
        std::ostream &os,
        bool quiet)
{
  if (!quiet)
    os << "vector grid pairlist algorithm\n";
  if (conf.boundary_type != math::rectangular){
    io::messages.add("VGrid Pairlist Algorithm",
                     "only implemented for rectangular boundary conditions?",
                     io::message::error);
    return 1;
  }

  if (sim.param().perturbation.perturbation ||
      sim.param().pairlist.atomic_cutoff){
    io::messages.add("VGrid pairlist algorithm: no perturbation / atomic cutoff allowed",
		     "vgrid_pairlist_algorithm",
		     io::message::error);
  }
  return 0;
};

/**
 * prepare the pairlist
 */
int interaction::VGrid_Pairlist_Algorithm::
prepare(topology::Topology & topo,
	configuration::Configuration & conf,
	simulation::Simulation & sim)
{
  DEBUG(7, "vgrid pairlist algorithm : prepare");

  grid_prepare(topo, conf, sim, *m_param);

  return 0;
}

void interaction::VGrid_Pairlist_Algorithm::
update(topology::Topology & topo,
       configuration::Configuration & conf,
       simulation::Simulation & sim,
       interaction::Storage & storage,
       interaction::Pairlist & pairlist,
       unsigned int begin, unsigned int end,
       unsigned int stride)
{
  grid_update(topo, conf, sim);
  m_timing += grid_timing;
}

void interaction::VGrid_Pairlist_Algorithm::
update_perturbed(topology::Topology & topo,
		 configuration::Configuration & conf,
		 simulation::Simulation & sim,
		 interaction::Storage & storage,
		 interaction::Pairlist & pairlist,
		 interaction::Pairlist & perturbed_pairlist,
		 unsigned int begin, unsigned int end,
		 unsigned int stride)
{
  io::messages.add("no perturbed vgrid pairlist implemented",
		   "vgrid_pairlist_algorithm",
		   io::message::error);
}
