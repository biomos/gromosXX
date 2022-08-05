/**
 * @file dfunct_interaction.cc
 * dfunct
 * 
 */

#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"
#include "../../interaction/interaction.h"

#include "../../math/periodicity.h"
#include "../../math/gmath.h"
#include "../../interaction/interaction_types.h"

#include "../../interaction/special/dfunct_interaction.h"

#include "../../util/template_split.h"
#include "../../util/debug.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE special

template<math::boundary_enum B>
static int _calculate_dfunct_substitution_form(topology::Topology& topo, 
																					     configuration::Configuration& conf, 
																					     simulation::Simulation& sim) {
	// shorten the code 
	int atom_i = sim.param().dfunct.atom_i; 
	int atom_j = sim.param().dfunct.atom_j;
	int atom_k = sim.param().dfunct.atom_k;
	int atom_l = sim.param().dfunct.atom_l;
	double r_0 = sim.param().dfunct.r_0;
	int d = sim.param().dfunct.d;
	double force = sim.param().dfunct.force;
	DEBUG(10, "DFUNCT Calculating subtitution type potential")
	DEBUG(10, "DFUNCT atom_i " << atom_i << math::v2s(conf.current().pos(atom_i)));
	DEBUG(10, "DFUNCT atom_j " << atom_j << math::v2s(conf.current().pos(atom_j)));
	DEBUG(10, "DFUNCT atom_k " << atom_k << math::v2s(conf.current().pos(atom_k)));
	DEBUG(10, "DFUNCT atom_l " << atom_l << math::v2s(conf.current().pos(atom_l)));
	math::Vec dist_vec_ji, dist_vec_lk;
	math::Periodicity<B> periodicity(conf.current().box);
	periodicity.nearest_image(conf.current().pos(atom_j), conf.current().pos(atom_i), dist_vec_ji);
	periodicity.nearest_image(conf.current().pos(atom_l), conf.current().pos(atom_k), dist_vec_lk);
	double dist_ji   = math::abs(dist_vec_ji);
	double dist_lk   = math::abs(dist_vec_lk);
	DEBUG(30, "DFUNCT dist_vec_ji " << math::v2s(dist_vec_ji));
	DEBUG(30, "DFUNCT dist_vec_lk " << math::v2s(dist_vec_lk));
	DEBUG(30, "DFUNCT dist_ji " << dist_ji);
	DEBUG(30, "DFUNCT dist_lk " << dist_lk);

	// compute forces on atoms i, j, k, l and the combined biasing potential
	math::Vec prefactor_i =     force * (dist_vec_ji / dist_ji);
	math::Vec prefactor_k = d * force * (dist_vec_lk / dist_lk);
	double force_term = (dist_ji + d * dist_lk - r_0);
	math::Vec force_i =  prefactor_i  * force_term;
	math::Vec force_j = -force_i;
	math::Vec force_k =  prefactor_k * force_term;
	math::Vec force_l = -force_k;
	double V_bias = 0.5 * force * (dist_ji + d * dist_lk - r_0) * (dist_ji + d * dist_lk - r_0);
	DEBUG(30, "DFUNCT Prefactor i " << math::v2s(prefactor_i));
	DEBUG(30, "DFUNCT Prefactor k " << math::v2s(prefactor_k));
	DEBUG(10, "DFUNCT Force on i " << math::v2s(force_i));
  DEBUG(10, "DFUNCT Force on j " << math::v2s(force_j));
  DEBUG(10, "DFUNCT Force on k " << math::v2s(force_k));
  DEBUG(10, "DFUNCT Force on l " << math::v2s(force_l));
	DEBUG(10, "DFUNCT V_bias " << V_bias);
	
	// store forces
	conf.current().force(atom_i) += force_i;
	conf.current().force(atom_j) += force_j;
	conf.current().force(atom_k) += force_k;
	conf.current().force(atom_l) += force_l;

	// distribute potential over all participating atoms
	conf.current().energies.distanceres_energy[topo.atom_energy_group()[atom_i]] += 0.25 * V_bias;
	conf.current().energies.distanceres_energy[topo.atom_energy_group()[atom_j]] += 0.25 * V_bias;
	conf.current().energies.distanceres_energy[topo.atom_energy_group()[atom_k]] += 0.25 * V_bias;
	conf.current().energies.distanceres_energy[topo.atom_energy_group()[atom_l]] += 0.25 * V_bias;
	return 0;
}

template<math::boundary_enum B>
static int _calculate_dfunct_cycloaddition_form(topology::Topology& topo, 
																					     configuration::Configuration& conf, 
																					     simulation::Simulation& sim) {
	// shorten the code 
	int atom_i = sim.param().dfunct.atom_i; 
	int atom_j = sim.param().dfunct.atom_j;
	int atom_k = sim.param().dfunct.atom_k;
	int atom_l = sim.param().dfunct.atom_l;
	double r_0 = sim.param().dfunct.r_0;
	int d = sim.param().dfunct.d;
	double force = sim.param().dfunct.force;
	DEBUG(10, "DFUNCT Calculating cycloaddition-type potential")
	DEBUG(10, "DFUNCT atom_i " << atom_i << math::v2s(conf.current().pos(atom_i)));
	DEBUG(10, "DFUNCT atom_j " << atom_j << math::v2s(conf.current().pos(atom_j)));
	DEBUG(10, "DFUNCT atom_k " << atom_k << math::v2s(conf.current().pos(atom_k)));
	DEBUG(10, "DFUNCT atom_l " << atom_l << math::v2s(conf.current().pos(atom_l)));
	math::Vec dist_vec_ji, dist_vec_lk, dist_vec_lkji;
	math::Periodicity<B> periodicity(conf.current().box);
	// algorithm:
	// (1) find mean from i to j and k to l, respectively
	// (2) add r_i and r_k, respectively
	// (3) find mean between these new vectors
	// note: this procedure can be simplified as demonstrated below
	periodicity.nearest_image(conf.current().pos(atom_j), -1.0 * conf.current().pos(atom_i), dist_vec_ji);
	periodicity.nearest_image(conf.current().pos(atom_l), -1.0 * conf.current().pos(atom_k), dist_vec_lk);
	double dist_ji = math::abs(dist_vec_ji);
	double dist_lk = math::abs(dist_vec_lk);
	DEBUG(30, "DFUNCT dist_vec_ji " << math::v2s(dist_vec_ji));
	DEBUG(30, "DFUNCT dist_vec_lk " << math::v2s(dist_vec_lk));
	DEBUG(30, "DFUNCT dist_ji " << dist_ji);
	DEBUG(30, "DFUNCT dist_lk " << dist_lk);
	math::Vec dist_vec_lk_halfs = 0.5 * dist_vec_lk;
	math::Vec dist_vec_ji_halfs = 0.5 * dist_vec_ji;
	periodicity.nearest_image(dist_vec_lk_halfs, dist_vec_ji_halfs, dist_vec_lkji);
	double dist_ljki = math::abs(dist_vec_lkji);
	DEBUG(30, "DFUNCT dist_vec_lkji " << math::v2s(dist_vec_lkji));
	DEBUG(30, "DFUNCT dist_ljki " << dist_ljki);
	
	// compute forces
	math::Vec force_i = force * (dist_vec_lkji / dist_ljki) * (dist_ljki - r_0);
	math::Vec force_j = force_i;
	math::Vec force_k = -1.0 * force_i;
	math::Vec force_l = force_k;

	// compute potential
	double V_bias = 0.5 * force * (dist_ljki - r_0) * (dist_ljki - r_0);

	// store forces
	conf.current().force(atom_i) += force_i;
	conf.current().force(atom_j) += force_j;
	conf.current().force(atom_k) += force_k;
	conf.current().force(atom_l) += force_l;

	DEBUG(10, "DFUNCT Force on i " << math::v2s(force_i));
  DEBUG(10, "DFUNCT Force on j " << math::v2s(force_j));
  DEBUG(10, "DFUNCT Force on k " << math::v2s(force_k));
  DEBUG(10, "DFUNCT Force on l " << math::v2s(force_l));
	DEBUG(10, "DFUNCT V_bias " << V_bias);

	// distribute potential over all participating atoms
	conf.current().energies.distanceres_energy[topo.atom_energy_group()[atom_i]] += 0.25 * V_bias;
	conf.current().energies.distanceres_energy[topo.atom_energy_group()[atom_j]] += 0.25 * V_bias;
	conf.current().energies.distanceres_energy[topo.atom_energy_group()[atom_k]] += 0.25 * V_bias;
	conf.current().energies.distanceres_energy[topo.atom_energy_group()[atom_l]] += 0.25 * V_bias;
	return 0;
}

int interaction::DFunct_Interaction::init(topology::Topology& topo,
		     																  configuration::Configuration& conf,
		     																  simulation::Simulation& sim,
		     																  std::ostream& os,
		     																  bool quiet) {
  
	return 0;
}

int interaction::DFunct_Interaction::calculate_interactions(topology::Topology & topo,
				                                                    configuration::Configuration& conf,
				                                                    simulation::Simulation& sim) {
  // atomic distances expressed as vectors
	// in GROMOS vector r_ji is defined as the vector from point i to point j (r_j - r_i)
	// find nearest periodic copies
	switch (sim.param().dfunct.dfunct) {
		case simulation::dfunct_substitution:
		  // restrain the distances of an incoming (i) and outgoing (l) atom to a central atom (j,k), respectively
			SPLIT_BOUNDARY(_calculate_dfunct_substitution_form, topo, conf, sim);
			break;

		case simulation::dfunct_cycloaddition:
		  // restrain the distance between two cycloaddition reaction partners (i, j) and (k, l), respectively
			SPLIT_BOUNDARY(_calculate_dfunct_cycloaddition_form, topo, conf, sim);
			break;

		default:
      io::messages.add("DFUNCT functional not implemented", "DFunct_Interaction", io::message::critical);
      break;
		}
	
	return 0;
}