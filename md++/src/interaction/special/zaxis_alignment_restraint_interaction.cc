/**
 * @file zaxis_orientation_bias_interaction.cc
 * template methods of Zaxis_Orientation_Bias_Interaction
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

#include "../../interaction/special/zaxis_orientation_bias_interaction.h"

#include "../../util/template_split.h"
#include "../../util/debug.h"

#include <stdlib.h>

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE special

/**
 * calculate z-axis orientation bias interactions
 */
template<math::boundary_enum B, math::virial_enum V>
static int _calculate_zaxis_orientation_bias_interactions
(topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation & sim, double exponential_term,
std::map<int,math::Vec> &rah_map, int &error)
{
  // loop over the z-axis orientation bias
  std::vector<topology::zaxisori_restraint_struct>::const_iterator
    it = topo.zaxisori_restraints().begin(),
    to = topo.zaxisori_restraints().end();

 std::vector<double>::iterator ave_it = conf.special().zaxisoribias.av.begin();
 std::vector<double>::iterator ene_it = conf.special().zaxisoribias.energy.begin();
 std::vector<double>::iterator d_it = conf.special().zaxisoribias.d.begin();

  // math::VArray &pos   = conf.current().pos;
  // math::VArray &force = conf.current().force;
  math::Vec v, f;
  double energy;

  math::Periodicity<B> periodicity(conf.current().box);

  // for(int i=0; it != to; ++it, ++i){
  for(; it != to; ++it, ++ave_it, ++ene_it, ++d_it) {
    periodicity.nearest_image(it->v1.pos(conf,topo), it->v2.pos(conf,topo), v);

#ifndef NDEBUG
    for(int i=0;i<it->v1.size();i++){
      DEBUG(10, "v1 (atom " << i+1 << "): " << it->v1.atom(i)+1);
    }
    for(int i=0;i<it->v2.size();i++){
      DEBUG(10, "v2 (atom " << i+1 << "): " << it->v2.atom(i)+1);
    }
#endif
    DEBUG(10, "pos(v1) = " << math::v2s(it->v1.pos(conf,topo)));
    DEBUG(10, "pos(v2) = " << math::v2s(it->v2.pos(conf,topo)));

    DEBUG(9, "ZAXISORIBIAS v : " << math::v2s(v));
    int rah=it->rah;
    bool found=false;
    if(abs(rah)<=1) found=true;
    if(!found){
       std::ostringstream msg;
       msg << "Do not know how to handle " << it->rah << " for RAH in z-axis orientation bias" << std::endl;
       io::messages.add(msg.str(),
                         "Zaxisoribias_restraints",
                         io::message::error);
       error = 1;
       return error;
    }
    double theta = acos(v[2]/math::abs(v));

    double dist = math::abs(v);
    double vabs = dist;

    (*d_it) = dist;


    DEBUG(9, "ZAXISORIBIAS angle : " << theta << " a0 : " << it->a0);
    double force_scale = 1.0;

    (*ave_it) = 0.0;
    if(rah*theta < rah * it->a0){
      DEBUG(9, "ZAXISORIBIAS  : restraint fulfilled");
      energy = 0;
      f=0;
    } else if (sim.param().zaxisoribias.zaxisoribias == -2) {
        double z2 = v[2] * v[2];
        double vabs2 = vabs * vabs;
        energy = 0.5 * sim.param().zaxisoribias.K * (z2 / vabs2);
        double vabs4 = vabs2 * vabs2;
        f[0] = + sim.param().zaxisoribias.K * z2 / vabs4 * v[0];
        f[1] = + sim.param().zaxisoribias.K * z2 / vabs4 * v[1];
        f[2] = + sim.param().zaxisoribias.K * (z2 - vabs2) / vabs4 * v[2];
    } else if (sim.param().zaxisoribias.zaxisoribias == -1) {
        double zrcos = (v[2] / vabs) - cos(it->a0);
        energy = 0.5 * sim.param().zaxisoribias.K * zrcos * zrcos;
        double vabs2 = vabs * vabs;
        double vabs3 = vabs2 * vabs;
        f[0] = sim.param().zaxisoribias.K * zrcos * v[2] * v[0] / vabs3;
        f[1] = sim.param().zaxisoribias.K * zrcos * v[2] * v[1] / vabs3;
        f[2] = sim.param().zaxisoribias.K * zrcos * (v[2] * v[2]
                                                - vabs2) / vabs3;
    } else if (sim.param().zaxisoribias.zaxisoribias > 0){
        energy = 0.5 * sim.param().zaxisoribias.K *
                        (theta - it->a0) * (theta - it->a0);
        double vabs2 = vabs * vabs;
        double z2 = v[2] * v[2];
        double prefactor = - (sim.param().zaxisoribias.K *
                (theta - it->a0) / sqrt(vabs2 - z2));
        f[0] = prefactor * v[2] * v[0] / vabs2;
        f[1] = prefactor * v[2] * v[1] / vabs2;
        f[2] = prefactor * (z2 - vabs2) / vabs2;
    }


    if(sim.param().zaxisoribias.zaxisoribias == 2){
      energy *= it->w0;
      f *= it->w0;
    }

    f *= force_scale;
    DEBUG(9, "Zaxisoribiasres force : " << math::v2s(f));

    it->v1.force(conf, topo,  f);
    it->v2.force(conf, topo, -f);

    (*ene_it) = energy;

    DEBUG(9, "Zaxisoribiasres energy : " << energy);

    conf.current().energies.zaxisoribias_energy[topo.atom_energy_group()
					    [it->v1.atom(0)]] += energy;

    if (sim.param().zaxisoribias.virial) {
      for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
          conf.current().virial_tensor(a, b) += v(a) * f(b);
        }
      }
    }
  }
  return 0;
}

int interaction::Zaxis_Orientation_Bias_Interaction
::calculate_interactions(topology::Topology &topo,
			 configuration::Configuration &conf,
			 simulation::Simulation &sim)
{
  int error =0;
  SPLIT_VIRIAL_BOUNDARY(_calculate_zaxis_orientation_bias_interactions,
			topo, conf, sim, exponential_term, rah_map, error);

  return error;
}

/**
 * initiate z-axis orientation bias interactions
 */
template<math::boundary_enum B>
static void _init_averages
(topology::Topology & topo,
 configuration::Configuration & conf)
{
  math::Periodicity<B> periodicity(conf.current().box);
  math::Vec v;

  for(std::vector<topology::zaxisori_restraint_struct>::const_iterator
        it = topo.zaxisori_restraints().begin(),
        to = topo.zaxisori_restraints().end(); it != to; ++it) {
    periodicity.nearest_image(it->v1.pos(conf,topo), it->v2.pos(conf,topo), v);
    conf.special().zaxisoribias.d.push_back(math::abs(v));
    conf.special().zaxisoribias.energy.push_back(0.0);
    conf.special().zaxisoribias.av.push_back(pow(math::abs(v), -3.0));
  }
}

/**
 * initiate z-axis orientation bias interactions
 */
template<math::boundary_enum B>
static void _init_zalres_data
(topology::Topology & topo,
 configuration::Configuration & conf)
{
  math::Periodicity<B> periodicity(conf.current().box);
  math::Vec v;

  for(std::vector<topology::zaxisori_restraint_struct>::const_iterator
        it = topo.zaxisori_restraints().begin(),
        to = topo.zaxisori_restraints().end(); it != to; ++it) {
    periodicity.nearest_image(it->v1.pos(conf,topo), it->v2.pos(conf,topo), v);
    conf.special().zaxisoribias.d.push_back(math::abs(v));
    conf.special().zaxisoribias.energy.push_back(0.0);
  }
}

int interaction::Zaxis_Orientation_Bias_Interaction::init(topology::Topology &topo,
		     configuration::Configuration &conf,
		     simulation::Simulation &sim,
		     std::ostream &os,
		     bool quiet)
{

  SPLIT_BOUNDARY(_init_averages, topo, conf);

  if (!quiet) {
    os << "Z-axis orientation bias interaction";
    os << std::endl;
  }
  return 0;
}
