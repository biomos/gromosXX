/**
 * @file distance.cc
 * template methods of Distance_Colvar
 */

#include <limits>
#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"
#include "../../interaction/interaction.h"

#include "../../math/periodicity.h"

// special interactions
#include "../../interaction/interaction_types.h"

#include "../../interaction/special/colvar/colvar.h"
#include "../../interaction/special/colvar/distance.h"

#include "../../util/template_split.h"
#include "../../util/debug.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE special


/**
 * calculate distance restraint interactions
 */
template<math::boundary_enum B, math::virial_enum V>
static int _calculate_distance_colvar
(topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation & sim,
 math::VArray &derivatives, 
 topology::distance_restraint_struct_colvar *params,
 double &ct)
 {
  math::Periodicity<B> periodicity(conf.current().box);
  ct=0;
  math::Vec v;
  std::fill(derivatives.begin(), derivatives.end(), math::Vec(0));

  // Calculate vector between atoms or centers of mass:
  math::Vec pos1 = (*params).atoms1.pos(conf, topo);
  math::Vec pos2 = (*params).atoms2.pos(conf, topo);

  periodicity.nearest_image(pos1, pos2, v);
  DEBUG(8, "v1 atom: " <<  (*params).atoms1.atom(0)+1);
  DEBUG(8, "v2 atom: " <<  (*params).atoms2.atom(0)+1);
  DEBUG(8, "pos(v1) = " <<  math::v2s((*params).atoms1.pos(conf,topo)));
  DEBUG(8, "pos(v2) = " <<  math::v2s((*params).atoms2.pos(conf,topo)));
  double dist = math::abs(v);
  ct = dist;

  // derivative: gradient of distance wrt atom positions
  math::Vec d = v / dist;

  // For single atoms:
  derivatives[0] = d;
  derivatives[1] = -d;

  return 0;
}

int interaction::Distance_Colvar
::calculate(topology::Topology &topo,
			 configuration::Configuration &conf,
			 simulation::Simulation &sim)
{
  SPLIT_VIRIAL_BOUNDARY(_calculate_distance_colvar,
			topo, conf, sim, derivatives, params, ct);
  
  return 0;
}

/**
 * initiate distance restraint interactions
 */
int interaction::Distance_Colvar::init(topology::Topology &topo, 
		     configuration::Configuration &conf,
		     simulation::Simulation &sim,
		     std::ostream &os,
		     bool quiet) 
{

  targetvalue=(*params).r0;
  rah=(*params).rah;
  d0=(*params).d0;
  w0=(*params).w0;
  
  // atoms will be one concatenated list of pointers to atoms1 and atoms2 from
  // the restraint specification that can be given to the colvar restraint interaction
  atoms.push_back(&params->atoms1);
  atoms.push_back(&params->atoms2);
  
  derivatives.resize(2);
  
  if (!quiet) {
    os << "Distance restraint interaction";
    os << std::endl;
  }
  return 0;
}
