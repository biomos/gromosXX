/**
 * @file contactnum.cc
 * template methods of Contactnum_Colvar
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
#include "../../interaction/special/colvar/contactnum.h"

#include "../../util/template_split.h"
#include "../../util/debug.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE special

/**
 * calculate contactnum restraint interactions
 */
template<math::boundary_enum B, math::virial_enum V>
static int _calculate_contactnum_colvar
(topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation & sim, math::VArray &derivatives, topology::contactnum_restraint_struct *params, double &ct)
{
  //m_timer.start("calculate contacts");
  math::Periodicity<B> periodicity(conf.current().box);
  
  // calculate current number of contacts and derivatives with respect to the
  // position according to formula
  ct=0;
  std::fill(derivatives.begin(), derivatives.end(), math::Vec(0));
  std::cout << "# " << (*params).atoms1.size() << "atoms in group 1,"<< (*params).atoms2.size() << " in 2"<< std::endl;
  for (int i=0; i< (*params).atoms1.size(); i++) {
    for (int j=0; j < (*params).atoms2.size(); j++) {
       math::Vec v;
       double dfunc;
       double func; 
       periodicity.nearest_image((*params).atoms1[i].pos(conf,topo), (*params).atoms2[j].pos(conf,topo), v);
       double rdist=math::abs(v);
      // if(rdist<1.5) {
       func=interaction::switchingfunction(rdist/(*params).rcut,dfunc, (*params).nn, (*params).mm);
       ct+=func;
      // } else {
      // double func=0;
      // double dfunc=0;
      // }
       //if (sim.param().colvarres.write && ((sim.steps() - 1) % sim.param().colvarres.write) == 0) {
       //if (func > 0.1) std::cout << "ct " << i << " " << j <<" gromosnum " << (*params).atoms1[i].atom(0) << " " << (*params).atoms2[j].atom(0) <<" " << func <<" " << ct << " rdist " << rdist << std::endl;
      // }
       math::Vec d = dfunc*v/rdist;
       derivatives[i]+=d;
       derivatives[(*params).atoms1.size()+j]-=d;
    }
  }  
  //m_timer.stop("calculate contacts");
  return 0;
}

int interaction::Contactnum_Colvar
::calculate(topology::Topology &topo,
			 configuration::Configuration &conf,
			 simulation::Simulation &sim)
{
  SPLIT_VIRIAL_BOUNDARY(_calculate_contactnum_colvar,
			topo, conf, sim, derivatives,params,ct);
  
  return 0;
}



/**
 * initiate contactnum restraint interactions
 */
int interaction::Contactnum_Colvar::init(topology::Topology &topo, 
		     configuration::Configuration &conf,
		     simulation::Simulation &sim,
		     std::ostream &os,
		     bool quiet) 
{

  targetvalue=(*params).cont0;
  rcut=(*params).rcut;
  mm=(*params).mm;
  nn=(*params).nn;
  w0=(*params).w0;
  
  // atoms will be one concatenated list of pointers to atoms1 and atoms2 from
  // the restraint specification that can be given to the colvar restraint interaction
  for (std::vector<util::Virtual_Atom >::iterator it = (*params).atoms1.begin(),
        to = (*params).atoms1.end(); it != to; ++it) {
    atoms.push_back(&(*it));
  }
  for (std::vector<util::Virtual_Atom >::iterator it = (*params).atoms2.begin(),
        to = (*params).atoms2.end(); it != to; ++it) {
    atoms.push_back(&(*it));
  }
  
  derivatives.resize((*params).atoms1.size()+(*params).atoms2.size());
  

  
  if (!quiet) {
    os << "Contactnum restraint interaction";
    os << std::endl;
  }
  return 0;
}
