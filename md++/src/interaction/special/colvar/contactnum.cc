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

namespace {
  double fastpow(double base, int exp) {
    if (exp < 0) {
      exp = -exp;
      base = 1.0 / base;
    }
    double result = 1.0;
    while (exp) {
      if (exp & 1)
        result *= base;
      exp >>= 1;
      base *= base;
    }
    return result;
  }

  double switchingfunction(double rdist, double& dfunc, int nn, int mm) {
    const double epsilon(std::numeric_limits<double>::epsilon());
    double result;

    if (2 * nn == mm) {
      double rNdist = fastpow(rdist, nn - 1);
      double iden = 1.0 / (1.0 + rNdist * rdist);
      dfunc = -nn * rNdist * iden * iden;
      result = iden;
    } else {
      if (rdist > (1.0 - 100.0 * epsilon) && rdist < (1.0 + 100.0 * epsilon)) {
        result = static_cast<double>(nn) / mm;
        dfunc = 0.5 * nn * (nn - mm) / mm;
      } else {
        double rNdist = fastpow(rdist, nn - 1);
        double rMdist = fastpow(rdist, mm - 1);
        double num = 1.0 - rNdist * rdist;
        double iden = 1.0 / (1.0 - rMdist * rdist);
        double func = num * iden;
        result = func;
        dfunc = ((-nn * rNdist * iden) + (func * (iden * mm) * rMdist));
      }
    }
    return result;
  }
} // unnamed namespace


/**
 * calculate contactnum restraint interactions
 */
template<math::boundary_enum B, math::virial_enum V>
static int _calculate_contactnum_colvar
(topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation & sim,
 math::VArray &derivatives,
 topology::contactnum_restraint_struct *params,
 double &ct)
{
  math::Periodicity<B> periodicity(conf.current().box);
  
  // calculate current number of contacts and derivatives with respect to the
  // position according to formula
  ct=0;
  std::fill(derivatives.begin(), derivatives.end(), math::Vec(0));
  for (int i=0; i< (*params).atoms1.size(); i++) {
    for (int j=0; j < (*params).atoms2.size(); j++) {
       math::Vec v;
       double dfunc;
       double func; 
       periodicity.nearest_image((*params).atoms1[i].pos(conf,topo), (*params).atoms2[j].pos(conf,topo), v);
       double rdist=math::abs(v);
       //func=interaction::switchingfunction(rdist/(*params).rcut, dfunc, (*params).nn, (*params).mm);
       func = switchingfunction(rdist / params->rcut, dfunc, params->nn, params->mm);

       ct+=func;
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
