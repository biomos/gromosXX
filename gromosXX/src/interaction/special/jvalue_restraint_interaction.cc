/**
 * @file jvalue_restraint_interaction.cc
 * template methods of Jvalue_Restraint_Interaction
 */

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE special

#include <util/debug.h>

/**
 * calculate position restraint interactions
 */
template<math::boundary_enum b, typename t_interaction_spec>
static int _calculate_jvalue_restraint_interactions
(topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation & sim)
{
  // loop over the jvalue restraints
  std::vector<topology::jvalue_restraint_struct>::const_iterator 
    it = topo.jvalue_restraints().begin(),
    to = topo.jvalue_restraints().end();

  math::VArray &pos   = conf.current().pos;
  math::VArray &force = conf.current().force;
  math::Vec rij, rkj, rkl, rmj, rnk, rim, rln;

  math::Vec v, f;

  double dkj2, dmj2, dnk2, dim, dln;

  double energy;

  math::Periodicity<b> periodicity(conf.current().box);

  for( ; it != to; ++it){

    //get nearest image, calculate rij, rkj, rkl
    periodicity.nearest_image(pos(it->i), pos(it->j), rij);
    periodicity.nearest_image(pos(it->k), pos(it->j), rkj);
    periodicity.nearest_image(pos(it->k), pos(it->l), rkl);
    
    //calculate phi, cross- and dot-products
    rmj = cross(rij, rkj);
    rnk = cross(rkj, rkl);
    dkj2 = dot(rkj, rkj);
    dmj2 = dot(rmj, rmj);
    dnk2 = dot(rnk, rnk);
    
    const double frim = dot(rij, rkj) / dkj2;
    const double frln = dot(rkl, rkj) / dkj2;

    rim = rij - frim * rkj;
    rln = frln * rkj - rkl;
    dim = sqrt(dot(rim, rim));
    dln = sqrt(dot(rln, rln));
    
    const double ip = dot(rim, rln);
    const double cosphi = ip / (dim*dln);
    
    // double phi = acos(dot(rmj, rnk)/(sqrt(dmj2)*sqrt(dnk2)));
    double phi = acos(cosphi);
    
    const double sign = dot(rij, rnk);

    if(sign < 0){ 
      phi = -phi;
    }
    
    DEBUG(10, "JVAL phi: " << phi);
    
    const double cos_phi_delta = cos(phi + it->delta);
    const double sin_phi_delta = sin(phi + it->delta);

    v = 0.0;
    f = 0.0;
    energy = 0.0;

    if (t_interaction_spec::do_virial == math::atomic_virial){
      for(int d0=0; d0<3; ++d0)
	for(int d1=0; d1<3; ++d1)
	  conf.current().virial_tensor(d0, d1) += 
	    v(d0) * f(d1);
      DEBUG(7, "\tatomic virial done");
    }

    conf.current().energies.jvalue_energy[topo.atom_energy_group()
					  [it->i]] += energy;
  }

  return 0;
}

template<typename t_interaction_spec>
int interaction::Jvalue_Restraint_Interaction<t_interaction_spec>
::calculate_interactions(topology::Topology &topo,
			 configuration::Configuration &conf,
			 simulation::Simulation &sim)
{
  switch(conf.boundary_type){
    case math::vacuum :
      return _calculate_jvalue_restraint_interactions<math::vacuum, t_interaction_spec>
	(topo, conf, sim);
      break;
    case math::triclinic :
      return _calculate_jvalue_restraint_interactions<math::triclinic, t_interaction_spec>
	(topo, conf, sim);
      break;
    case math::rectangular :
      return _calculate_jvalue_restraint_interactions<math::rectangular, t_interaction_spec>
	(topo, conf, sim);
      break;
    default:
      throw std::string("wrong boundary type");
  }
  
}
