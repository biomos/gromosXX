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
  math::Vec rij, rkj, rkl, rmj, rnk, rim, rln;

  // math::Vec v, f;

  double dkj2, dmj2, dnk2, dim, dln;

  math::Periodicity<b> periodicity(conf.current().box);

  int n = 0;
  for( ; it != to; ++it, ++n){

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

    ////////////////////////////////////////////////////////////////////////////////

    double exp_term;
    double memory_decay;

    //decide on time averaging
    if (sim.param().jvalue.mode != simulation::restr_inst){
 
      // time averaging 
      exp_term = exp(-sim.time_step_size() / sim.param().jvalue.tau);
      memory_decay = 1 - exp(-sim.time_step_size() / sim.param().jvalue.tau);
    }
    else{
      
      // instantaneous J-values
      exp_term = 0;
      memory_decay = 1;
    }

    DEBUG(10, "JVAL exp_term " << exp_term);
    DEBUG(10, "JVAL memory_decay " << memory_decay);
        
    // calculate J-value
    assert(conf.special().jvalue_av.size() > unsigned(n));
    
    const double Jav =   it->a * memory_decay * cos_phi_delta * cos_phi_delta 
      + it->b * memory_decay * cos_phi_delta    
      + it->c * memory_decay
      + conf.special().jvalue_av[n] * exp_term;
 
    const double Jcurr = it->a * cos_phi_delta * cos_phi_delta
      + it->b * cos_phi_delta
      + it->c;

    // save current J
    // d_it->J = Jcurr;
    
    DEBUG(10, "JDATA time: " << sim.time() << "   Jcurr: " << Jcurr
	  << "\tJav: " << Jav);
    
    //write new average and current value
    conf.special().jvalue_av[n] = Jav;
    conf.special().jvalue_curr[n] = Jcurr;
    
    const double dV_dphi = 
      _calculate_derivative(topo, conf, sim.param(),
			    it,
			    Jcurr, Jav,
			    cos_phi_delta, sin_phi_delta);
    
    //calculate forces 		 
    const math::Vec dphi_dri =  (sqrt(dkj2)/dmj2)*rmj;
    const math::Vec dphi_drl = -(sqrt(dkj2)/dnk2)*rnk;			
    const math::Vec dphi_drj = (frim -1)*dphi_dri - frln*dphi_drl;
    const math::Vec dphi_drk = -1.0*dphi_dri - dphi_drj - dphi_drl;			
 
    const math::Vec fi = - dV_dphi * dphi_dri;
    const math::Vec fj = - dV_dphi * dphi_drj;
    const math::Vec fk = - dV_dphi * dphi_drk; 
    const math::Vec fl = - dV_dphi * dphi_drl;
    
    DEBUG(10, "JVAL Force on i " << math::v2s(fi) );
    DEBUG(10, "JVAL Force on j " << math::v2s(fj) );
    DEBUG(10, "JVAL Force on k " << math::v2s(fk) );
    DEBUG(10, "JVAL Force on l " << math::v2s(fl) );
 
  
    conf.current().force(it->i) += fi;
    conf.current().force(it->j) += fj;
    conf.current().force(it->k) += fk;
    conf.current().force(it->l) += fl;

    if (t_interaction_spec::do_virial == math::atomic_virial){
      math::Vec rlj;
      periodicity.nearest_image(pos(it->l), pos(it->j), rlj);

      for(int d0=0; d0<3; ++d0)
	for(int d1=0; d1<3; ++d1)
	  conf.current().virial_tensor(d0, d1) += 
	    rij(d0) * fi(d1) +
	    rkj(d0) * fk(d1) +
	    rlj(d0) * fl(d1);
      
      DEBUG(7, "\tatomic virial done");
    }

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

static double _calculate_derivative(topology::Topology & topo,
				    configuration::Configuration &conf,
				    simulation::Parameter const & param,
				    std::vector<topology::jvalue_restraint_struct>::const_iterator it,
				    double Jcurr, double Jav,
				    double cos_phi_delta, double sin_phi_delta)
{
  // instantaneous / time averaged
  if (param.jvalue.mode == simulation::restr_inst ||
      param.jvalue.mode == simulation::restr_av){

    // check for half - harmoic functional forms
    if ( (it->H == topology::repulsive && Jav - it->J0 > 0) || 
	 (it->H == topology::attractive && Jav - it->J0 <= 0) ){
      return 0;
    }
    else{
      // calculate derivatives + energy	
      // Jav == Jcurr for instantaneous...
      const double dV_dJ = it->K * (Jav - it->J0);
      
      //memory_decay factor is omitted for practical reasons.
      const double dJ_dphi = - (2 * it->a * cos_phi_delta * sin_phi_delta + it->b * sin_phi_delta);
      
      const double energy = 0.5 * it->K * (Jav - it->J0) * (Jav - it->J0);
      conf.current().energies.jvalue_energy[topo.atom_energy_group()
					   [it->i]]
	+= energy;

      return dV_dJ * dJ_dphi;
    }
  }
  else if (param.jvalue.mode == simulation::restr_biq){

    // check for half - harmoic functional forms
    if ( (it->H == topology::repulsive && (Jcurr - it->J0 > 0 || Jav - it->J0 > 0)) || 
	 (it->H == topology::attractive && (Jcurr - it->J0 <= 0 || Jav - it->J0 <= 0)) ){
      return 0;
    }
    else{
      // calculate derivatives + energy	

      const double dV_dJ = it->K * (Jcurr - it->J0) * (Jav - it->J0) * (Jav - it->J0);
      //memory_decay factor is omitted for practical reasons.
      const double dJ_dphi = - (2 * it->a * cos_phi_delta * sin_phi_delta + it->b * sin_phi_delta);
      

      const double dV_dJav = it->K * (Jcurr - it->J0) * (Jcurr - it->J0) * (Jav - it->J0);
      //memory_decay factor is omitted for practical reasons.
      const double dJav_dphi = - (2 * it->a * cos_phi_delta * sin_phi_delta + it->b * sin_phi_delta);
      
      const double energy = 0.5 * it->K
	* (Jcurr - it->J0) * (Jcurr - it->J0) * (Jav - it->J0) * (Jav - it->J0);

      conf.current().energies.jvalue_energy[topo.atom_energy_group()[it->i]]
	+= energy;

      return dV_dJ * dJ_dphi + dV_dJav * dJav_dphi;
    }
  }

  throw std::string("JValue restraints: derivative not implemented");
  return 0;
  
}