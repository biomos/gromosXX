/**
 * @file dihedral_interaction.tcc
 * template methods of Dihedral_interaction.
 */

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE interaction

#include <util/debug.h>

/**
 * calculate dihedral forces and energies.
 */
template<math::boundary_enum b, typename t_interaction_spec>
static int _calculate_dihedral_interactions(topology::Topology & topo,
					    configuration::Configuration & conf,
					    simulation::Simulation & sim,
					    std::vector<interaction::dihedral_type_struct> 
					    const & param)
{
  // loop over the improper dihedrals
  std::vector<topology::four_body_term_struct>::iterator d_it =
    topo.solute().dihedrals().begin(),
    d_to = topo.solute().dihedrals().end();

  math::VArray &pos   = conf.current().pos;
  math::VArray &force = conf.current().force;
  math::Vec rij, rkj, rkl, rlj, rim, rln, rmj, rnk, fi, fj, fk, fl;
  double dkj2, dim, dln, ip;
  double energy;
  
  math::Periodicity<b> periodicity(conf.current().box);

  for(int n =0; d_it != d_to; ++d_it, ++n){
    periodicity.nearest_image(pos(d_it->i), pos(d_it->j), rij);
    periodicity.nearest_image(pos(d_it->k), pos(d_it->j), rkj);
    periodicity.nearest_image(pos(d_it->k), pos(d_it->l), rkl);

    rmj = cross(rij, rkj);
    rnk = cross(rkj, rkl);
    
    dkj2 = dot(rkj, rkj);

    double frim = dot(rij, rkj)/dkj2;
    double frln = dot(rkl, rkj)/dkj2;

    rim = rij - frim * rkj;
    rln = frln * rkj - rkl;
    dim = sqrt(dot(rim, rim));
    dln = sqrt(dot(rln, rln));
    
    ip = dot(rim, rln);
    double cosphi = ip / (dim*dln);
    
    double cosphi2 = cosphi  * cosphi;
    double cosphi3 = cosphi2 * cosphi;
    double cosphi4 = cosphi3 * cosphi;

    assert(unsigned(d_it->type) < param.size());
    
    double dcosmphi = 0;
    double cosmphi = 0;
    
    switch(param[d_it->type].m){
      case 0:
	cosmphi = 0.0;
	dcosmphi = 0.0;
	break;
      case 1:
	cosmphi = cosphi;
	dcosmphi = 1;
	break;
      case 2:
	cosmphi =  2*cosphi2 -1;
	dcosmphi = 4*cosphi;
	break;
      case 3:
	cosmphi  = 4*cosphi3 - 3*cosphi;
	dcosmphi = 12*cosphi2 - 3;
	break;
      case 4:
	cosmphi  = 8*cosphi4 - 8*cosphi2 + 1;
	dcosmphi = 32*cosphi3-16*cosphi;
	break;
      case 5:
	cosmphi  = 16*cosphi4*cosphi - 20*cosphi3 + 5*cosphi;
	dcosmphi = 80*cosphi4-60*cosphi2+5;
	break;
      case 6:
	cosmphi  = 32*cosphi4*cosphi2 - 48*cosphi4 + 18*cosphi2 -1;
	dcosmphi = 192*cosphi4*cosphi-192*cosphi3+36*cosphi;
	break;
      default:
	//io::messages.add("dihedral function not implemented for m>6", 
	//		 "dihedral_interaction", io::message::error);
	throw std::runtime_error("dihedral type for m=6 not implemented");
	
    }
    double     K = param[d_it->type].K;
    double delta = param[d_it->type].pd;

    DEBUG(10, "dihedral K=" << K << " delta=" << delta << " dcos=" << dcosmphi);

    double ki = -K * delta * dcosmphi / dim;
    double kl = -K * delta * dcosmphi / dln;
    double kj1 = frim - 1.0;
    double kj2 = frln;
    
    fi = ki * (rln / dln - rim / dim * cosphi);
    fl = kl * (rim / dim - rln / dln * cosphi);
    fj = kj1 * fi - kj2 * fl;
    fk = -1.0 * (fi + fj + fl);
    
    force(d_it->i) += fi;
    force(d_it->j) += fj;
    force(d_it->k) += fk;
    force(d_it->l) += fl;

    if (t_interaction_spec::do_virial == math::atomic_virial){
      periodicity.nearest_image(pos(d_it->l), pos(d_it->j), rlj);

      for(int a=0; a<3; ++a)
	for(int bb=0; bb<3; ++bb)
	  conf.current().virial_tensor(a, bb) += 
	    rij(a) * fi(bb) +
	    rkj(a) * fk(bb) +
	    rlj(a) * fl(bb);

      DEBUG(11, "\tatomic virial done");
    }

    energy = K * (1 + delta * cosmphi);
    conf.current().energies.dihedral_energy
      [topo.atom_energy_group()[d_it->i]] += energy;
    
    // dihedral angle monitoring.
    if(sim.param().print.monitor_dihedrals){
      DEBUG(8, "monitoring dihedrals");
      
      double phi = acos(cosphi);
      ip = dot(rij, rnk);
      if(ip < 0) phi*=-1.0;
      DEBUG(11, "dihedral angle: " << phi
	    << " previous minimum: " <<conf.special().dihedral_angle_minimum[n]);
      
      if(fabs(conf.special().dihedral_angle_minimum[n] - phi) > 
	 2*math::Pi / param[d_it->type].m){
	double old_min=conf.special().dihedral_angle_minimum[n];
	conf.special().dihedral_angle_minimum[n] = 
	  _calculate_nearest_minimum(phi, param[d_it->type].m, delta);
	// ugly check to see that it is not the first...
	if(old_min != 4*math::Pi){
	  // could be written to a separate file or by a separate function
	  // should at least be more descriptive.
	  std::cout << "transition from: " << old_min << " to " 
		    << conf.special().dihedral_angle_minimum[n] << "\n";
	}
      }
    }
  }
  
  return 0;
  
}

template<typename t_interaction_spec>
int interaction::Dihedral_Interaction<t_interaction_spec>
::calculate_interactions(topology::Topology &topo,
			 configuration::Configuration &conf,
			 simulation::Simulation &sim)
{
  switch(conf.boundary_type){
    case math::vacuum :
      return _calculate_dihedral_interactions<math::vacuum, t_interaction_spec>
	(topo, conf, sim, m_parameter);
      break;
    case math::triclinic :
      return _calculate_dihedral_interactions<math::triclinic, t_interaction_spec>
	(topo, conf, sim, m_parameter);
      break;
    case math::rectangular :
      return _calculate_dihedral_interactions<math::rectangular, t_interaction_spec>
	(topo, conf, sim, m_parameter);
      break;
    default:
      throw std::string("Wrong boundary type");
  }
  
}

/**
 * calculate nearest minimum
 */
double _calculate_nearest_minimum(double phi, int m, double pd)
{
  // copy from gromos++ nearest_minimum function
  double a_minimum = 0.5*math::Pi*(3.0 - pd)/ m;
  double delta_phi = 2*math::Pi / m;
  double nearest_min = a_minimum - int(rint((a_minimum - phi)/delta_phi))*delta_phi;
  if(nearest_min >= 2*math::Pi - math::epsilon) nearest_min -= 2*math::Pi;
  
  return nearest_min;
}

