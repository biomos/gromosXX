/**
 * @file dihedral_interaction.tcc
 * template methods of Dihedral_interaction.
 */

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE interaction

#include "../../debug.h"

/**
 * Constructor.
 */
template<typename t_simulation, typename t_interaction_spec>
inline interaction::Dihedral_interaction<t_simulation, t_interaction_spec>
::Dihedral_interaction()
  : Interaction<t_simulation, t_interaction_spec>("Dihedral")
{
}

/**
 * Destructor.
 */
template<typename t_simulation, typename t_interaction_spec>
inline interaction::Dihedral_interaction<t_simulation, t_interaction_spec>
::~Dihedral_interaction()
{
}

/**
 * calculate dihedral forces and energies.
 */
template<typename t_simulation, typename t_interaction_spec>
inline void interaction::Dihedral_interaction<t_simulation, t_interaction_spec>
::calculate_interactions(t_simulation &sim)
{
  // loop over the improper dihedrals
  std::vector<simulation::Dihedral>::iterator d_it =
    sim.topology().solute().dihedrals().begin(),
    d_to = sim.topology().solute().dihedrals().end();

  math::VArray &pos   = sim.system().pos();
  math::VArray &force = sim.system().force();
  math::Vec rij, rkj, rkl, rlj, rim, rln, rmj, rnk, fi, fj, fk, fl;
  double dkj2, dim, dln, ip;
  double energy;
  
  for( ; d_it != d_to; ++d_it){
    sim.system().periodicity().
      nearest_image(pos(d_it->i), pos(d_it->j), rij);
    sim.system().periodicity().
      nearest_image(pos(d_it->k), pos(d_it->j), rkj);
    sim.system().periodicity().
      nearest_image(pos(d_it->k), pos(d_it->l), rkl);

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

    assert(unsigned(d_it->type) < m_dihedral_parameter.size());
    
    double dcosmphi = 0;
    double cosmphi = 0;
    
    switch(m_dihedral_parameter[d_it->type].m){
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
    double      K = m_dihedral_parameter[d_it->type].K;
    double delta = m_dihedral_parameter[d_it->type].pd;

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

    if (t_interaction_spec::do_virial == atomic_virial){
      sim.system().periodicity().
	nearest_image(pos(d_it->l), pos(d_it->j), rlj);

      for(int a=0; a<3; ++a)
	for(int b=0; b<3; ++b)
	  sim.system().virial()(a, b) += 
	    rij(a) * fi(b) +
	    rkj(a) * fk(b) +
	    rlj(a) * fl(b);

      DEBUG(7, "\tatomic virial done");
    }

    energy = K * (1 + delta * cosmphi);
    sim.system().energies().dihedral_energy
      [sim.topology().atom_energy_group()[d_it->i]] += energy;
    
  }
}

/**
 * add dihedral type.
 */
template<typename t_simulation, typename t_interaction_spec>
inline void interaction::Dihedral_interaction<t_simulation, t_interaction_spec>
::add(dihedral_type_struct s)
{
  m_dihedral_parameter.push_back(s);
}

/**
 * add dihedral type.
 */
template<typename t_simulation, typename t_interaction_spec>
inline void interaction::Dihedral_interaction<t_simulation, t_interaction_spec>
::add(double K, double pd, int m)
{
  dihedral_type_struct s;
  s.K = K;
  s.pd = pd;
  s.m = m;
  add(s);
}

/**
 * access dihedral parameter.
 */
template<typename t_simulation, typename t_interaction_spec>
inline std::vector<interaction::dihedral_type_struct> const &
interaction::Dihedral_interaction<t_simulation, t_interaction_spec>
::parameter()const
{
  return m_dihedral_parameter;
}
