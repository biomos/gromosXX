/**
 * @file dihedral_interaction.tcc
 * template methods of Dihedral_interaction.
 */

/**
 * Destructor.
 */
template<typename t_simulation>
inline interaction::Dihedral_interaction<t_simulation>
::~Dihedral_interaction()
{
}

/**
 * calculate dihedral forces and energies.
 */
template<typename t_simulation>
inline void interaction::Dihedral_interaction<t_simulation>
::calculate_interactions(t_simulation &simu)
{
  // loop over the improper dihedrals
  simulation::Dihedral::iterator d_it =
    simu.topology().solute().dihedrals().begin();

  math::VArray &pos   = simu.system().pos();
  math::VArray &force = simu.system().force();
  math::Vec rij, rkj, rkl, rim, rln, rmj, rnk, fi, fj, fk, fl;
  double dkj2, dim, dln, ip;
  
  for( ; !d_it.eol(); ++d_it){
    rkj = pos(d_it.k()) - pos(d_it.l());
    rij = pos(d_it.i()) - pos(d_it.j());
    rkl = pos(d_it.k()) - pos(d_it.l());
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

    assert(unsigned(d_it.type()) < m_dihedral_parameter.size());
    
    double dcosmphi = 0;
    switch(m_dihedral_parameter[d_it.type()].m){
      case 0:
	dcosmphi = 0.0;
	break;
      case 1:
	dcosmphi = 1;
	break;
      case 2:
	dcosmphi = 4*cosphi;
	break;
      case 3:
	dcosmphi = 12*cosphi2 - 3;
	break;
      case 4:
	dcosmphi = 32*cosphi3-16*cosphi;
	break;
      case 5:
	dcosmphi = 80*cosphi4-60*cosphi2+5;
	break;
      case 6:
	dcosmphi = 192*cosphi4*cosphi-192*cosphi3+36*cosphi;
	break;
      default:
	//io::messages.add("dihedral function not implemented for m>6", 
	//		 "dihedral_interaction", io::message::error);
	;
	
    }
    double      K = m_dihedral_parameter[d_it.type()].K;
    double delta = m_dihedral_parameter[d_it.type()].pd;


    double ki = -K * delta * dcosmphi / dim;
    double kl = -K * delta * dcosmphi / dln;
    double kj1 = frim - 1.0;
    double kj2 = frln;
    
    fi = ki * (rln / dln - rim / dim * cosphi);
    fl = kl * (rim / dim - rln / dln * cosphi);
    fj = kj1 * fi - kj2 * fl;
    fk = -1.0 * (fi + fj + fl);
    
    force(d_it.i()) += fi;
    force(d_it.j()) += fj;
    force(d_it.k()) += fk;
    force(d_it.l()) += fl;
  }
}

/**
 * add dihedral type.
 */
template<typename t_simulation>
inline void interaction::Dihedral_interaction<t_simulation>
::add(dihedral_type_struct s)
{
  m_dihedral_parameter.push_back(s);
}

/**
 * add dihedral type.
 */
template<typename t_simulation>
inline void interaction::Dihedral_interaction<t_simulation>
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
template<typename t_simulation>
inline std::vector<interaction::dihedral_type_struct> const &
interaction::Dihedral_interaction<t_simulation>
::parameter()const
{
  return m_dihedral_parameter;
}
