/**
 * @file nonbonded_innerloop.cc
 * template methods of Nonbonded_Innerloop
 */

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE nonbonded

template<typename t_nonbonded_spec>
inline void 
interaction::Nonbonded_Innerloop::lj_crf_innerloop
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 unsigned int i,
 unsigned int j,
 Storage & storage,
 math::Periodicity<t_nonbonded_spec::boundary_type> const & periodicity
 )
{
  DEBUG(8, "\tpair\t" << i << "\t" << j);
  
  math::Vec r;
  double f;
  double e_lj, e_crf;
  
  periodicity.nearest_image(conf.current().pos(i), 
			    conf.current().pos(j), r);
  DEBUG(10, "\tni r " << r(0) << " / " << r(1) << " / " << r(2));
  
  const lj_parameter_struct &lj = 
    m_param->lj_parameter(topo.iac(i),
			  topo.iac(j));
  
  DEBUG(11, "\tlj-parameter c6=" << lj.c6 << " c12=" << lj.c12);
  DEBUG(11, "\tcharge i=" << topo.charge()(i) << " j=" << topo.charge()(j));
  
  switch(t_nonbonded_spec::interaction_func){
    case simulation::lj_crf_func :
      lj_crf_interaction(r, lj.c6, lj.c12,
			 topo.charge()(i) * 
			 topo.charge()(j),
			 f, e_lj, e_crf);
      break;
    default:
      io::messages.add("Nonbonded_Innerloop",
		       "interaction function not implemented",
		       io::message::critical);
  }
  
  // most common case
  if (t_nonbonded_spec::do_virial == math::molecular_virial){
    math::Vec rf = f * r;
    storage.force(i) += rf;
    storage.force(j) -= rf;
    
    for(int b=0; b<3; ++b){
      const double rr = r(b) - conf.special().rel_mol_com_pos(i)(b) 
	+ conf.special().rel_mol_com_pos(j)(b);
      for(int a=0; a<3; ++a){
	storage.virial_tensor(b, a) += rr * rf(a);
      }
    }
  }
  else{
    for (int a=0; a<3; ++a){
      
      const double term = f * r(a);
      storage.force(i)(a) += term;
      storage.force(j)(a) -= term;
      
      if (t_nonbonded_spec::do_virial == math::atomic_virial){
	for(int b=0; b<3; ++b){
	  storage.virial_tensor(b, a) += 
	    r(b) * term;
	}
      }
    }
  }
  
  // energy
  DEBUG(11, "\tenergy group i " << topo.atom_energy_group(i)
	<< " j " << topo.atom_energy_group(j));

  assert(storage.energies.lj_energy.size() > 
	 topo.atom_energy_group(i));
  assert(storage.energies.lj_energy.size() >
	 topo.atom_energy_group(j));
  
  storage.energies.lj_energy[topo.atom_energy_group(i)]
    [topo.atom_energy_group(j)] += e_lj;
  
  storage.energies.crf_energy[topo.atom_energy_group(i)]
    [topo.atom_energy_group(j)] += e_crf;
  
}

inline void 
interaction::Nonbonded_Innerloop::lj_crf_innerloop_central
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 unsigned int i,
 unsigned int j,
 Storage & storage
 )
{
  DEBUG(8, "\tpair\t" << i << "\t" << j);
  
  double f;
  double e_lj, e_crf;
  
  const math::Vec r = conf.current().pos(i) - conf.current().pos(j);
  
  const lj_parameter_struct &lj = 
    m_param->lj_parameter(topo.iac(i),
			  topo.iac(j));
  
  DEBUG(11, "\tlj-parameter c6=" << lj.c6 << " c12=" << lj.c12);
  DEBUG(11, "\tcharge i=" << topo.charge()(i) << " j=" << topo.charge()(j));
  
  lj_crf_interaction(r, lj.c6, lj.c12,
		     topo.charge()(i) * 
		     topo.charge()(j),
		     f, e_lj, e_crf);
  
  // most common case
  // if (t_nonbonded_spec::do_virial == math::molecular_virial){
  math::Vec rf = f * r;
  storage.force(i) += rf;
  storage.force(j) -= rf;
  
  for(int b=0; b<3; ++b){
    const double rr = r(b) - conf.special().rel_mol_com_pos(i)(b) + conf.special().rel_mol_com_pos(j)(b);
    for(int a=0; a<3; ++a){
      storage.virial_tensor(b, a) += rr * rf(a);
    }
  }

  /*
  }
  else{
    for (int a=0; a<3; ++a){
      
      const double term = f * r(a);
      storage.force(i)(a) += term;
      storage.force(j)(a) -= term;
      
      if (t_nonbonded_spec::do_virial == math::atomic_virial){
	for(int b=0; b<3; ++b){
	  storage.virial_tensor(b, a) += 
	    r(b) * term;
	}
      }
    }
  }
  */
  
  // energy
  DEBUG(11, "\tenergy group i " << topo.atom_energy_group(i)
	<< " j " << topo.atom_energy_group(j));

  assert(storage.energies.lj_energy.size() > 
	 topo.atom_energy_group(i));
  assert(storage.energies.lj_energy.size() >
	 topo.atom_energy_group(j));
  
  storage.energies.lj_energy[topo.atom_energy_group(i)]
    [topo.atom_energy_group(j)] += e_lj;
  
  storage.energies.crf_energy[topo.atom_energy_group(i)]
    [topo.atom_energy_group(j)] += e_crf;
  
}

inline void 
interaction::Nonbonded_Innerloop::lj_crf_innerloop_shift
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 unsigned int i,
 unsigned int j,
 Storage & storage,
 math::Vec const & shift
 )
{
  DEBUG(8, "\tpair\t" << i << "\t" << j);
  
  double f;
  double e_lj, e_crf;
  
  math::Vec r = conf.current().pos(i) - conf.current().pos(j);
  r -= shift;
  
  const lj_parameter_struct &lj = 
    m_param->lj_parameter(topo.iac(i),
			  topo.iac(j));
  
  DEBUG(11, "\tlj-parameter c6=" << lj.c6 << " c12=" << lj.c12);
  DEBUG(11, "\tcharge i=" << topo.charge()(i) << " j=" << topo.charge()(j));
  
  lj_crf_interaction(r, lj.c6, lj.c12,
		     topo.charge()(i) * 
		     topo.charge()(j),
		     f, e_lj, e_crf);
  
  // most common case
  // if (t_nonbonded_spec::do_virial == math::molecular_virial){
  math::Vec rf = f * r;
  storage.force(i) += rf;
  storage.force(j) -= rf;
  
  for(int b=0; b<3; ++b){
    const double rr = r(b) - conf.special().rel_mol_com_pos(i)(b) + conf.special().rel_mol_com_pos(j)(b);
    for(int a=0; a<3; ++a){
      storage.virial_tensor(b, a) += rr * rf(a);
    }
  }

  /*
  }
  else{
    for (int a=0; a<3; ++a){
      
      const double term = f * r(a);
      storage.force(i)(a) += term;
      storage.force(j)(a) -= term;
      
      if (t_nonbonded_spec::do_virial == math::atomic_virial){
	for(int b=0; b<3; ++b){
	  storage.virial_tensor(b, a) += 
	    r(b) * term;
	}
      }
    }
  }
  */
  
  // energy
  DEBUG(11, "\tenergy group i " << topo.atom_energy_group(i)
	<< " j " << topo.atom_energy_group(j));

  assert(storage.energies.lj_energy.size() > 
	 topo.atom_energy_group(i));
  assert(storage.energies.lj_energy.size() >
	 topo.atom_energy_group(j));
  
  storage.energies.lj_energy[topo.atom_energy_group(i)]
    [topo.atom_energy_group(j)] += e_lj;
  
  storage.energies.crf_energy[topo.atom_energy_group(i)]
    [topo.atom_energy_group(j)] += e_crf;
  
}


template<typename t_nonbonded_spec>
void interaction::Nonbonded_Innerloop::one_four_interaction_innerloop
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 int i,
 int j,
 math::Periodicity<t_nonbonded_spec::boundary_type> const & periodicity)
{
  DEBUG(8, "\t1,4-pair\t" << i << "\t" << j);
  
  math::Vec r;
  double f, e_lj, e_crf;
  
  periodicity.nearest_image(conf.current().pos(i), 
			    conf.current().pos(j), r);

  const lj_parameter_struct &lj = 
    m_param->lj_parameter(topo.iac(i),
			  topo.iac(j));
  
  DEBUG(11, "\tlj-parameter cs6=" << lj.cs6 << " cs12=" << lj.cs12);

  switch(t_nonbonded_spec::interaction_func){
    case simulation::lj_crf_func :
      lj_crf_interaction(r, lj.cs6, lj.cs12,
			 topo.charge()(i) * 
			 topo.charge()(j),
			 f, e_lj, e_crf);
      break;
    default:
      io::messages.add("Nonbonded_Innerloop",
		       "interaction function not implemented",
		       io::message::critical);
  }

  for (int a=0; a<3; ++a){

    const double term = f * r(a);
    // storage.force(i)(a) += term;
    // storage.force(j)(a) -= term;
    conf.current().force(i)(a) += term;
    conf.current().force(j)(a) -= term;
    
    if (t_nonbonded_spec::do_virial == math::atomic_virial){
      for(int b=0; b<3; ++b){
	conf.current().virial_tensor(b, a) += 
	  r(b) * term;
      }
    }
  }
  
  // energy
  conf.current().energies.lj_energy[topo.atom_energy_group(i)]
    [topo.atom_energy_group(j)] += e_lj;
    
  conf.current().energies.crf_energy[topo.atom_energy_group(i)]
    [topo.atom_energy_group(j)] += e_crf;
    
  DEBUG(11, "\tenergy group i " << topo.atom_energy_group(i)
	<< " j " << topo.atom_energy_group(j));
    
}

template<typename t_nonbonded_spec>
inline void 
interaction::Nonbonded_Innerloop::RF_excluded_interaction_innerloop
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 int i,
 math::Periodicity<t_nonbonded_spec::boundary_type> const & periodicity)
{
  math::Vec r, f;
  double e_crf;
  
  math::VArray &pos   = conf.current().pos;
  math::VArray &force = conf.current().force;

  std::set<int>::const_iterator it, to;
  it = topo.exclusion(i).begin();
  to = topo.exclusion(i).end();
  
  DEBUG(8, "\tself-term " << i );
  r=0;
  
  // this will only contribute in the energy, the force should be zero.
  rf_interaction(r,topo.charge()(i) * topo.charge()(i),
		 f, e_crf);
  conf.current().energies.crf_energy[topo.atom_energy_group(i)]
    [topo.atom_energy_group(i)] += 0.5 * e_crf;
  DEBUG(11, "\tcontribution " << 0.5*e_crf);
  
  for( ; it != to; ++it){
    
    DEBUG(11, "\texcluded pair " << i << " - " << *it);
    
    periodicity.nearest_image(pos(i), pos(*it), r);
    
    
    rf_interaction(r, topo.charge()(i) * 
		   topo.charge()(*it),
		   f, e_crf);
    
    force(i) += f;
    force(*it) -= f;
    
    if (t_nonbonded_spec::do_virial == math::atomic_virial){
      for(int a=0; a<3; ++a)
	for(int b=0; b<3; ++b)
	  conf.current().virial_tensor(a, b) += 
	    r(a) * f(b);
      DEBUG(11, "\tatomic virial done");
    }

    // energy
    conf.current().energies.crf_energy[topo.atom_energy_group(i)]
      [topo.atom_energy_group(*it)] += e_crf;
    
    DEBUG(11, "\tcontribution " << e_crf);
    
  } // loop over excluded pairs
  
}

template<typename t_nonbonded_spec>
inline void 
interaction::Nonbonded_Innerloop::RF_solvent_interaction_innerloop
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 topology::Chargegroup_Iterator const & cg_it,
 math::Periodicity<t_nonbonded_spec::boundary_type> const & periodicity
 )
{
  math::Vec r;
  double e_crf;
  
  math::VArray &pos   = conf.current().pos;

  // loop over the atoms
  topology::Atom_Iterator at_it = cg_it.begin(),
    at_to = cg_it.end();
  
  for ( ; at_it != at_to; ++at_it){
    DEBUG(11, "\tsolvent self term " << *at_it);
    // no solvent self term. The distance dependent part and the forces
    // are zero. The distance independent part should add up to zero 
    // for the energies and is left out.
    
    for(topology::Atom_Iterator at2_it=at_it+1; at2_it!=at_to; ++at2_it){
      
      DEBUG(11, "\tsolvent " << *at_it << " - " << *at2_it);
      periodicity.nearest_image(pos(*at_it), 
				pos(*at2_it), r);
      
      // for solvent, we don't calculate internal forces (rigid molecules)
      // and the distance independent parts should go to zero
      e_crf = -topo.charge()(*at_it) * 
	topo.charge()(*at2_it) * 
	math::four_pi_eps_i * 
	crf_2cut3i() * abs2(r);
      
      // energy
      conf.current().energies.crf_energy
	[topo.atom_energy_group(*at_it) ]
	[topo.atom_energy_group(*at2_it)] += e_crf;
    } // loop over at2_it
  } // loop over at_it
  
}

template<typename t_nonbonded_spec>
inline void 
interaction::Nonbonded_Innerloop::spc_innerloop
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 int i,
 int j,
 Storage & storage,
 math::Periodicity<t_nonbonded_spec::boundary_type> const & periodicity
 )
{
  math::Vec r;
  
  double x[4], y[4], z[4], r2[4], r2i[4], ri[4], ff[4], tx, ty, tz, fx, fy, fz, rx, ry, rz;

  // only one energy group
  const int egroup = topo.atom_energy_group(topo.num_solute_atoms());
  DEBUG(8, "\tspc pair\t" << i << "\t" << j << " egroup " << egroup);
  
  const int ii = topo.num_solute_atoms() + i * 3;
  const int jj = topo.num_solute_atoms() + j * 3;
  DEBUG(9, "ii = " << ii << " jj = " << jj);
  
  math::Vec const * const pos_i = &conf.current().pos(ii);
  math::Vec const * const pos_j = &conf.current().pos(jj);

  math::Vec * const force_i = &storage.force(ii);
  math::Vec * const force_j = &storage.force(jj);

  double dist6i, e_lj, e_crf, f;

  // O - O

  periodicity.nearest_image(*pos_i, *pos_j, r);

  tx = r(0) - (*pos_i)(0) + (*pos_j)(0);
  ty = r(1) - (*pos_i)(1) + (*pos_j)(1);
  tz = r(2) - (*pos_i)(2) + (*pos_j)(2);
    
  assert(abs2(r) != 0);

  r2[0] = abs2(r);
  r2i[0] = 1.0 / r2[0];
  dist6i = r2i[0] * r2i[0] * r2i[0];
  ri[0] = sqrt(r2i[0]);
  
  e_lj = (2.634129E-6 * dist6i - 2.617346E-3) * dist6i;
  e_crf = 0.82 * 0.82 * 138.935 * (ri[0] - m_crf_2cut3i * r2[0] - m_crf_cut);

  f = (12 * 2.634129E-6 * dist6i - 6 * 2.617346E-3) * dist6i * r2i[0] +
    0.82 * 0.82 * 138.935 * (ri[0] * r2i[0] + m_crf_cut3i);

  fx = f * r(0);
  fy = f * r(1);
  fz = f * r(2);
    
  (*force_i)(0) += fx;
  (*force_j)(0) -= fx;
  (*force_i)(1) += fy;
  (*force_j)(1) -= fy;
  (*force_i)(2) += fz;
  (*force_j)(2) -= fz;

  rx = r(0) -
    conf.special().rel_mol_com_pos(ii)(0) +
    conf.special().rel_mol_com_pos(jj)(0);
  ry = r(1) -
    conf.special().rel_mol_com_pos(ii)(1) +
    conf.special().rel_mol_com_pos(jj)(1);
  rz = r(2) -
    conf.special().rel_mol_com_pos(ii)(2) +
    conf.special().rel_mol_com_pos(jj)(2);

  storage.virial_tensor(0, 0) += rx * fx;
  storage.virial_tensor(0, 1) += rx * fy;
  storage.virial_tensor(0, 2) += rx * fz;
  storage.virial_tensor(1, 0) += ry * fx;
  storage.virial_tensor(1, 1) += ry * fy;
  storage.virial_tensor(1, 2) += ry * fz;
  storage.virial_tensor(2, 0) += rz * fx;
  storage.virial_tensor(2, 1) += rz * fy;
  storage.virial_tensor(2, 2) += rz * fz;
    
  storage.energies.lj_energy[egroup][egroup] += e_lj;
  
  // O - H interactions...

  x[0] = (*pos_i)(0) - (*(pos_j+1))(0) + tx;
  y[0] = (*pos_i)(1) - (*(pos_j+1))(1) + ty;
  z[0] = (*pos_i)(2) - (*(pos_j+1))(2) + tz;

  x[1] = (*pos_i)(0) - (*(pos_j+2))(0) + tx;
  y[1] = (*pos_i)(1) - (*(pos_j+2))(1) + ty;
  z[1] = (*pos_i)(2) - (*(pos_j+2))(2) + tz;

  x[2] = (*(pos_i+1))(0) - (*(pos_j))(0) + tx;
  y[2] = (*(pos_i+1))(1) - (*(pos_j))(1) + ty;
  z[2] = (*(pos_i+1))(2) - (*(pos_j))(2) + tz;

  x[3] = (*(pos_i+2))(0) - (*(pos_j))(0) + tx;
  y[3] = (*(pos_i+2))(1) - (*(pos_j))(1) + ty;
  z[3] = (*(pos_i+2))(2) - (*(pos_j))(2) + tz;

  r2[0] = x[0]*x[0] + y[0]*y[0] + z[0]*z[0];
  r2[1] = x[1]*x[1] + y[1]*y[1] + z[1]*z[1];
  r2[2] = x[2]*x[2] + y[2]*y[2] + z[2]*z[2];
  r2[3] = x[3]*x[3] + y[3]*y[3] + z[3]*z[3];
  
  r2i[0] = 1.0 / r2[0];
  r2i[1] = 1.0 / r2[1];
  r2i[2] = 1.0 / r2[2];
  r2i[3] = 1.0 / r2[3];

  ri[0] = sqrt(r2i[0]);
  ri[1] = sqrt(r2i[1]);
  ri[2] = sqrt(r2i[2]);
  ri[3] = sqrt(r2i[3]);

  e_crf -= 0.82 * 0.41 * 138.935 * (ri[0] + ri[1] + ri[2] + ri[3] -
				    m_crf_2cut3i * (r2[0] + r2[1] + r2[2] + r2[3]) - 4 * m_crf_cut);

  ff[0] = -0.82 * 0.41 * 138.935 * (ri[0] * r2i[0] + m_crf_cut3i);
  ff[1] = -0.82 * 0.41 * 138.935 * (ri[1] * r2i[1] + m_crf_cut3i);
  ff[2] = -0.82 * 0.41 * 138.935 * (ri[2] * r2i[2] + m_crf_cut3i);
  ff[3] = -0.82 * 0.41 * 138.935 * (ri[3] * r2i[3] + m_crf_cut3i);
  
  fx = ff[0] * x[0];
  fy = ff[0] * y[0];
  fz = ff[0] * z[0];
  
  (*force_i)(0) += fx;
  (*(force_j+1))(0) -= fx;
  (*force_i)(1) += fy;
  (*(force_j+1))(1) -= fy;
  (*force_i)(2) += fz;
  (*(force_j+1))(2) -= fz;

  rx = x[0] -
    conf.special().rel_mol_com_pos(ii)(0) +
    conf.special().rel_mol_com_pos(jj+1)(0);
  ry = y[0] - 
    conf.special().rel_mol_com_pos(ii)(1) +
    conf.special().rel_mol_com_pos(jj+1)(1);
  rz = z[0] -
    conf.special().rel_mol_com_pos(ii)(2) +
    conf.special().rel_mol_com_pos(jj+1)(2);

  storage.virial_tensor(0, 0) += rx * fx;
  storage.virial_tensor(0, 1) += rx * fy;
  storage.virial_tensor(0, 2) += rx * fz;
  storage.virial_tensor(1, 0) += ry * fx;
  storage.virial_tensor(1, 1) += ry * fy;
  storage.virial_tensor(1, 2) += ry * fz;
  storage.virial_tensor(2, 0) += rz * fx;
  storage.virial_tensor(2, 1) += rz * fy;
  storage.virial_tensor(2, 2) += rz * fz;

  fx = ff[1] * x[1];
  fy = ff[1] * y[1];
  fz = ff[1] * z[1];
  
  (*force_i)(0) += fx;
  (*(force_j+2))(0) -= fx;
  (*force_i)(1) += fy;
  (*(force_j+2))(1) -= fy;
  (*force_i)(2) += fz;
  (*(force_j+2))(2) -= fz;

  rx = x[1] -
    conf.special().rel_mol_com_pos(ii)(0) +
    conf.special().rel_mol_com_pos(jj+2)(0);
  ry = y[1] - 
    conf.special().rel_mol_com_pos(ii)(1) +
    conf.special().rel_mol_com_pos(jj+2)(1);
  rz = z[1] -
    conf.special().rel_mol_com_pos(ii)(2) +
    conf.special().rel_mol_com_pos(jj+2)(2);

  storage.virial_tensor(0, 0) += rx * fx;
  storage.virial_tensor(0, 1) += rx * fy;
  storage.virial_tensor(0, 2) += rx * fz;
  storage.virial_tensor(1, 0) += ry * fx;
  storage.virial_tensor(1, 1) += ry * fy;
  storage.virial_tensor(1, 2) += ry * fz;
  storage.virial_tensor(2, 0) += rz * fx;
  storage.virial_tensor(2, 1) += rz * fy;
  storage.virial_tensor(2, 2) += rz * fz;

  fx = ff[2] * x[2];
  fy = ff[2] * y[2];
  fz = ff[2] * z[2];
  
  (*(force_i+1))(0) += fx;
  (*(force_j))(0) -= fx;
  (*(force_i+1))(1) += fy;
  (*(force_j))(1) -= fy;
  (*(force_i+1))(2) += fz;
  (*(force_j))(2) -= fz;

  rx = x[2] -
    conf.special().rel_mol_com_pos(ii+1)(0) +
    conf.special().rel_mol_com_pos(jj)(0);
  ry = y[2] - 
    conf.special().rel_mol_com_pos(ii+1)(1) +
    conf.special().rel_mol_com_pos(jj)(1);
  rz = z[2] -
    conf.special().rel_mol_com_pos(ii+1)(2) +
    conf.special().rel_mol_com_pos(jj)(2);

  storage.virial_tensor(0, 0) += rx * fx;
  storage.virial_tensor(0, 1) += rx * fy;
  storage.virial_tensor(0, 2) += rx * fz;
  storage.virial_tensor(1, 0) += ry * fx;
  storage.virial_tensor(1, 1) += ry * fy;
  storage.virial_tensor(1, 2) += ry * fz;
  storage.virial_tensor(2, 0) += rz * fx;
  storage.virial_tensor(2, 1) += rz * fy;
  storage.virial_tensor(2, 2) += rz * fz;

  fx = ff[3] * x[3];
  fy = ff[3] * y[3];
  fz = ff[3] * z[3];
  
  (*(force_i+2))(0) += fx;
  (*(force_j))(0) -= fx;
  (*(force_i+2))(1) += fy;
  (*(force_j))(1) -= fy;
  (*(force_i+2))(2) += fz;
  (*(force_j))(2) -= fz;

  rx = x[3] -
    conf.special().rel_mol_com_pos(ii+2)(0) +
    conf.special().rel_mol_com_pos(jj)(0);
  ry = y[3] - 
    conf.special().rel_mol_com_pos(ii+2)(1) +
    conf.special().rel_mol_com_pos(jj)(1);
  rz = z[3] -
    conf.special().rel_mol_com_pos(ii+2)(2) +
    conf.special().rel_mol_com_pos(jj)(2);

  storage.virial_tensor(0, 0) += rx * fx;
  storage.virial_tensor(0, 1) += rx * fy;
  storage.virial_tensor(0, 2) += rx * fz;
  storage.virial_tensor(1, 0) += ry * fx;
  storage.virial_tensor(1, 1) += ry * fy;
  storage.virial_tensor(1, 2) += ry * fz;
  storage.virial_tensor(2, 0) += rz * fx;
  storage.virial_tensor(2, 1) += rz * fy;
  storage.virial_tensor(2, 2) += rz * fz;

  // H - H interactions...

  x[0] = (*(pos_i+1))(0) - (*(pos_j+1))(0) + tx;
  y[0] = (*(pos_i+1))(1) - (*(pos_j+1))(1) + ty;
  z[0] = (*(pos_i+1))(2) - (*(pos_j+1))(2) + tz;

  x[1] = (*(pos_i+1))(0) - (*(pos_j+2))(0) + tx;
  y[1] = (*(pos_i+1))(1) - (*(pos_j+2))(1) + ty;
  z[1] = (*(pos_i+1))(2) - (*(pos_j+2))(2) + tz;

  x[2] = (*(pos_i+2))(0) - (*(pos_j+1))(0) + tx;
  y[2] = (*(pos_i+2))(1) - (*(pos_j+1))(1) + ty;
  z[2] = (*(pos_i+2))(2) - (*(pos_j+1))(2) + tz;

  x[3] = (*(pos_i+2))(0) - (*(pos_j+2))(0) + tx;
  y[3] = (*(pos_i+2))(1) - (*(pos_j+2))(1) + ty;
  z[3] = (*(pos_i+2))(2) - (*(pos_j+2))(2) + tz;

  r2[0] = x[0]*x[0] + y[0]*y[0] + z[0]*z[0];
  r2[1] = x[1]*x[1] + y[1]*y[1] + z[1]*z[1];
  r2[2] = x[2]*x[2] + y[2]*y[2] + z[2]*z[2];
  r2[3] = x[3]*x[3] + y[3]*y[3] + z[3]*z[3];
  
  r2i[0] = 1.0 / r2[0];
  r2i[1] = 1.0 / r2[1];
  r2i[2] = 1.0 / r2[2];
  r2i[3] = 1.0 / r2[3];

  ri[0] = sqrt(r2i[0]);
  ri[1] = sqrt(r2i[1]);
  ri[2] = sqrt(r2i[2]);
  ri[3] = sqrt(r2i[3]);

  e_crf += 0.41 * 0.41 * 138.935 * (ri[0] + ri[1] + ri[2] + ri[3] -
				    m_crf_2cut3i * (r2[0] + r2[1] + r2[2] + r2[3]) - 4 * m_crf_cut);

  ff[0] = 0.41 * 0.41 * 138.935 * (ri[0] * r2i[0] + m_crf_cut3i);
  ff[1] = 0.41 * 0.41 * 138.935 * (ri[1] * r2i[1] + m_crf_cut3i);
  ff[2] = 0.41 * 0.41 * 138.935 * (ri[2] * r2i[2] + m_crf_cut3i);
  ff[3] = 0.41 * 0.41 * 138.935 * (ri[3] * r2i[3] + m_crf_cut3i);
  
  fx = ff[0] * x[0];
  fy = ff[0] * y[0];
  fz = ff[0] * z[0];
  
  (*(force_i+1))(0) += fx;
  (*(force_j+1))(0) -= fx;
  (*(force_i+1))(1) += fy;
  (*(force_j+1))(1) -= fy;
  (*(force_i+1))(2) += fz;
  (*(force_j+1))(2) -= fz;

  rx = x[0] -
    conf.special().rel_mol_com_pos(ii+1)(0) +
    conf.special().rel_mol_com_pos(jj+1)(0);
  ry = y[0] - 
    conf.special().rel_mol_com_pos(ii+1)(1) +
    conf.special().rel_mol_com_pos(jj+1)(1);
  rz = z[0] -
    conf.special().rel_mol_com_pos(ii+1)(2) +
    conf.special().rel_mol_com_pos(jj+1)(2);

  storage.virial_tensor(0, 0) += rx * fx;
  storage.virial_tensor(0, 1) += rx * fy;
  storage.virial_tensor(0, 2) += rx * fz;
  storage.virial_tensor(1, 0) += ry * fx;
  storage.virial_tensor(1, 1) += ry * fy;
  storage.virial_tensor(1, 2) += ry * fz;
  storage.virial_tensor(2, 0) += rz * fx;
  storage.virial_tensor(2, 1) += rz * fy;
  storage.virial_tensor(2, 2) += rz * fz;

  fx = ff[1] * x[1];
  fy = ff[1] * y[1];
  fz = ff[1] * z[1];
  
  (*(force_i+1))(0) += fx;
  (*(force_j+2))(0) -= fx;
  (*(force_i+1))(1) += fy;
  (*(force_j+2))(1) -= fy;
  (*(force_i+1))(2) += fz;
  (*(force_j+2))(2) -= fz;

  rx = x[1] -
    conf.special().rel_mol_com_pos(ii+1)(0) +
    conf.special().rel_mol_com_pos(jj+2)(0);
  ry = y[1] - 
    conf.special().rel_mol_com_pos(ii+1)(1) +
    conf.special().rel_mol_com_pos(jj+2)(1);
  rz = z[1] -
    conf.special().rel_mol_com_pos(ii+1)(2) +
    conf.special().rel_mol_com_pos(jj+2)(2);

  storage.virial_tensor(0, 0) += rx * fx;
  storage.virial_tensor(0, 1) += rx * fy;
  storage.virial_tensor(0, 2) += rx * fz;
  storage.virial_tensor(1, 0) += ry * fx;
  storage.virial_tensor(1, 1) += ry * fy;
  storage.virial_tensor(1, 2) += ry * fz;
  storage.virial_tensor(2, 0) += rz * fx;
  storage.virial_tensor(2, 1) += rz * fy;
  storage.virial_tensor(2, 2) += rz * fz;

  fx = ff[2] * x[2];
  fy = ff[2] * y[2];
  fz = ff[2] * z[2];
  
  (*(force_i+2))(0) += fx;
  (*(force_j+1))(0) -= fx;
  (*(force_i+2))(1) += fy;
  (*(force_j+1))(1) -= fy;
  (*(force_i+2))(2) += fz;
  (*(force_j+1))(2) -= fz;

  rx = x[2] -
    conf.special().rel_mol_com_pos(ii+2)(0) +
    conf.special().rel_mol_com_pos(jj+1)(0);
  ry = y[2] - 
    conf.special().rel_mol_com_pos(ii+2)(1) +
    conf.special().rel_mol_com_pos(jj+1)(1);
  rz = z[2] -
    conf.special().rel_mol_com_pos(ii+2)(2) +
    conf.special().rel_mol_com_pos(jj+1)(2);

  storage.virial_tensor(0, 0) += rx * fx;
  storage.virial_tensor(0, 1) += rx * fy;
  storage.virial_tensor(0, 2) += rx * fz;
  storage.virial_tensor(1, 0) += ry * fx;
  storage.virial_tensor(1, 1) += ry * fy;
  storage.virial_tensor(1, 2) += ry * fz;
  storage.virial_tensor(2, 0) += rz * fx;
  storage.virial_tensor(2, 1) += rz * fy;
  storage.virial_tensor(2, 2) += rz * fz;

  fx = ff[3] * x[3];
  fy = ff[3] * y[3];
  fz = ff[3] * z[3];
  
  (*(force_i+2))(0) += fx;
  (*(force_j+2))(0) -= fx;
  (*(force_i+2))(1) += fy;
  (*(force_j+2))(1) -= fy;
  (*(force_i+2))(2) += fz;
  (*(force_j+2))(2) -= fz;

  rx = x[3] -
    conf.special().rel_mol_com_pos(ii+2)(0) +
    conf.special().rel_mol_com_pos(jj+2)(0);
  ry = y[3] - 
    conf.special().rel_mol_com_pos(ii+2)(1) +
    conf.special().rel_mol_com_pos(jj+2)(1);
  rz = z[3] -
    conf.special().rel_mol_com_pos(ii+2)(2) +
    conf.special().rel_mol_com_pos(jj+2)(2);

  storage.virial_tensor(0, 0) += rx * fx;
  storage.virial_tensor(0, 1) += rx * fy;
  storage.virial_tensor(0, 2) += rx * fz;
  storage.virial_tensor(1, 0) += ry * fx;
  storage.virial_tensor(1, 1) += ry * fy;
  storage.virial_tensor(1, 2) += ry * fz;
  storage.virial_tensor(2, 0) += rz * fx;
  storage.virial_tensor(2, 1) += rz * fy;
  storage.virial_tensor(2, 2) += rz * fz;
    
  storage.energies.crf_energy[egroup][egroup] += e_crf;
}
