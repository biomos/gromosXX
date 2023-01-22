/**
 * @file solvent_innerloop.cc
 * special loops for fast treatment of solvent
 */

template<typename t_nonbonded_spec>
inline void 
interaction::Nonbonded_Innerloop<t_nonbonded_spec>::solvent_innerloop
(
 topology::Topology & topo,
 interaction::Nonbonded_Innerloop<t_nonbonded_spec>::solvent_pair_parameter * pair_parameter,
 configuration::Configuration & conf,
 const unsigned int num_solvent_atoms,
 const int i,
 const int j,
 Storage & storage,
 math::Periodicity<t_nonbonded_spec::boundary_type> const & periodicity,
 unsigned int eps)
{
  math::Vec r;
  
  // only one energy group
  const int egroup = topo.atom_energy_group(topo.num_solute_atoms());
  DEBUG(8, "\tspc pair\t" << i << "\t" << j << " egroup " << egroup);
  DEBUG(1, "Fast solvent");
  math::Vec const * const pos_i = &conf.current().pos(i);
  math::Vec * const force_i = &storage.force(i);

  math::Vec const * const pos_j = &conf.current().pos(j);
  math::Vec * const force_j = &storage.force(j);
  DEBUG(9, "i = " << i << " j = " << j);    

  // first atom vs. first atom
  periodicity.nearest_image(*pos_i, *pos_j, r);
  const double tx = r(0) - (*pos_i)(0) + (*pos_j)(0);
  const double ty = r(1) - (*pos_i)(1) + (*pos_j)(1);
  const double tz = r(2) - (*pos_i)(2) + (*pos_j)(2);
  
  double e_lj = 0.0, e_crf = 0.0;
  
  // loop over atom pairs
  for(unsigned int param = 0, atom_i = 0; atom_i < num_solvent_atoms; ++atom_i) {
    const double xi = (*(pos_i+atom_i))(0) + tx;
    const double yi = (*(pos_i+atom_i))(1) + ty;
    const double zi = (*(pos_i+atom_i))(2) + tz;
    for(unsigned int atom_j = 0; atom_j < num_solvent_atoms; ++atom_j, ++param) {
      DEBUG(15, "\tatoms: i: " << atom_i << " j: " << atom_j);
      const double x = xi - (*(pos_j+atom_j))(0);
      const double y = yi - (*(pos_j+atom_j))(1);
      const double z = zi - (*(pos_j+atom_j))(2);
      
      const double r2 = x * x + y * y + z * z;
      DEBUG(15, "\tr2: " << r2);
      assert(r2 != 0.0);
      const double r2i = 1.0 / r2;
      const double ri = sqrt(r2i);
      const double dist6i = r2i * r2i * r2i;
      const double dist6i_c12 = pair_parameter[param].c12 * dist6i;

      e_lj += (dist6i_c12 - pair_parameter[param].c6) * dist6i;
      e_crf += pair_parameter[param].q * (ri - m_crf_2cut3i[eps] * r2 - m_crf_cut[eps]);

      const double f = (dist6i_c12 + dist6i_c12 - pair_parameter[param].c6) * 6.0 
              * dist6i * r2i + pair_parameter[param].q * (ri * r2i + m_crf_cut3i[eps]);

      const double fx = f * x;
      const double fy = f * y;
      const double fz = f * z;

      (*(force_i+atom_i))(0) += fx;
      (*(force_j+atom_j))(0) -= fx;
      (*(force_i+atom_i))(1) += fy;
      (*(force_j+atom_j))(1) -= fy;
      (*(force_i+atom_i))(2) += fz;
      (*(force_j+atom_j))(2) -= fz;

      storage.virial_tensor(0, 0) += x * fx;
      storage.virial_tensor(0, 1) += x * fy;
      storage.virial_tensor(0, 2) += x * fz;
      storage.virial_tensor(1, 0) += y * fx;
      storage.virial_tensor(1, 1) += y * fy;
      storage.virial_tensor(1, 2) += y * fz;
      storage.virial_tensor(2, 0) += z * fx;
      storage.virial_tensor(2, 1) += z * fy;
      storage.virial_tensor(2, 2) += z * fz;
    }
  }
  storage.energies.lj_energy[egroup][egroup] += e_lj;
  storage.energies.crf_energy[egroup][egroup] += e_crf;
}


template<typename t_nonbonded_spec>
inline void 
interaction::Nonbonded_Innerloop<t_nonbonded_spec>::spc_innerloop
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 int i,
 int j,
 Storage & storage,
 math::Periodicity<t_nonbonded_spec::boundary_type> const & periodicity,
 unsigned int eps)
{
  math::Vec r;
  double x[4], y[4], z[4], r2[4], r2i[4], ri[4], ff[4], tx = 0.0, ty = 0.0, tz = 0.0, fx = 0.0, fy = 0.0, fz = 0.0;
  // , rx, ry, rz;

  // only one energy group
  const int egroup = topo.atom_energy_group(topo.num_solute_atoms());
  DEBUG(8, "\tspc pair\t" << i << "\t" << j << " egroup " << egroup);
  
  double dist6i = 0.0, e_lj = 0.0, e_crf = 0.0, f = 0.0;
  
  const int ii = i;
  math::Vec const * const pos_i = &conf.current().pos(ii);
  math::Vec * const force_i = &storage.force(ii);

  const int jj = j;
  math::Vec const * const pos_j = &conf.current().pos(jj);
  math::Vec * const force_j = &storage.force(jj);
  DEBUG(9, "ii = " << ii << " jj = " << jj);    

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
  e_crf = 0.82 * 0.82 * 138.9354 * (ri[0] - m_crf_2cut3i[eps] * r2[0] - m_crf_cut[eps]);
      
  f = (12 * 2.634129E-6 * dist6i - 6 * 2.617346E-3) * dist6i * r2i[0] +
    0.82 * 0.82 * 138.9354 * (ri[0] * r2i[0] + m_crf_cut3i[eps]);
  
  DEBUG(10, "r: " << sqrt(r2[0]) << " r2: " << r2[0]);
  DEBUG(10, "e_lj: " << e_lj << " c_crf: " << e_crf << " f: " << f);
      
  fx = f * r(0);
  fy = f * r(1);
  fz = f * r(2);
      
  (*force_i)(0) += fx;
  (*force_j)(0) -= fx;
  (*force_i)(1) += fy;
  (*force_j)(1) -= fy;
  (*force_i)(2) += fz;
  (*force_j)(2) -= fz;
      
  /*
  rx = r(0) -
    conf.special().rel_mol_com_pos(ii)(0) +
    conf.special().rel_mol_com_pos(jj)(0);
  ry = r(1) -
    conf.special().rel_mol_com_pos(ii)(1) +
    conf.special().rel_mol_com_pos(jj)(1);
  rz = r(2) -
    conf.special().rel_mol_com_pos(ii)(2) +
    conf.special().rel_mol_com_pos(jj)(2);
  */
  
  storage.virial_tensor(0, 0) += r(0) * fx;
  storage.virial_tensor(0, 1) += r(0) * fy;
  storage.virial_tensor(0, 2) += r(0) * fz;
  storage.virial_tensor(1, 0) += r(1) * fx;
  storage.virial_tensor(1, 1) += r(1) * fy;
  storage.virial_tensor(1, 2) += r(1) * fz;
  storage.virial_tensor(2, 0) += r(2) * fx;
  storage.virial_tensor(2, 1) += r(2) * fy;
  storage.virial_tensor(2, 2) += r(2) * fz;
      
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
      
  e_crf -= 0.82 * 0.41 * 138.9354 * (ri[0] + ri[1] + ri[2] + ri[3] -
				    m_crf_2cut3i[eps] * (r2[0] + r2[1] + r2[2] + r2[3]) - 4 * m_crf_cut[eps]);
  
  DEBUG(10, "e_crf: " << e_crf);
      
  ff[0] = -0.82 * 0.41 * 138.9354 * (ri[0] * r2i[0] + m_crf_cut3i[eps]);
  ff[1] = -0.82 * 0.41 * 138.9354 * (ri[1] * r2i[1] + m_crf_cut3i[eps]);
  ff[2] = -0.82 * 0.41 * 138.9354 * (ri[2] * r2i[2] + m_crf_cut3i[eps]);
  ff[3] = -0.82 * 0.41 * 138.9354 * (ri[3] * r2i[3] + m_crf_cut3i[eps]);
  
  DEBUG(10, "f: " << ff[0]);
  DEBUG(10, "f: " << ff[1]);
  DEBUG(10, "f: " << ff[2]);
  DEBUG(10, "f: " << ff[3]);
      
  fx = ff[0] * x[0];
  fy = ff[0] * y[0];
  fz = ff[0] * z[0];
      
  (*force_i)(0) += fx;
  (*(force_j+1))(0) -= fx;
  (*force_i)(1) += fy;
  (*(force_j+1))(1) -= fy;
  (*force_i)(2) += fz;
  (*(force_j+1))(2) -= fz;
      
  /*
  rx = x[0] -
    conf.special().rel_mol_com_pos(ii)(0) +
    conf.special().rel_mol_com_pos(jj+1)(0);
  ry = y[0] - 
    conf.special().rel_mol_com_pos(ii)(1) +
    conf.special().rel_mol_com_pos(jj+1)(1);
  rz = z[0] -
    conf.special().rel_mol_com_pos(ii)(2) +
    conf.special().rel_mol_com_pos(jj+1)(2);
  */    
  
  storage.virial_tensor(0, 0) += x[0] * fx;
  storage.virial_tensor(0, 1) += x[0] * fy;
  storage.virial_tensor(0, 2) += x[0] * fz;
  storage.virial_tensor(1, 0) += y[0] * fx;
  storage.virial_tensor(1, 1) += y[0] * fy;
  storage.virial_tensor(1, 2) += y[0] * fz;
  storage.virial_tensor(2, 0) += z[0] * fx;
  storage.virial_tensor(2, 1) += z[0] * fy;
  storage.virial_tensor(2, 2) += z[0] * fz;
      
  fx = ff[1] * x[1];
  fy = ff[1] * y[1];
  fz = ff[1] * z[1];
      
  (*force_i)(0) += fx;
  (*(force_j+2))(0) -= fx;
  (*force_i)(1) += fy;
  (*(force_j+2))(1) -= fy;
  (*force_i)(2) += fz;
  (*(force_j+2))(2) -= fz;
    
  /*
  rx = x[1] -
    conf.special().rel_mol_com_pos(ii)(0) +
    conf.special().rel_mol_com_pos(jj+2)(0);
  ry = y[1] - 
    conf.special().rel_mol_com_pos(ii)(1) +
    conf.special().rel_mol_com_pos(jj+2)(1);
  rz = z[1] -
    conf.special().rel_mol_com_pos(ii)(2) +
    conf.special().rel_mol_com_pos(jj+2)(2);
  */
  
  storage.virial_tensor(0, 0) += x[1] * fx;
  storage.virial_tensor(0, 1) += x[1] * fy;
  storage.virial_tensor(0, 2) += x[1] * fz;
  storage.virial_tensor(1, 0) += y[1] * fx;
  storage.virial_tensor(1, 1) += y[1] * fy;
  storage.virial_tensor(1, 2) += y[1] * fz;
  storage.virial_tensor(2, 0) += z[1] * fx;
  storage.virial_tensor(2, 1) += z[1] * fy;
  storage.virial_tensor(2, 2) += z[1] * fz;
      
  fx = ff[2] * x[2];
  fy = ff[2] * y[2];
  fz = ff[2] * z[2];
      
  (*(force_i+1))(0) += fx;
  (*(force_j))(0) -= fx;
  (*(force_i+1))(1) += fy;
  (*(force_j))(1) -= fy;
  (*(force_i+1))(2) += fz;
  (*(force_j))(2) -= fz;
      
  /*
  rx = x[2] -
    conf.special().rel_mol_com_pos(ii+1)(0) +
    conf.special().rel_mol_com_pos(jj)(0);
  ry = y[2] - 
    conf.special().rel_mol_com_pos(ii+1)(1) +
    conf.special().rel_mol_com_pos(jj)(1);
  rz = z[2] -
    conf.special().rel_mol_com_pos(ii+1)(2) +
    conf.special().rel_mol_com_pos(jj)(2);
  */
  
  storage.virial_tensor(0, 0) += x[2] * fx;
  storage.virial_tensor(0, 1) += x[2] * fy;
  storage.virial_tensor(0, 2) += x[2] * fz;
  storage.virial_tensor(1, 0) += y[2] * fx;
  storage.virial_tensor(1, 1) += y[2] * fy;
  storage.virial_tensor(1, 2) += y[2] * fz;
  storage.virial_tensor(2, 0) += z[2] * fx;
  storage.virial_tensor(2, 1) += z[2] * fy;
  storage.virial_tensor(2, 2) += z[2] * fz;
      
  fx = ff[3] * x[3];
  fy = ff[3] * y[3];
  fz = ff[3] * z[3];
      
  (*(force_i+2))(0) += fx;
  (*(force_j))(0) -= fx;
  (*(force_i+2))(1) += fy;
  (*(force_j))(1) -= fy;
  (*(force_i+2))(2) += fz;
  (*(force_j))(2) -= fz;

  /*
  rx = x[3] -
    conf.special().rel_mol_com_pos(ii+2)(0) +
    conf.special().rel_mol_com_pos(jj)(0);
  ry = y[3] - 
    conf.special().rel_mol_com_pos(ii+2)(1) +
    conf.special().rel_mol_com_pos(jj)(1);
  rz = z[3] -
    conf.special().rel_mol_com_pos(ii+2)(2) +
    conf.special().rel_mol_com_pos(jj)(2);
  */
  
  storage.virial_tensor(0, 0) += x[3] * fx;
  storage.virial_tensor(0, 1) += x[3] * fy;
  storage.virial_tensor(0, 2) += x[3] * fz;
  storage.virial_tensor(1, 0) += y[3] * fx;
  storage.virial_tensor(1, 1) += y[3] * fy;
  storage.virial_tensor(1, 2) += y[3] * fz;
  storage.virial_tensor(2, 0) += z[3] * fx;
  storage.virial_tensor(2, 1) += z[3] * fy;
  storage.virial_tensor(2, 2) += z[3] * fz;
      
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

  e_crf += 0.41 * 0.41 * 138.9354 * (ri[0] + ri[1] + ri[2] + ri[3] -
				    m_crf_2cut3i[eps] * (r2[0] + r2[1] + r2[2] + r2[3]) - 4 * m_crf_cut[eps]);
  
  DEBUG(10, "e_crf: " << e_crf);

  ff[0] = 0.41 * 0.41 * 138.9354 * (ri[0] * r2i[0] + m_crf_cut3i[eps]);
  ff[1] = 0.41 * 0.41 * 138.9354 * (ri[1] * r2i[1] + m_crf_cut3i[eps]);
  ff[2] = 0.41 * 0.41 * 138.9354 * (ri[2] * r2i[2] + m_crf_cut3i[eps]);
  ff[3] = 0.41 * 0.41 * 138.9354 * (ri[3] * r2i[3] + m_crf_cut3i[eps]);
  
  DEBUG(10, "f: " << ff[0]);
  DEBUG(10, "f: " << ff[1]);
  DEBUG(10, "f: " << ff[2]);
  DEBUG(10, "f: " << ff[3]);
  
  fx = ff[0] * x[0];
  fy = ff[0] * y[0];
  fz = ff[0] * z[0];
  
  (*(force_i+1))(0) += fx;
  (*(force_j+1))(0) -= fx;
  (*(force_i+1))(1) += fy;
  (*(force_j+1))(1) -= fy;
  (*(force_i+1))(2) += fz;
  (*(force_j+1))(2) -= fz;

  /*
  rx = x[0] -
    conf.special().rel_mol_com_pos(ii+1)(0) +
    conf.special().rel_mol_com_pos(jj+1)(0);
  ry = y[0] - 
    conf.special().rel_mol_com_pos(ii+1)(1) +
    conf.special().rel_mol_com_pos(jj+1)(1);
  rz = z[0] -
    conf.special().rel_mol_com_pos(ii+1)(2) +
    conf.special().rel_mol_com_pos(jj+1)(2);
  */

  storage.virial_tensor(0, 0) += x[0] * fx;
  storage.virial_tensor(0, 1) += x[0] * fy;
  storage.virial_tensor(0, 2) += x[0] * fz;
  storage.virial_tensor(1, 0) += y[0] * fx;
  storage.virial_tensor(1, 1) += y[0] * fy;
  storage.virial_tensor(1, 2) += y[0] * fz;
  storage.virial_tensor(2, 0) += z[0] * fx;
  storage.virial_tensor(2, 1) += z[0] * fy;
  storage.virial_tensor(2, 2) += z[0] * fz;

  fx = ff[1] * x[1];
  fy = ff[1] * y[1];
  fz = ff[1] * z[1];
  
  (*(force_i+1))(0) += fx;
  (*(force_j+2))(0) -= fx;
  (*(force_i+1))(1) += fy;
  (*(force_j+2))(1) -= fy;
  (*(force_i+1))(2) += fz;
  (*(force_j+2))(2) -= fz;

  /*
  rx = x[1] -
    conf.special().rel_mol_com_pos(ii+1)(0) +
    conf.special().rel_mol_com_pos(jj+2)(0);
  ry = y[1] - 
    conf.special().rel_mol_com_pos(ii+1)(1) +
    conf.special().rel_mol_com_pos(jj+2)(1);
  rz = z[1] -
    conf.special().rel_mol_com_pos(ii+1)(2) +
    conf.special().rel_mol_com_pos(jj+2)(2);
  */

  storage.virial_tensor(0, 0) += x[1] * fx;
  storage.virial_tensor(0, 1) += x[1] * fy;
  storage.virial_tensor(0, 2) += x[1] * fz;
  storage.virial_tensor(1, 0) += y[1] * fx;
  storage.virial_tensor(1, 1) += y[1] * fy;
  storage.virial_tensor(1, 2) += y[1] * fz;
  storage.virial_tensor(2, 0) += z[1] * fx;
  storage.virial_tensor(2, 1) += z[1] * fy;
  storage.virial_tensor(2, 2) += z[1] * fz;

  fx = ff[2] * x[2];
  fy = ff[2] * y[2];
  fz = ff[2] * z[2];
  
  (*(force_i+2))(0) += fx;
  (*(force_j+1))(0) -= fx;
  (*(force_i+2))(1) += fy;
  (*(force_j+1))(1) -= fy;
  (*(force_i+2))(2) += fz;
  (*(force_j+1))(2) -= fz;

  /*
  rx = x[2] -
    conf.special().rel_mol_com_pos(ii+2)(0) +
    conf.special().rel_mol_com_pos(jj+1)(0);
  ry = y[2] - 
    conf.special().rel_mol_com_pos(ii+2)(1) +
    conf.special().rel_mol_com_pos(jj+1)(1);
  rz = z[2] -
    conf.special().rel_mol_com_pos(ii+2)(2) +
    conf.special().rel_mol_com_pos(jj+1)(2);
  */

  storage.virial_tensor(0, 0) += x[2] * fx;
  storage.virial_tensor(0, 1) += x[2] * fy;
  storage.virial_tensor(0, 2) += x[2] * fz;
  storage.virial_tensor(1, 0) += y[2] * fx;
  storage.virial_tensor(1, 1) += y[2] * fy;
  storage.virial_tensor(1, 2) += y[2] * fz;
  storage.virial_tensor(2, 0) += z[2] * fx;
  storage.virial_tensor(2, 1) += z[2] * fy;
  storage.virial_tensor(2, 2) += z[2] * fz;

  fx = ff[3] * x[3];
  fy = ff[3] * y[3];
  fz = ff[3] * z[3];
  
  (*(force_i+2))(0) += fx;
  (*(force_j+2))(0) -= fx;
  (*(force_i+2))(1) += fy;
  (*(force_j+2))(1) -= fy;
  (*(force_i+2))(2) += fz;
  (*(force_j+2))(2) -= fz;

  /*
  rx = x[3] -
    conf.special().rel_mol_com_pos(ii+2)(0) +
    conf.special().rel_mol_com_pos(jj+2)(0);
  ry = y[3] - 
    conf.special().rel_mol_com_pos(ii+2)(1) +
    conf.special().rel_mol_com_pos(jj+2)(1);
  rz = z[3] -
    conf.special().rel_mol_com_pos(ii+2)(2) +
    conf.special().rel_mol_com_pos(jj+2)(2);
  */

  storage.virial_tensor(0, 0) += x[3] * fx;
  storage.virial_tensor(0, 1) += x[3] * fy;
  storage.virial_tensor(0, 2) += x[3] * fz;
  storage.virial_tensor(1, 0) += y[3] * fx;
  storage.virial_tensor(1, 1) += y[3] * fy;
  storage.virial_tensor(1, 2) += y[3] * fz;
  storage.virial_tensor(2, 0) += z[3] * fx;
  storage.virial_tensor(2, 1) += z[3] * fy;
  storage.virial_tensor(2, 2) += z[3] * fz;
    
  storage.energies.crf_energy[egroup][egroup] += e_crf;
      
}


template<typename t_nonbonded_spec>
inline void 
interaction::Nonbonded_Innerloop<t_nonbonded_spec>::spc_innerloop
(
  double &e_lj,
  double &e_crf,
  double dist6i,
  double f[9],
  double r2[9],
  double r2i[9],
  double ri[9],
  unsigned int eps
)
{  
  // O - O
      
  e_lj += (2.634129E-6 * dist6i - 2.617346E-3) * dist6i;
  e_crf += 0.82 * 0.82 * 138.9354 * (ri[0] - m_crf_2cut3i[eps] * r2[0] - m_crf_cut[eps]);
      
  f[0] = (12 * 2.634129E-6 * dist6i - 6 * 2.617346E-3) * dist6i * r2i[0] +
    0.82 * 0.82 * 138.9354 * (ri[0] * r2i[0] + m_crf_cut3i[eps]);  
  
  DEBUG(10, "e_lj: " << e_lj << " c_crf: " << e_crf << " f: " << f[0]);
  
  // O - H interactions...
      
  e_crf -= 0.82 * 0.41 * 138.9354 * (ri[1] + ri[2] + ri[3] + ri[4] -
				    m_crf_2cut3i[eps] * (r2[1] + r2[2] + r2[3] + r2[4]) - 4 * m_crf_cut[eps]);
  
  DEBUG(10, "e_crf: " << e_crf);
      
  f[1] = -0.82 * 0.41 * 138.9354 * (ri[1] * r2i[1] + m_crf_cut3i[eps]);
  f[2] = -0.82 * 0.41 * 138.9354 * (ri[2] * r2i[2] + m_crf_cut3i[eps]);
  f[3] = -0.82 * 0.41 * 138.9354 * (ri[3] * r2i[3] + m_crf_cut3i[eps]);
  f[4] = -0.82 * 0.41 * 138.9354 * (ri[4] * r2i[4] + m_crf_cut3i[eps]);
  
  DEBUG(10, "f: " << f[1]);
  DEBUG(10, "f: " << f[2]);
  DEBUG(10, "f: " << f[3]);
  DEBUG(10, "f: " << f[4]);
      
  // H - H interactions...

  e_crf += 0.41 * 0.41 * 138.9354 * (ri[5] + ri[6] + ri[7] + ri[8] -
				    m_crf_2cut3i[eps] * (r2[5] + r2[6] + r2[7] + r2[8]) - 4 * m_crf_cut[eps]);
  
  DEBUG(10, "e_crf: " << e_crf);

  f[5] = 0.41 * 0.41 * 138.9354 * (ri[5] * r2i[5] + m_crf_cut3i[eps]);
  f[6] = 0.41 * 0.41 * 138.9354 * (ri[6] * r2i[6] + m_crf_cut3i[eps]);
  f[7] = 0.41 * 0.41 * 138.9354 * (ri[7] * r2i[7] + m_crf_cut3i[eps]);
  f[8] = 0.41 * 0.41 * 138.9354 * (ri[8] * r2i[8] + m_crf_cut3i[eps]);
  
  DEBUG(10, "f: " << f[5]);
  DEBUG(10, "f: " << f[6]);
  DEBUG(10, "f: " << f[7]);
  DEBUG(10, "f: " << f[8]);
  
}
 
template<typename t_nonbonded_spec>
inline void 
interaction::Nonbonded_Innerloop<t_nonbonded_spec>::shortrange_spc_table_innerloop
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
  unsigned int r2_int[4];
  double x[4], y[4], z[4], r2[4], r2_tab[4], ff[4], tx = 0.0, ty = 0.0, tz = 0.0, fx = 0.0, fy = 0.0, fz = 0.0;

  // only one energy group
  const int egroup = topo.atom_energy_group(topo.num_solute_atoms());
  DEBUG(8, "\ttable shortrange");
  DEBUG(8, "\tspc pair\t" << i << "\t" << j << " egroup " << egroup);
  
  double e_lj = 0.0, e_crf = 0.0;
  
  const int ii = i;
  math::Vec const * const pos_i = &conf.current().pos(ii);
  math::Vec * const force_i = &storage.force(ii);

  const int jj = j;
  math::Vec const * const pos_j = &conf.current().pos(jj);
  math::Vec * const force_j = &storage.force(jj);
  DEBUG(9, "ii = " << ii << " jj = " << jj);    

  // O - O
  periodicity.nearest_image(*pos_i, *pos_j, r);
      
  tx = r(0) - (*pos_i)(0) + (*pos_j)(0);
  ty = r(1) - (*pos_i)(1) + (*pos_j)(1);
  tz = r(2) - (*pos_i)(2) + (*pos_j)(2);
      
  assert(abs2(r) != 0);
      
  r2[0] = abs2(r);
  
  struct lj_crf_values {
    double e_lj, e_lj_diff, e_crf, e_crf_diff, f, f_diff;
  };
  
  struct crf_values {
    double e_crf, e_crf_diff, f, f_diff;
  };
  
  #define SHORTRANGE "shortrange"
  #include "spc_table.h"
  #undef SHORTRANGE

  static const double to_table = table_size / data_range;
  const crf_values * values[4];
  
  // calculate the grid point
  r2_tab[0] = r2[0] * to_table;    
  // round the grid point to the lower integer
  r2_int[0] = int(r2_tab[0]);
  assert(r2_int[0] < table_size);
  // safe the offset from the grid point for linear interpolation
  r2_tab[0] -= r2_int[0];
  DEBUG(10, "r: " << sqrt(r2[0]) << " r2: " << r2[0] << " r2_tab: " << r2_tab[0] << " r2_int: " << r2_int[0]);
  
  // the the data from the grid point
  const lj_crf_values & val = OO_table[r2_int[0]];
  
  // the _diff variables hold the gradient for the interpolation.
  e_lj = val.e_lj + val.e_lj_diff * r2_tab[0];
  e_crf = val.e_crf + val.e_crf_diff * r2_tab[0];
  ff[0] = val.f + val.f_diff * r2_tab[0];
  
  DEBUG(10, "e_lj: " << e_lj << " c_crf: " << e_crf << " f: " << ff[0]);

  fx = ff[0] * r(0);
  fy = ff[0] * r(1);
  fz = ff[0] * r(2);
      
  (*force_i)(0) += fx;
  (*force_j)(0) -= fx;
  (*force_i)(1) += fy;
  (*force_j)(1) -= fy;
  (*force_i)(2) += fz;
  (*force_j)(2) -= fz;
  
  storage.virial_tensor(0, 0) += r(0) * fx;
  storage.virial_tensor(0, 1) += r(0) * fy;
  storage.virial_tensor(0, 2) += r(0) * fz;
  storage.virial_tensor(1, 0) += r(1) * fx;
  storage.virial_tensor(1, 1) += r(1) * fy;
  storage.virial_tensor(1, 2) += r(1) * fz;
  storage.virial_tensor(2, 0) += r(2) * fx;
  storage.virial_tensor(2, 1) += r(2) * fy;
  storage.virial_tensor(2, 2) += r(2) * fz;
      
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

  r2_tab[0] = r2[0] * to_table;
  r2_tab[1] = r2[1] * to_table;
  r2_tab[2] = r2[2] * to_table;
  r2_tab[3] = r2[3] * to_table;
 
  r2_int[0] = int(r2_tab[0]);
  r2_int[1] = int(r2_tab[1]);
  r2_int[2] = int(r2_tab[2]);
  r2_int[3] = int(r2_tab[3]);
  
  DEBUG(12, "\t\t0: r2: " << r2[0] << " r2_tab: " << r2_tab[0]);
  DEBUG(12, "\t\t1: r2: " << r2[1] << " r2_tab: " << r2_tab[1]);
  DEBUG(12, "\t\t2: r2: " << r2[2] << " r2_tab: " << r2_tab[2]);
  DEBUG(12, "\t\t3: r2: " << r2[3] << " r2_tab: " << r2_tab[3]);
  
  assert(r2_int[0] < table_size);
  assert(r2_int[1] < table_size);
  assert(r2_int[2] < table_size);
  assert(r2_int[3] < table_size);

  r2_tab[0] -= r2_int[0];
  r2_tab[1] -= r2_int[1];
  r2_tab[2] -= r2_int[2];
  r2_tab[3] -= r2_int[3];
  
  values[0] = &OH_table[r2_int[0]];
  values[1] = &OH_table[r2_int[1]];
  values[2] = &OH_table[r2_int[2]];
  values[3] = &OH_table[r2_int[3]];

  e_crf += values[0]->e_crf + values[0]->e_crf_diff * r2_tab[0];
  e_crf += values[1]->e_crf + values[1]->e_crf_diff * r2_tab[1];
  e_crf += values[2]->e_crf + values[2]->e_crf_diff * r2_tab[2];
  e_crf += values[3]->e_crf + values[3]->e_crf_diff * r2_tab[3];
  
  DEBUG(10, "e_crf: " << e_crf);  

  ff[0] = values[0]->f + values[0]->f_diff * r2_tab[0];
  ff[1] = values[1]->f + values[1]->f_diff * r2_tab[1];
  ff[2] = values[2]->f + values[2]->f_diff * r2_tab[2];
  ff[3] = values[3]->f + values[3]->f_diff * r2_tab[3];
  
  DEBUG(10, "f: " << ff[0]);
  DEBUG(10, "f: " << ff[1]);
  DEBUG(10, "f: " << ff[2]);
  DEBUG(10, "f: " << ff[3]);

  fx = ff[0] * x[0];
  fy = ff[0] * y[0];
  fz = ff[0] * z[0];
      
  (*force_i)(0) += fx;
  (*(force_j+1))(0) -= fx;
  (*force_i)(1) += fy;
  (*(force_j+1))(1) -= fy;
  (*force_i)(2) += fz;
  (*(force_j+1))(2) -= fz;  
  
  storage.virial_tensor(0, 0) += x[0] * fx;
  storage.virial_tensor(0, 1) += x[0] * fy;
  storage.virial_tensor(0, 2) += x[0] * fz;
  storage.virial_tensor(1, 0) += y[0] * fx;
  storage.virial_tensor(1, 1) += y[0] * fy;
  storage.virial_tensor(1, 2) += y[0] * fz;
  storage.virial_tensor(2, 0) += z[0] * fx;
  storage.virial_tensor(2, 1) += z[0] * fy;
  storage.virial_tensor(2, 2) += z[0] * fz;
      
  fx = ff[1] * x[1];
  fy = ff[1] * y[1];
  fz = ff[1] * z[1];
      
  (*force_i)(0) += fx;
  (*(force_j+2))(0) -= fx;
  (*force_i)(1) += fy;
  (*(force_j+2))(1) -= fy;
  (*force_i)(2) += fz;
  (*(force_j+2))(2) -= fz;
    
  storage.virial_tensor(0, 0) += x[1] * fx;
  storage.virial_tensor(0, 1) += x[1] * fy;
  storage.virial_tensor(0, 2) += x[1] * fz;
  storage.virial_tensor(1, 0) += y[1] * fx;
  storage.virial_tensor(1, 1) += y[1] * fy;
  storage.virial_tensor(1, 2) += y[1] * fz;
  storage.virial_tensor(2, 0) += z[1] * fx;
  storage.virial_tensor(2, 1) += z[1] * fy;
  storage.virial_tensor(2, 2) += z[1] * fz;
      
  fx = ff[2] * x[2];
  fy = ff[2] * y[2];
  fz = ff[2] * z[2];
      
  (*(force_i+1))(0) += fx;
  (*(force_j))(0) -= fx;
  (*(force_i+1))(1) += fy;
  (*(force_j))(1) -= fy;
  (*(force_i+1))(2) += fz;
  (*(force_j))(2) -= fz;
  
  storage.virial_tensor(0, 0) += x[2] * fx;
  storage.virial_tensor(0, 1) += x[2] * fy;
  storage.virial_tensor(0, 2) += x[2] * fz;
  storage.virial_tensor(1, 0) += y[2] * fx;
  storage.virial_tensor(1, 1) += y[2] * fy;
  storage.virial_tensor(1, 2) += y[2] * fz;
  storage.virial_tensor(2, 0) += z[2] * fx;
  storage.virial_tensor(2, 1) += z[2] * fy;
  storage.virial_tensor(2, 2) += z[2] * fz;
      
  fx = ff[3] * x[3];
  fy = ff[3] * y[3];
  fz = ff[3] * z[3];
      
  (*(force_i+2))(0) += fx;
  (*(force_j))(0) -= fx;
  (*(force_i+2))(1) += fy;
  (*(force_j))(1) -= fy;
  (*(force_i+2))(2) += fz;
  (*(force_j))(2) -= fz;

  storage.virial_tensor(0, 0) += x[3] * fx;
  storage.virial_tensor(0, 1) += x[3] * fy;
  storage.virial_tensor(0, 2) += x[3] * fz;
  storage.virial_tensor(1, 0) += y[3] * fx;
  storage.virial_tensor(1, 1) += y[3] * fy;
  storage.virial_tensor(1, 2) += y[3] * fz;
  storage.virial_tensor(2, 0) += z[3] * fx;
  storage.virial_tensor(2, 1) += z[3] * fy;
  storage.virial_tensor(2, 2) += z[3] * fz;
      
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
 
  r2_tab[0] = r2[0] * to_table; 
  r2_tab[1] = r2[1] * to_table; 
  r2_tab[2] = r2[2] * to_table; 
  r2_tab[3] = r2[3] * to_table; 
 
  r2_int[0] = int(r2_tab[0]);
  r2_int[1] = int(r2_tab[1]);
  r2_int[2] = int(r2_tab[2]);
  r2_int[3] = int(r2_tab[3]);
  
  DEBUG(12, "\t\t0: r2: " << r2[0] << " r2_tab: " << r2_tab[0]);
  DEBUG(12, "\t\t1: r2: " << r2[1] << " r2_tab: " << r2_tab[1]);
  DEBUG(12, "\t\t2: r2: " << r2[2] << " r2_tab: " << r2_tab[2]);
  DEBUG(12, "\t\t3: r2: " << r2[3] << " r2_tab: " << r2_tab[3]); 
  
  assert(r2_int[0] < table_size);
  assert(r2_int[1] < table_size);
  assert(r2_int[2] < table_size);
  assert(r2_int[3] < table_size);

  r2_tab[0] -= r2_int[0];
  r2_tab[1] -= r2_int[1];
  r2_tab[2] -= r2_int[2];
  r2_tab[3] -= r2_int[3];
  
  values[0] = &HH_table[r2_int[0]];
  values[1] = &HH_table[r2_int[1]];
  values[2] = &HH_table[r2_int[2]];
  values[3] = &HH_table[r2_int[3]];

  e_crf += values[0]->e_crf + values[0]->e_crf_diff * r2_tab[0];
  e_crf += values[1]->e_crf + values[1]->e_crf_diff * r2_tab[1];
  e_crf += values[2]->e_crf + values[2]->e_crf_diff * r2_tab[2];
  e_crf += values[3]->e_crf + values[3]->e_crf_diff * r2_tab[3];
  
  DEBUG(10, "e_crf: " << e_crf);

  ff[0] = values[0]->f + values[0]->f_diff * r2_tab[0];
  ff[1] = values[1]->f + values[1]->f_diff * r2_tab[1];
  ff[2] = values[2]->f + values[2]->f_diff * r2_tab[2];
  ff[3] = values[3]->f + values[3]->f_diff * r2_tab[3];
  
  DEBUG(10, "f: " << ff[0]);
  DEBUG(10, "f: " << ff[1]);
  DEBUG(10, "f: " << ff[2]);
  DEBUG(10, "f: " << ff[3]);
  
  fx = ff[0] * x[0];
  fy = ff[0] * y[0];
  fz = ff[0] * z[0];
  
  (*(force_i+1))(0) += fx;
  (*(force_j+1))(0) -= fx;
  (*(force_i+1))(1) += fy;
  (*(force_j+1))(1) -= fy;
  (*(force_i+1))(2) += fz;
  (*(force_j+1))(2) -= fz;

  storage.virial_tensor(0, 0) += x[0] * fx;
  storage.virial_tensor(0, 1) += x[0] * fy;
  storage.virial_tensor(0, 2) += x[0] * fz;
  storage.virial_tensor(1, 0) += y[0] * fx;
  storage.virial_tensor(1, 1) += y[0] * fy;
  storage.virial_tensor(1, 2) += y[0] * fz;
  storage.virial_tensor(2, 0) += z[0] * fx;
  storage.virial_tensor(2, 1) += z[0] * fy;
  storage.virial_tensor(2, 2) += z[0] * fz;

  fx = ff[1] * x[1];
  fy = ff[1] * y[1];
  fz = ff[1] * z[1];
  
  (*(force_i+1))(0) += fx;
  (*(force_j+2))(0) -= fx;
  (*(force_i+1))(1) += fy;
  (*(force_j+2))(1) -= fy;
  (*(force_i+1))(2) += fz;
  (*(force_j+2))(2) -= fz;

  storage.virial_tensor(0, 0) += x[1] * fx;
  storage.virial_tensor(0, 1) += x[1] * fy;
  storage.virial_tensor(0, 2) += x[1] * fz;
  storage.virial_tensor(1, 0) += y[1] * fx;
  storage.virial_tensor(1, 1) += y[1] * fy;
  storage.virial_tensor(1, 2) += y[1] * fz;
  storage.virial_tensor(2, 0) += z[1] * fx;
  storage.virial_tensor(2, 1) += z[1] * fy;
  storage.virial_tensor(2, 2) += z[1] * fz;

  fx = ff[2] * x[2];
  fy = ff[2] * y[2];
  fz = ff[2] * z[2];
  
  (*(force_i+2))(0) += fx;
  (*(force_j+1))(0) -= fx;
  (*(force_i+2))(1) += fy;
  (*(force_j+1))(1) -= fy;
  (*(force_i+2))(2) += fz;
  (*(force_j+1))(2) -= fz;

  storage.virial_tensor(0, 0) += x[2] * fx;
  storage.virial_tensor(0, 1) += x[2] * fy;
  storage.virial_tensor(0, 2) += x[2] * fz;
  storage.virial_tensor(1, 0) += y[2] * fx;
  storage.virial_tensor(1, 1) += y[2] * fy;
  storage.virial_tensor(1, 2) += y[2] * fz;
  storage.virial_tensor(2, 0) += z[2] * fx;
  storage.virial_tensor(2, 1) += z[2] * fy;
  storage.virial_tensor(2, 2) += z[2] * fz;

  fx = ff[3] * x[3];
  fy = ff[3] * y[3];
  fz = ff[3] * z[3];
  
  (*(force_i+2))(0) += fx;
  (*(force_j+2))(0) -= fx;
  (*(force_i+2))(1) += fy;
  (*(force_j+2))(1) -= fy;
  (*(force_i+2))(2) += fz;
  (*(force_j+2))(2) -= fz;

  storage.virial_tensor(0, 0) += x[3] * fx;
  storage.virial_tensor(0, 1) += x[3] * fy;
  storage.virial_tensor(0, 2) += x[3] * fz;
  storage.virial_tensor(1, 0) += y[3] * fx;
  storage.virial_tensor(1, 1) += y[3] * fy;
  storage.virial_tensor(1, 2) += y[3] * fz;
  storage.virial_tensor(2, 0) += z[3] * fx;
  storage.virial_tensor(2, 1) += z[3] * fy;
  storage.virial_tensor(2, 2) += z[3] * fz;
    
  storage.energies.crf_energy[egroup][egroup] += e_crf;
}

template<typename t_nonbonded_spec>
inline void 
interaction::Nonbonded_Innerloop<t_nonbonded_spec>::longrange_spc_table_innerloop
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
  unsigned int r2_int[4];
  double x[4], y[4], z[4], r2[4], r2_tab[4], ff[4], tx = 0.0, ty = 0.0, tz = 0.0, fx = 0.0, fy = 0.0, fz = 0.0;

  // only one energy group
  const int egroup = topo.atom_energy_group(topo.num_solute_atoms());
  DEBUG(8, "\ttable longrange");
  DEBUG(8, "\tspc pair\t" << i << "\t" << j << " egroup " << egroup);
  
  double e_lj = 0.0, e_crf = 0.0;
  
  const int ii = i;
  math::Vec const * const pos_i = &conf.current().pos(ii);
  math::Vec * const force_i = &storage.force(ii);

  const int jj = j;
  math::Vec const * const pos_j = &conf.current().pos(jj);
  math::Vec * const force_j = &storage.force(jj);
  DEBUG(9, "ii = " << ii << " jj = " << jj);    

  // O - O
  periodicity.nearest_image(*pos_i, *pos_j, r);
      
  tx = r(0) - (*pos_i)(0) + (*pos_j)(0);
  ty = r(1) - (*pos_i)(1) + (*pos_j)(1);
  tz = r(2) - (*pos_i)(2) + (*pos_j)(2);
      
  assert(abs2(r) != 0);
      
  r2[0] = abs2(r);
  
  struct lj_crf_values {
    double e_lj, e_lj_diff, e_crf, e_crf_diff, f, f_diff;
  };
  
  struct crf_values {
    double e_crf, e_crf_diff, f, f_diff;
  };
  
  #define LONGRANGE "longrange"
  #include "spc_table.h"
  #undef LONGRANGE

  static const double to_table = table_size / data_range;
  
  const crf_values * values[4];
  
  r2_tab[0] = (r2[0] - data_start) * to_table;    
  r2_int[0] = int(r2_tab[0]);
  assert(r2_int[0] < table_size);
  r2_tab[0] -= r2_int[0];
  DEBUG(10, "r: " << sqrt(r2[0]) << " r2: " << r2[0] << " r2_tab: " << r2_tab[0] << " r2_int: " << r2_int[0]);
  const lj_crf_values & val = OO_table[r2_int[0]];
  
  e_lj = val.e_lj + val.e_lj_diff * r2_tab[0];
  e_crf = val.e_crf + val.e_crf_diff * r2_tab[0];
  ff[0] = val.f + val.f_diff * r2_tab[0];
  
  DEBUG(10, "e_lj: " << e_lj << " c_crf: " << e_crf << " f: " << ff[0]);

  fx = ff[0] * r(0);
  fy = ff[0] * r(1);
  fz = ff[0] * r(2);
      
  (*force_i)(0) += fx;
  (*force_j)(0) -= fx;
  (*force_i)(1) += fy;
  (*force_j)(1) -= fy;
  (*force_i)(2) += fz;
  (*force_j)(2) -= fz;
  
  storage.virial_tensor(0, 0) += r(0) * fx;
  storage.virial_tensor(0, 1) += r(0) * fy;
  storage.virial_tensor(0, 2) += r(0) * fz;
  storage.virial_tensor(1, 0) += r(1) * fx;
  storage.virial_tensor(1, 1) += r(1) * fy;
  storage.virial_tensor(1, 2) += r(1) * fz;
  storage.virial_tensor(2, 0) += r(2) * fx;
  storage.virial_tensor(2, 1) += r(2) * fy;
  storage.virial_tensor(2, 2) += r(2) * fz;
      
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

  r2_tab[0] = (r2[0] - data_start) * to_table;
  r2_tab[1] = (r2[1] - data_start) * to_table;
  r2_tab[2] = (r2[2] - data_start) * to_table;
  r2_tab[3] = (r2[3] - data_start) * to_table;
 
  r2_int[0] = int(r2_tab[0]);
  r2_int[1] = int(r2_tab[1]);
  r2_int[2] = int(r2_tab[2]);
  r2_int[3] = int(r2_tab[3]);
  
  DEBUG(12, "\t\t0: r2: " << r2[0] << " r2_tab: " << r2_tab[0]);
  DEBUG(12, "\t\t1: r2: " << r2[1] << " r2_tab: " << r2_tab[1]);
  DEBUG(12, "\t\t2: r2: " << r2[2] << " r2_tab: " << r2_tab[2]);
  DEBUG(12, "\t\t3: r2: " << r2[3] << " r2_tab: " << r2_tab[3]);
  
  assert(r2_int[0] < table_size);
  assert(r2_int[1] < table_size);
  assert(r2_int[2] < table_size);
  assert(r2_int[3] < table_size);

  r2_tab[0] -= r2_int[0];
  r2_tab[1] -= r2_int[1];
  r2_tab[2] -= r2_int[2];
  r2_tab[3] -= r2_int[3];
  
  values[0] = &OH_table[r2_int[0]];
  values[1] = &OH_table[r2_int[1]];
  values[2] = &OH_table[r2_int[2]];
  values[3] = &OH_table[r2_int[3]];

  e_crf += values[0]->e_crf + values[0]->e_crf_diff * r2_tab[0];
  e_crf += values[1]->e_crf + values[1]->e_crf_diff * r2_tab[1];
  e_crf += values[2]->e_crf + values[2]->e_crf_diff * r2_tab[2];
  e_crf += values[3]->e_crf + values[3]->e_crf_diff * r2_tab[3];

  DEBUG(10, "e_crf: " << e_crf);
  
  ff[0] = values[0]->f + values[0]->f_diff * r2_tab[0];
  ff[1] = values[1]->f + values[1]->f_diff * r2_tab[1];
  ff[2] = values[2]->f + values[2]->f_diff * r2_tab[2];
  ff[3] = values[3]->f + values[3]->f_diff * r2_tab[3];
  
  DEBUG(10, "f: " << ff[0]);
  DEBUG(10, "f: " << ff[1]);
  DEBUG(10, "f: " << ff[2]);
  DEBUG(10, "f: " << ff[3]);

  fx = ff[0] * x[0];
  fy = ff[0] * y[0];
  fz = ff[0] * z[0];
      
  (*force_i)(0) += fx;
  (*(force_j+1))(0) -= fx;
  (*force_i)(1) += fy;
  (*(force_j+1))(1) -= fy;
  (*force_i)(2) += fz;
  (*(force_j+1))(2) -= fz;  
  
  storage.virial_tensor(0, 0) += x[0] * fx;
  storage.virial_tensor(0, 1) += x[0] * fy;
  storage.virial_tensor(0, 2) += x[0] * fz;
  storage.virial_tensor(1, 0) += y[0] * fx;
  storage.virial_tensor(1, 1) += y[0] * fy;
  storage.virial_tensor(1, 2) += y[0] * fz;
  storage.virial_tensor(2, 0) += z[0] * fx;
  storage.virial_tensor(2, 1) += z[0] * fy;
  storage.virial_tensor(2, 2) += z[0] * fz;
      
  fx = ff[1] * x[1];
  fy = ff[1] * y[1];
  fz = ff[1] * z[1];
      
  (*force_i)(0) += fx;
  (*(force_j+2))(0) -= fx;
  (*force_i)(1) += fy;
  (*(force_j+2))(1) -= fy;
  (*force_i)(2) += fz;
  (*(force_j+2))(2) -= fz;
    
  storage.virial_tensor(0, 0) += x[1] * fx;
  storage.virial_tensor(0, 1) += x[1] * fy;
  storage.virial_tensor(0, 2) += x[1] * fz;
  storage.virial_tensor(1, 0) += y[1] * fx;
  storage.virial_tensor(1, 1) += y[1] * fy;
  storage.virial_tensor(1, 2) += y[1] * fz;
  storage.virial_tensor(2, 0) += z[1] * fx;
  storage.virial_tensor(2, 1) += z[1] * fy;
  storage.virial_tensor(2, 2) += z[1] * fz;
      
  fx = ff[2] * x[2];
  fy = ff[2] * y[2];
  fz = ff[2] * z[2];
      
  (*(force_i+1))(0) += fx;
  (*(force_j))(0) -= fx;
  (*(force_i+1))(1) += fy;
  (*(force_j))(1) -= fy;
  (*(force_i+1))(2) += fz;
  (*(force_j))(2) -= fz;
  
  storage.virial_tensor(0, 0) += x[2] * fx;
  storage.virial_tensor(0, 1) += x[2] * fy;
  storage.virial_tensor(0, 2) += x[2] * fz;
  storage.virial_tensor(1, 0) += y[2] * fx;
  storage.virial_tensor(1, 1) += y[2] * fy;
  storage.virial_tensor(1, 2) += y[2] * fz;
  storage.virial_tensor(2, 0) += z[2] * fx;
  storage.virial_tensor(2, 1) += z[2] * fy;
  storage.virial_tensor(2, 2) += z[2] * fz;
      
  fx = ff[3] * x[3];
  fy = ff[3] * y[3];
  fz = ff[3] * z[3];
      
  (*(force_i+2))(0) += fx;
  (*(force_j))(0) -= fx;
  (*(force_i+2))(1) += fy;
  (*(force_j))(1) -= fy;
  (*(force_i+2))(2) += fz;
  (*(force_j))(2) -= fz;

  storage.virial_tensor(0, 0) += x[3] * fx;
  storage.virial_tensor(0, 1) += x[3] * fy;
  storage.virial_tensor(0, 2) += x[3] * fz;
  storage.virial_tensor(1, 0) += y[3] * fx;
  storage.virial_tensor(1, 1) += y[3] * fy;
  storage.virial_tensor(1, 2) += y[3] * fz;
  storage.virial_tensor(2, 0) += z[3] * fx;
  storage.virial_tensor(2, 1) += z[3] * fy;
  storage.virial_tensor(2, 2) += z[3] * fz;
      
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
 
  r2_tab[0] = (r2[0] - data_start) * to_table; 
  r2_tab[1] = (r2[1] - data_start) * to_table; 
  r2_tab[2] = (r2[2] - data_start) * to_table; 
  r2_tab[3] = (r2[3] - data_start) * to_table; 
 
  r2_int[0] = int(r2_tab[0]);
  r2_int[1] = int(r2_tab[1]);
  r2_int[2] = int(r2_tab[2]);
  r2_int[3] = int(r2_tab[3]);
  
  DEBUG(12, "\t\t0: r2: " << r2[0] << " r2_tab: " << r2_tab[0]);
  DEBUG(12, "\t\t1: r2: " << r2[1] << " r2_tab: " << r2_tab[1]);
  DEBUG(12, "\t\t2: r2: " << r2[2] << " r2_tab: " << r2_tab[2]);
  DEBUG(12, "\t\t3: r2: " << r2[3] << " r2_tab: " << r2_tab[3]);  
  
  assert(r2_int[0] < table_size);
  assert(r2_int[1] < table_size);
  assert(r2_int[2] < table_size);
  assert(r2_int[3] < table_size);

  r2_tab[0] -= r2_int[0];
  r2_tab[1] -= r2_int[1];
  r2_tab[2] -= r2_int[2];
  r2_tab[3] -= r2_int[3];
  
  values[0] = &HH_table[r2_int[0]];
  values[1] = &HH_table[r2_int[1]];
  values[2] = &HH_table[r2_int[2]];
  values[3] = &HH_table[r2_int[3]];

  e_crf += values[0]->e_crf + values[0]->e_crf_diff * r2_tab[0];
  e_crf += values[1]->e_crf + values[1]->e_crf_diff * r2_tab[1];
  e_crf += values[2]->e_crf + values[2]->e_crf_diff * r2_tab[2];
  e_crf += values[3]->e_crf + values[3]->e_crf_diff * r2_tab[3];

  DEBUG(10, "e_crf: " << e_crf);         
  
  ff[0] = values[0]->f + values[0]->f_diff * r2_tab[0];
  ff[1] = values[1]->f + values[1]->f_diff * r2_tab[1];
  ff[2] = values[2]->f + values[2]->f_diff * r2_tab[2];
  ff[3] = values[3]->f + values[3]->f_diff * r2_tab[3];
  
  DEBUG(10, "f: " << ff[0]);
  DEBUG(10, "f: " << ff[1]);
  DEBUG(10, "f: " << ff[2]);
  DEBUG(10, "f: " << ff[3]);
  
  fx = ff[0] * x[0];
  fy = ff[0] * y[0];
  fz = ff[0] * z[0];
  
  (*(force_i+1))(0) += fx;
  (*(force_j+1))(0) -= fx;
  (*(force_i+1))(1) += fy;
  (*(force_j+1))(1) -= fy;
  (*(force_i+1))(2) += fz;
  (*(force_j+1))(2) -= fz;

  storage.virial_tensor(0, 0) += x[0] * fx;
  storage.virial_tensor(0, 1) += x[0] * fy;
  storage.virial_tensor(0, 2) += x[0] * fz;
  storage.virial_tensor(1, 0) += y[0] * fx;
  storage.virial_tensor(1, 1) += y[0] * fy;
  storage.virial_tensor(1, 2) += y[0] * fz;
  storage.virial_tensor(2, 0) += z[0] * fx;
  storage.virial_tensor(2, 1) += z[0] * fy;
  storage.virial_tensor(2, 2) += z[0] * fz;

  fx = ff[1] * x[1];
  fy = ff[1] * y[1];
  fz = ff[1] * z[1];
  
  (*(force_i+1))(0) += fx;
  (*(force_j+2))(0) -= fx;
  (*(force_i+1))(1) += fy;
  (*(force_j+2))(1) -= fy;
  (*(force_i+1))(2) += fz;
  (*(force_j+2))(2) -= fz;

  storage.virial_tensor(0, 0) += x[1] * fx;
  storage.virial_tensor(0, 1) += x[1] * fy;
  storage.virial_tensor(0, 2) += x[1] * fz;
  storage.virial_tensor(1, 0) += y[1] * fx;
  storage.virial_tensor(1, 1) += y[1] * fy;
  storage.virial_tensor(1, 2) += y[1] * fz;
  storage.virial_tensor(2, 0) += z[1] * fx;
  storage.virial_tensor(2, 1) += z[1] * fy;
  storage.virial_tensor(2, 2) += z[1] * fz;

  fx = ff[2] * x[2];
  fy = ff[2] * y[2];
  fz = ff[2] * z[2];
  
  (*(force_i+2))(0) += fx;
  (*(force_j+1))(0) -= fx;
  (*(force_i+2))(1) += fy;
  (*(force_j+1))(1) -= fy;
  (*(force_i+2))(2) += fz;
  (*(force_j+1))(2) -= fz;

  storage.virial_tensor(0, 0) += x[2] * fx;
  storage.virial_tensor(0, 1) += x[2] * fy;
  storage.virial_tensor(0, 2) += x[2] * fz;
  storage.virial_tensor(1, 0) += y[2] * fx;
  storage.virial_tensor(1, 1) += y[2] * fy;
  storage.virial_tensor(1, 2) += y[2] * fz;
  storage.virial_tensor(2, 0) += z[2] * fx;
  storage.virial_tensor(2, 1) += z[2] * fy;
  storage.virial_tensor(2, 2) += z[2] * fz;

  fx = ff[3] * x[3];
  fy = ff[3] * y[3];
  fz = ff[3] * z[3];
  
  (*(force_i+2))(0) += fx;
  (*(force_j+2))(0) -= fx;
  (*(force_i+2))(1) += fy;
  (*(force_j+2))(1) -= fy;
  (*(force_i+2))(2) += fz;
  (*(force_j+2))(2) -= fz;

  storage.virial_tensor(0, 0) += x[3] * fx;
  storage.virial_tensor(0, 1) += x[3] * fy;
  storage.virial_tensor(0, 2) += x[3] * fz;
  storage.virial_tensor(1, 0) += y[3] * fx;
  storage.virial_tensor(1, 1) += y[3] * fy;
  storage.virial_tensor(1, 2) += y[3] * fz;
  storage.virial_tensor(2, 0) += z[3] * fx;
  storage.virial_tensor(2, 1) += z[3] * fy;
  storage.virial_tensor(2, 2) += z[3] * fz;
    
  storage.energies.crf_energy[egroup][egroup] += e_crf;
}


template<typename t_nonbonded_spec>
inline void 
interaction::Nonbonded_Innerloop<t_nonbonded_spec>::shortrange_spc_table_innerloop
(
  double &e_lj,
  double &e_crf,
  double f[9],
  double r2[9]
)
{
  unsigned int r2_int[4];
  double r2_tab[4];

  DEBUG(8, "\ttable shortrange");

  // O - O
  struct lj_crf_values {
    double e_lj, e_lj_diff, e_crf, e_crf_diff, f, f_diff;
  };
  
  struct crf_values {
    double e_crf, e_crf_diff, f, f_diff;
  };
  
  #define SHORTRANGE "shortrange"
  #include "spc_table.h"
  #undef SHORTRANGE

  static const double to_table = table_size / data_range;
  
  const crf_values * values[4];
  
  // calculate the grid point
  r2_tab[0] = r2[0] * to_table;    
  // round the grid point to the lower integer
  r2_int[0] = int(r2_tab[0]);
  assert(r2_int[0] < table_size);
  // safe the offset from the grid point for linear interpolation
  r2_tab[0] -= r2_int[0];
  DEBUG(10, "r: " << sqrt(r2[0]) << " r2: " << r2[0] << " r2_tab: " << r2_tab[0] << " r2_int: " << r2_int[0]);
  
  // the the data from the grid point
  const lj_crf_values & val = OO_table[r2_int[0]];
  
  // the _diff variables hold the gradient for the interpolation.
  e_lj  += val.e_lj + val.e_lj_diff * r2_tab[0];
  e_crf += val.e_crf + val.e_crf_diff * r2_tab[0];
  f[0] = val.f + val.f_diff * r2_tab[0];
  
  DEBUG(10, "e_lj: " << e_lj << " c_crf: " << e_crf << " f: " << f[0]);
      
  // O - H interactions...
  r2_tab[0] = r2[1] * to_table;
  r2_tab[1] = r2[2] * to_table;
  r2_tab[2] = r2[3] * to_table;
  r2_tab[3] = r2[4] * to_table;
 
  r2_int[0] = int(r2_tab[0]);
  r2_int[1] = int(r2_tab[1]);
  r2_int[2] = int(r2_tab[2]);
  r2_int[3] = int(r2_tab[3]);
  
  DEBUG(12, "\t\t0: r2: " << r2[1] << " r2_tab: " << r2_tab[0]);
  DEBUG(12, "\t\t1: r2: " << r2[2] << " r2_tab: " << r2_tab[1]);
  DEBUG(12, "\t\t2: r2: " << r2[3] << " r2_tab: " << r2_tab[2]);
  DEBUG(12, "\t\t3: r2: " << r2[4] << " r2_tab: " << r2_tab[3]);
  
  assert(r2_int[0] < table_size);
  assert(r2_int[1] < table_size);
  assert(r2_int[2] < table_size);
  assert(r2_int[3] < table_size);

  r2_tab[0] -= r2_int[0];
  r2_tab[1] -= r2_int[1];
  r2_tab[2] -= r2_int[2];
  r2_tab[3] -= r2_int[3];
  
  values[0] = &OH_table[r2_int[0]];
  values[1] = &OH_table[r2_int[1]];
  values[2] = &OH_table[r2_int[2]];
  values[3] = &OH_table[r2_int[3]];

  e_crf += values[0]->e_crf + values[0]->e_crf_diff * r2_tab[0];
  e_crf += values[1]->e_crf + values[1]->e_crf_diff * r2_tab[1];
  e_crf += values[2]->e_crf + values[2]->e_crf_diff * r2_tab[2];
  e_crf += values[3]->e_crf + values[3]->e_crf_diff * r2_tab[3];

  DEBUG(10, "e_crf: " << e_crf);

  f[1] = values[0]->f + values[0]->f_diff * r2_tab[0];
  f[2] = values[1]->f + values[1]->f_diff * r2_tab[1];
  f[3] = values[2]->f + values[2]->f_diff * r2_tab[2];
  f[4] = values[3]->f + values[3]->f_diff * r2_tab[3];
  
  DEBUG(10, "f: " << f[1]);
  DEBUG(10, "f: " << f[2]);
  DEBUG(10, "f: " << f[3]);
  DEBUG(10, "f: " << f[4]);
      
  // H - H interactions...
 
  r2_tab[0] = r2[5] * to_table; 
  r2_tab[1] = r2[6] * to_table; 
  r2_tab[2] = r2[7] * to_table; 
  r2_tab[3] = r2[8] * to_table; 
 
  r2_int[0] = int(r2_tab[0]);
  r2_int[1] = int(r2_tab[1]);
  r2_int[2] = int(r2_tab[2]);
  r2_int[3] = int(r2_tab[3]);
  
  DEBUG(12, "\t\t0: r2: " << r2[5] << " r2_tab: " << r2_tab[0]);
  DEBUG(12, "\t\t1: r2: " << r2[6] << " r2_tab: " << r2_tab[1]);
  DEBUG(12, "\t\t2: r2: " << r2[7] << " r2_tab: " << r2_tab[2]);
  DEBUG(12, "\t\t3: r2: " << r2[8] << " r2_tab: " << r2_tab[3]);
  
  assert(r2_int[0] < table_size);
  assert(r2_int[1] < table_size);
  assert(r2_int[2] < table_size);
  assert(r2_int[3] < table_size);

  r2_tab[0] -= r2_int[0];
  r2_tab[1] -= r2_int[1];
  r2_tab[2] -= r2_int[2];
  r2_tab[3] -= r2_int[3];
  
  values[0] = &HH_table[r2_int[0]];
  values[1] = &HH_table[r2_int[1]];
  values[2] = &HH_table[r2_int[2]];
  values[3] = &HH_table[r2_int[3]];

  e_crf += values[0]->e_crf + values[0]->e_crf_diff * r2_tab[0];
  e_crf += values[1]->e_crf + values[1]->e_crf_diff * r2_tab[1];
  e_crf += values[2]->e_crf + values[2]->e_crf_diff * r2_tab[2];
  e_crf += values[3]->e_crf + values[3]->e_crf_diff * r2_tab[3];

  DEBUG(10, "e_crf: " << e_crf);

  f[5] = values[0]->f + values[0]->f_diff * r2_tab[0];
  f[6] = values[1]->f + values[1]->f_diff * r2_tab[1];
  f[7] = values[2]->f + values[2]->f_diff * r2_tab[2];
  f[8] = values[3]->f + values[3]->f_diff * r2_tab[3];
  
  DEBUG(10, "f: " << f[5]);
  DEBUG(10, "f: " << f[6]);
  DEBUG(10, "f: " << f[7]);
  DEBUG(10, "f: " << f[8]);
}

template<typename t_nonbonded_spec>
inline void 
interaction::Nonbonded_Innerloop<t_nonbonded_spec>::longrange_spc_table_innerloop
(
  double &e_lj,
  double &e_crf,
  double f[9],
  double r2[9]
)
{
  unsigned int r2_int[4];
  double r2_tab[4];

  DEBUG(8, "\ttable longrange");

  // O - O
  struct lj_crf_values {
    double e_lj, e_lj_diff, e_crf, e_crf_diff, f, f_diff;
  };
  
  struct crf_values {
    double e_crf, e_crf_diff, f, f_diff;
  };
  
  #define LONGRANGE "longrange"
  #include "spc_table.h"
  #undef LONGRANGE

  static const double to_table = table_size / data_range;
  
  const crf_values * values[4];
  
  // calculate the grid point
  r2_tab[0] = (r2[0] - data_start) * to_table;    
  // round the grid point to the lower integer
  r2_int[0] = int(r2_tab[0]);
  assert(r2_int[0] < table_size);
  // safe the offset from the grid point for linear interpolation
  r2_tab[0] -= r2_int[0];
  DEBUG(10, "r: " << sqrt(r2[0]) << " r2: " << r2[0] << " r2_tab: " << r2_tab[0] << " r2_int: " << r2_int[0]);
  
  // the the data from the grid point
  const lj_crf_values & val = OO_table[r2_int[0]];
  
  // the _diff variables hold the gradient for the interpolation.
  e_lj  += val.e_lj + val.e_lj_diff * r2_tab[0];
  e_crf += val.e_crf + val.e_crf_diff * r2_tab[0];
  f[0] = val.f + val.f_diff * r2_tab[0];
  
  DEBUG(10, "e_lj: " << e_lj << " c_crf: " << e_crf << " f: " << f[0]);
      
  // O - H interactions...
  r2_tab[0] = (r2[1] - data_start) * to_table;
  r2_tab[1] = (r2[2] - data_start) * to_table;
  r2_tab[2] = (r2[3] - data_start) * to_table;
  r2_tab[3] = (r2[4] - data_start) * to_table;
 
  r2_int[0] = int(r2_tab[0]);
  r2_int[1] = int(r2_tab[1]);
  r2_int[2] = int(r2_tab[2]);
  r2_int[3] = int(r2_tab[3]);
  
  DEBUG(12, "\t\t0: r2: " << r2[1] << " r2_tab: " << r2_tab[0]);
  DEBUG(12, "\t\t1: r2: " << r2[2] << " r2_tab: " << r2_tab[1]);
  DEBUG(12, "\t\t2: r2: " << r2[3] << " r2_tab: " << r2_tab[2]);
  DEBUG(12, "\t\t3: r2: " << r2[4] << " r2_tab: " << r2_tab[3]);
  
  assert(r2_int[0] < table_size);
  assert(r2_int[1] < table_size);
  assert(r2_int[2] < table_size);
  assert(r2_int[3] < table_size);

  r2_tab[0] -= r2_int[0];
  r2_tab[1] -= r2_int[1];
  r2_tab[2] -= r2_int[2];
  r2_tab[3] -= r2_int[3];
  
  values[0] = &OH_table[r2_int[0]];
  values[1] = &OH_table[r2_int[1]];
  values[2] = &OH_table[r2_int[2]];
  values[3] = &OH_table[r2_int[3]];

  e_crf += values[0]->e_crf + values[0]->e_crf_diff * r2_tab[0];
  e_crf += values[1]->e_crf + values[1]->e_crf_diff * r2_tab[1];
  e_crf += values[2]->e_crf + values[2]->e_crf_diff * r2_tab[2];
  e_crf += values[3]->e_crf + values[3]->e_crf_diff * r2_tab[3];

  DEBUG(10, "e_crf: " << e_crf);

  f[1] = values[0]->f + values[0]->f_diff * r2_tab[0];
  f[2] = values[1]->f + values[1]->f_diff * r2_tab[1];
  f[3] = values[2]->f + values[2]->f_diff * r2_tab[2];
  f[4] = values[3]->f + values[3]->f_diff * r2_tab[3];
  
  DEBUG(10, "f: " << f[1]);
  DEBUG(10, "f: " << f[2]);
  DEBUG(10, "f: " << f[3]);
  DEBUG(10, "f: " << f[4]);
      
  // H - H interactions...
 
  r2_tab[0] = (r2[5] - data_start) * to_table; 
  r2_tab[1] = (r2[6] - data_start) * to_table; 
  r2_tab[2] = (r2[7] - data_start) * to_table; 
  r2_tab[3] = (r2[8] - data_start) * to_table; 
 
  r2_int[0] = int(r2_tab[0]);
  r2_int[1] = int(r2_tab[1]);
  r2_int[2] = int(r2_tab[2]);
  r2_int[3] = int(r2_tab[3]);
  
  DEBUG(12, "\t\t0: r2: " << r2[5] << " r2_tab: " << r2_tab[0]);
  DEBUG(12, "\t\t1: r2: " << r2[6] << " r2_tab: " << r2_tab[1]);
  DEBUG(12, "\t\t2: r2: " << r2[7] << " r2_tab: " << r2_tab[2]);
  DEBUG(12, "\t\t3: r2: " << r2[8] << " r2_tab: " << r2_tab[3]);
  
  assert(r2_int[0] < table_size);
  assert(r2_int[1] < table_size);
  assert(r2_int[2] < table_size);
  assert(r2_int[3] < table_size);

  r2_tab[0] -= r2_int[0];
  r2_tab[1] -= r2_int[1];
  r2_tab[2] -= r2_int[2];
  r2_tab[3] -= r2_int[3];
  
  values[0] = &HH_table[r2_int[0]];
  values[1] = &HH_table[r2_int[1]];
  values[2] = &HH_table[r2_int[2]];
  values[3] = &HH_table[r2_int[3]];

  e_crf += values[0]->e_crf + values[0]->e_crf_diff * r2_tab[0];
  e_crf += values[1]->e_crf + values[1]->e_crf_diff * r2_tab[1];
  e_crf += values[2]->e_crf + values[2]->e_crf_diff * r2_tab[2];
  e_crf += values[3]->e_crf + values[3]->e_crf_diff * r2_tab[3];

  DEBUG(10, "e_crf: " << e_crf);

  f[5] = values[0]->f + values[0]->f_diff * r2_tab[0];
  f[6] = values[1]->f + values[1]->f_diff * r2_tab[1];
  f[7] = values[2]->f + values[2]->f_diff * r2_tab[2];
  f[8] = values[3]->f + values[3]->f_diff * r2_tab[3];
  
  DEBUG(10, "f: " << f[5]);
  DEBUG(10, "f: " << f[6]);
  DEBUG(10, "f: " << f[7]);
  DEBUG(10, "f: " << f[8]);
}
