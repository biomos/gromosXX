/**
 * @file perturbed_dihedral_constraint.cc
 * contains the perturbed dihedral constraint iteration
 * as its a template, the file is included in perturbed_shake.cc
 */

/**
 * Perturbed Dihedral Constraints
 * 
 * see: Sampling of rare events using hidden restraints
 * in preparation for: Journal of Physical Chemistry
 * Markus Christen, Anna-Pitschna E. Kunz and
 * Wilfred F. van Gunsteren
 * 2006
 * Appendix
 *
 */
template<math::boundary_enum B, math::virial_enum V>
int algorithm::Perturbed_Shake::perturbed_dih_constr_iteration
(
 topology::Topology const & topo,
 configuration::Configuration & conf,
 simulation::Simulation const & sim,
 bool & convergence,
 std::vector<bool> &skip_now,
 std::vector<bool> &skip_next,
 math::Periodicity<B> const & periodicity
 )
{
  const double tolerance = 1.0 / 180.0 * math::Pi;
  // const double delta = 45.0 / 180.0 * math::Pi;
  
  convergence = true;

  math::VArray &pos   = conf.current().pos;
  math::VArray &ref   = conf.old().pos;

  std::vector<topology::perturbed_dihedral_restraint_struct>::const_iterator 
    it = topo.perturbed_dihedral_restraints().begin(),
    to = topo.perturbed_dihedral_restraints().end();

  for( ; it != to; ++it){
	
    // check whether we can skip this constraint
    if (skip_now[it->i] && skip_now[it->j] &&
	skip_now[it->k] && skip_now[it->l]) continue;

    // atom i determines the energy group for the output. 
    // we use the same definition for the individual lambdas
    // we use the dihedral restaints lambda
    const double lambda = topo.individual_lambda(simulation::dihres_lambda)
      [topo.atom_energy_group()[it->i]][topo.atom_energy_group()[it->i]];
    const double lambda_derivative = topo.individual_lambda_derivative
      (simulation::dihres_lambda)
      [topo.atom_energy_group()[it->i]][topo.atom_energy_group()[it->i]];

    // first equation in Dihedral-angle constraints (Appendix)
    // is number 36.
    
    // calculate dihedral angle
    DEBUG(9, "dihedral angle " 
	  << it->i << "-" << it->j << "-" << it->k << "-" << it->l);

    math::Vec r12, r32, r34;
    periodicity.nearest_image(pos(it->i), pos(it->j), r12);
    periodicity.nearest_image(pos(it->k), pos(it->j), r32);
    periodicity.nearest_image(pos(it->k), pos(it->l), r34);

    // eq 37
    const math::Vec r52 = math::cross(r12, r32);
    const double d52 = math::abs(r52);
    // eq 38
    const math::Vec r63 = math::cross(r32, r34);
    const double d63 = math::abs(r63);
    // eq 39
    const int sign_phi = (math::dot(r12, r63) >= 0.0) ? 1 : -1;
    // eq 36
    const double cos_phi = math::dot(r52, r63) / (d52 * d63);
    double phi = sign_phi * acos(cos_phi);

    const double phi0 = (1.0 - lambda) * it->A_phi + lambda * it->B_phi;

    while(phi < phi0 - math::Pi)
      phi += 2 * math::Pi;
    while(phi > phi0 + math::Pi)
      phi -= 2 * math::Pi;
    
    // decide if constraint is fulfilled using phi
    // and not cos(phi)
    // and a relatvie criterion is not what we want here

    const double diff = phi - phi0;
    DEBUG(8, "phi=" << 180 * phi / math::Pi 
	  << "\tphi0=" << 180 * phi0 / math::Pi);
    
    if(fabs(diff) >= tolerance){
      // we have to shake
      DEBUG(10, "shaking");
      
      math::Vec ref12, ref32, ref34;
      periodicity.nearest_image(ref(it->i), ref(it->j), ref12);
      periodicity.nearest_image(ref(it->k), ref(it->j), ref32);
      periodicity.nearest_image(ref(it->k), ref(it->l), ref34);
      
      // eq 37
      const math::Vec ref52 = math::cross(ref12, ref32);
      const double dref52 = math::abs(ref52);
      // eq 38
      const math::Vec ref63 = math::cross(ref32, ref34);
      const double dref63 = math::abs(ref63);

      // eq 45
      const double dref32 = math::abs(ref32);

      const math::Vec a1 = 
	dref32/(dref52 * dref52) * ref52;

      const math::Vec a2 = 
	(math::dot(ref12, ref32) / (dref32 * dref32) - 1) *
	dref32 / (dref52 * dref52) * ref52 +
	math::dot(ref34, ref32) / (dref32 * dref32) *
	dref32 / (dref63 * dref63) * ref63;

      const math::Vec a3 =
	- ((math::dot(ref34, ref32) / (dref32 * dref32) - 1) *
	   dref32 / (dref63 * dref63) * ref63 +
	   math::dot(ref12, ref32) / (dref32 * dref32) *
	   dref32 / (dref52 * dref52) * ref52);
      
      const math::Vec a4 =
	- dref32 / (dref63 * dref63) * ref63;

      const double m1 = topo.mass(it->i);
      const double m2 = topo.mass(it->j);
      const double m3 = topo.mass(it->k);
      const double m4 = topo.mass(it->l);

      // eq 50, 51      
      const math::Vec b123 =
	math::cross(r12, a3 / m3 - a2 / m2) -
	math::cross(r32, a1 / m1 - a2 / m2);
      
      // eq 52, 53
      const math::Vec b234 =
	math::cross(r32, a3 / m3 - a4 / m4) -
	math::cross(r34, a3 / m3 - a2 / m2);
      
      // eq 54, 55
      const double c1234 =
	math::dot(
		  math::cross(r12, r32),
		  math::cross(r32, r34)
		  );
      
      const double d1234 =
	math::dot(
		  math::cross(r12, r32),
		  b234
		  ) +
	math::dot(
		  math::cross(r32, r34),
		  b123
		  );

      // Finally, eq 60
      const double dt = sim.time_step_size();
      const double dt2 = dt * dt;
      
      //////////////////////////////////////////////////
      // reference phi!
      const double sin_2phi = sin(2 * phi0);
      const double cos_phi2 = cos(phi0) * cos(phi0);
      
      const double nominator = 
	cos_phi2 *
	math::dot(math::cross(r12, r32), math::cross(r12, r32)) *
	math::dot(math::cross(r32, r34), math::cross(r32, r34)) -
	c1234 * c1234;

      if (cos_phi2 < math::epsilon){
	/*
	io::messages.add("nominator is zero (phi0)",
			 "Perturbed_Dihedral_Constraints",
			 io::message::error);
	*/
	std::cerr << "cos^2(phi0) < math::epsilon: " << cos_phi2 << std::endl;
      }
      if (fabs(sin_2phi) < math::epsilon){
	/*
	io::messages.add("sin(2 phi0) is zero",
			 "Perturbed_Dihedral_Constraints",
			 io::message::error);
	*/
	std::cerr << "sin(2 phi0) < math::epsilon: " << sin_2phi << std::endl;
      }
      
      const double denominator =
	2 * sin_2phi * dt2 *
	(
	 c1234 * d1234 - cos_phi2 *
	 (
	  math::dot(math::cross(r12, r32), math::cross(r12, r32)) *
	  math::dot(math::cross(r32, r34), b234) +
	  math::dot(math::cross(r32, r34), math::cross(r32, r34)) *
	  math::dot(math::cross(r12, r32), b123)
	  )
	 );

      if (fabs(denominator) < math::epsilon * 0.1){
	/*
	io::messages.add("denominator approaches zero!",
			 "Perturbed_Dihedral_Constraints",
			 io::message::error);
	*/
	std::cerr << "denominator < math::epsilon: " << denominator << std::endl;
      }

      const double l = nominator / denominator;
      //////////////////////////////////////////////////

      pos(it->i) += l * sin_2phi * dt2 * a1 / m1;
      pos(it->j) += l * sin_2phi * dt2 * a2 / m2;
      pos(it->k) += l * sin_2phi * dt2 * a3 / m3;
      pos(it->l) += l * sin_2phi * dt2 * a4 / m4;
      
      if (V == math::atomic_virial){
	io::messages.add("atomic virial not implemented. copy from dihedral angle interaction!",
			 "Dihedral Constraints",
			 io::message::error);
      }

      // dH/dl
      // eq 61
      const double dHdl = lambda_derivative * l * sin_2phi * 
	(it->B_phi - it->A_phi);
      conf.old().perturbed_energy_derivatives.
	constraints_energy[topo.atom_energy_group()[it->i]] += dHdl;
      
      convergence = false;

      // consider atoms in the next step
      skip_next[it->i] = false;
      skip_next[it->j] = false;
      skip_next[it->k] = false;
      skip_next[it->l] = false;
      
    } // we have to shake
  } // constraints
      
  return 0;
}    
