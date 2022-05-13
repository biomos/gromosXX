/**
 * @file perturbed_angle_constraint.cc
 * contains the perturbed angle constraint iteration
 * as it's a template, the file is included in shake.cc
 */

/**
 * Perturbed Angle Constraints
 * see: J Comput Chem. 2021;42:418–434.
 *
 */
#ifndef PERT_ANGLE_CONSTRAINT_CC
#define PERT_ANGLE_CONSTRAINT_CC

template <math::boundary_enum B, math::virial_enum V>
int algorithm::Perturbed_Shake::perturbed_ang_constr_iteration(
    topology::Topology const &topo,
    configuration::Configuration &conf,
    simulation::Simulation const &sim,
    bool &convergence,
    std::vector<bool> &skip_now,
    std::vector<bool> &skip_next,
    math::Periodicity<B> const &periodicity)
{
  const double tolerance = sim.param().angrest.tolerance;
  convergence = true;

  math::VArray &pos = conf.current().pos;
  math::VArray &ref = conf.old().pos;

  std::vector<topology::perturbed_angle_restraint_struct>::const_iterator
      it = topo.perturbed_angle_restraints().begin(),
      to = topo.perturbed_angle_restraints().end();

  for (; it != to; ++it)
  {

    // check whether we can skip this constraint
    if (skip_now[it->i] && skip_now[it->j] && skip_now[it->k])
      continue;

    // atom i determines the energy group for the output.
    // we use the same definition for the individual lambdas
    // we use the angle lambda
    const double lambda = topo.individual_lambda(simulation::angle_lambda)
                              [topo.atom_energy_group()[it->i]][topo.atom_energy_group()[it->i]];
    const double lambda_derivative = topo.individual_lambda_derivative(simulation::angle_lambda)
                                         [topo.atom_energy_group()[it->i]][topo.atom_energy_group()[it->i]];

    // calculate angle
    math::Vec r12, r32;
    periodicity.nearest_image(pos(it->i), pos(it->j), r12);
    periodicity.nearest_image(pos(it->k), pos(it->j), r32);

    double d12 = math::abs(r12);
    double d32 = math::abs(r32);

    const long double c1 = math::dot(r12, r32);
    const long double cost = c1 / (d12 * d32);

    double theta = 0.0;
    // cost can be >1 or <-1 because of precision limits
    if (cost > 1)
      theta = 0.0;
    else if (cost < -1)
      theta = math::Pi;
    else
      theta = acos(cost);

    const double theta0 = (1.0 - lambda) * it->A_theta + lambda * it->B_theta;
    // decide if constraint is fulfilled using theta
    // and not cos(theta)
    const double diff = theta - theta0;
    DEBUG(12, "angle constraint" << it->i << "-" << it->j << "-" << it->k << "theta="
                                 << std::setprecision(7) << 180 * theta / math::Pi
                                 << "\ttheta0=" << 180 * theta0 / math::Pi
                                 << std::setprecision(7) << "\tdiff=" << 180 * diff / math::Pi);

    if (fabs(diff) >= tolerance)
    {
      // we have to shake
      DEBUG(10, "shaking angle");

      math::Vec ref12, ref32;
      periodicity.nearest_image(ref(it->i), ref(it->j), ref12);
      periodicity.nearest_image(ref(it->k), ref(it->j), ref32);

      const double dref12 = math::abs(ref12);
      const double dref32 = math::abs(ref32);

      // eq. 18
      math::Vec a123 = (dref12 * dref12 * ref32 - math::dot(ref12, ref32) * ref12) / (dref12 * dref12 * dref12 * dref32);
      math::Vec a321 = (dref32 * dref32 * ref12 - math::dot(ref12, ref32) * ref32) / (dref12 * dref32 * dref32 * dref32);
      //math::Vec a123 = (ref32 - math::dot(ref12, ref32) * ref12 / (dref12 * dref12)) / (dref12 * dref32);
      //math::Vec a321 = (ref12 - math::dot(ref12, ref32) * ref32 / (dref32 * dref32)) / (dref12 * dref32);

      const double m1 = topo.mass(it->i);
      const double m2 = topo.mass(it->j);
      const double m3 = topo.mass(it->k);

      // eq. 28
      math::Vec b123 = a123 / m1 + (a123 + a321) / m2;

      math::Vec b321 = a321 / m3 + (a123 + a321) / m2;

      // eq. 39
      const long double c2 = math::dot(r12, b321) + math::dot(r32, b123);
      // eq. 40
      const double c3 = d12 * d32;
      // eq. 41
      const long double c4 = d12 / d32 * math::dot(r32, b321) + d32 / d12 * math::dot(r12, b123);

      // eq. 43
      double numerator = c1 - c3 * cos(theta0);
      double denom = (c2 - c4 * cos(theta0));

      double l_dt2 = numerator / denom;

      // eq. 14+17
      pos(it->i) -= l_dt2 * a123 / m1;
      pos(it->j) += l_dt2 * (a123 + a321) / m2;
      pos(it->k) -= l_dt2 * a321 / m3;

      /** ------------------------------------*/

      if (V == math::atomic_virial)
      {
        io::messages.add("atomic virial not implemented. copy from angle interaction!",
                         "Angle Constraints", io::message::error);
      }

      const double dt = sim.time_step_size();
      double l = l_dt2 / (dt * dt);

      // dH/dl, eq. 48
      const double dHdl = lambda_derivative * l * sin(theta0) * (it->B_theta - it->A_theta);
      conf.old().perturbed_energy_derivatives.constraints_energy[topo.atom_energy_group()[it->i]] += dHdl;

      convergence = false;

      // consider atoms in the next step
      skip_next[it->i] = false;
      skip_next[it->j] = false;
      skip_next[it->k] = false;

    } // we have to shake
  }   // constraints

  return 0;
}
#endif //PERTURBED_ANGLE_CONSTRAINT_CC
