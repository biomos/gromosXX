/**
 * @file perturbed_dihedral_restraint_interaction.cc
 * methods of Perturbed_Dihedral_Restraint_Interaction
 */

#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"
#include "../../interaction/interaction.h"

#include "../../math/periodicity.h"

// special interactions
#include "../../interaction/interaction_types.h"

#include "../../interaction/special/perturbed_dihedral_restraint_interaction.h"

#include "../../util/template_split.h"
#include "../../util/debug.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE special

/**
 * calculate dihedral restraint interactions
 */
template <math::boundary_enum B, math::virial_enum V>
static int _calculate_perturbed_dihedral_restraint_interactions(
    topology::Topology &topo,
    configuration::Configuration &conf,
    simulation::Simulation &sim)
{
  DEBUG(5, "perturbed dihedral angle restraint interaction");
  // loop over the dihedral restraints
  std::vector<topology::perturbed_dihedral_restraint_struct>::const_iterator
      it = topo.perturbed_dihedral_restraints().begin(),
      to = topo.perturbed_dihedral_restraints().end();

  std::vector<double>::iterator ene_it = conf.special().pertdihedralres.energy.begin();
  std::vector<double>::iterator d_it = conf.special().pertdihedralres.d.begin();

  math::VArray &pos = conf.current().pos;
  math::VArray &force = conf.current().force;
  math::Vec rij, rkj, rkl, rlj, rmj, rnk, fi, fj, fk, fl;
  double dkj2 = 0.0, dkj = 0.0, dmj2 = 0.0, dmj = 0.0, dnk2 = 0.0, dnk = 0.0, ip = 0.0, phi = 0.0;
  double lower_bound = 0.0, upper_bound = 0.0;
  double energy = 0.0, en_term = 0.0, f = 0.0, energy_derivative = 0.0, dlam_term = 0.0;

  math::Periodicity<B> periodicity(conf.current().box);

  for (; it != to; ++it, ++ene_it, ++d_it)
  {

    // atom i determines the energy group for the output.
    // we use the same definition for the individual lambdas
    const double l = topo.individual_lambda(simulation::dihres_lambda)
                         [topo.atom_energy_group()[it->i]]
                         [topo.atom_energy_group()[it->i]];
    const double l_deriv = topo.individual_lambda_derivative(simulation::dihres_lambda)
                               [topo.atom_energy_group()[it->i]]
                               [topo.atom_energy_group()[it->i]];

    periodicity.nearest_image(pos(it->k), pos(it->j), rkj);
    periodicity.nearest_image(pos(it->i), pos(it->j), rij);
    periodicity.nearest_image(pos(it->k), pos(it->l), rkl);

    rmj = cross(rij, rkj);
    rnk = cross(rkj, rkl);

    dkj2 = abs2(rkj);
    dmj2 = abs2(rmj);
    dnk2 = abs2(rnk);
    dkj = sqrt(dkj2);
    dmj = sqrt(dmj2);
    dnk = sqrt(dnk2);

    DEBUG(15, "dkj=" << dkj << " dmj=" << dmj << " dnk=" << dnk);

    assert(dmj != 0.0);
    assert(dnk != 0.0);

    ip = dot(rmj, rnk);

    double acs = ip / (dmj * dnk);
    if (acs > 1.0)
    {
      if (acs < 1.0 + math::epsilon)
      {
        acs = 1.0;
      }
      else
      {
        io::messages.add("improper dihedral",
                         "acs > 1.0",
                         io::message::critical);
      }
    }

    phi = acos(acs);

    DEBUG(10, "raw phi=" << 180.0 * phi / math::Pi);

    ip = dot(rij, rnk);

    if (ip < 0)
      phi *= -1.0;

    DEBUG(9, "uncorrected phi=" << 180.0 * phi / math::Pi);

    double phi0_A = it->A_phi;
    double phi0_B = it->B_phi;
    double phi0 = (1 - l) * phi0_A + l * phi0_B;

    upper_bound = phi0 + it->delta;
    lower_bound = upper_bound - 2 * math::Pi;

    // bring the calculated value of phi to the interval between upper_bound and lower_bound
    while (phi < lower_bound)
      phi += 2 * math::Pi;
    while (phi > upper_bound)
      phi -= 2 * math::Pi;

    (*d_it) = phi;

    // in case delta is larger than 2*Pi, we need to do the same for phi_0
    // to be sure, we also do this for phi0_A and phi0_B
    // this may go wrong if you put phi0_A and phi0_B more than 2Pi apart
    while (phi0_A < lower_bound)
      phi0_A += 2 * math::Pi;
    while (phi0_A > upper_bound)
      phi0_A -= 2 * math::Pi;

    while (phi0_B < lower_bound)
      phi0_B += 2 * math::Pi;
    while (phi0_B > upper_bound)
      phi0_B -= 2 * math::Pi;

    while (phi0 < lower_bound)
      phi0 += 2 * math::Pi;
    while (phi0 > upper_bound)
      phi0 -= 2 * math::Pi;

    DEBUG(9, "phi=" << 180 * phi / math::Pi << " phi0=" << 180 * phi0 / math::Pi
                    << " delta=" << 180 * it->delta / math::Pi);
    DEBUG(9, "lower_bound =" << 180 * lower_bound / math::Pi
                             << " upper_bound =" << 180 * upper_bound / math::Pi);

    double delta_phi = phi - phi0;
    double phi_lin = sim.param().dihrest.phi_lin;
    double K = sim.param().dihrest.K;
    double A_K = sim.param().dihrest.K;
    double B_K = sim.param().dihrest.K;

    if (sim.param().dihrest.dihrest == simulation::dihedral_restr_inst_weighted)
    {
      K *= (1.0 - l) * it->A_w0 + l * it->B_w0;
      A_K *= it->A_w0;
      B_K *= it->B_w0;
    }

    double prefactor = pow(2.0, it->m + it->n) * pow(l, it->n) * pow(1.0 - l, it->m);

    if (phi_lin >= 0.0 && fabs(delta_phi) > phi_lin)
    {
      // LINEAR
      double zeta = 1;
      if (delta_phi < 0) zeta = -1;
      en_term = K * (zeta * delta_phi - 0.5 * phi_lin) * phi_lin;
      dlam_term = phi_lin * ((B_K - A_K) * (zeta * delta_phi - 0.5 * phi_lin) + K * zeta * (phi0_A - phi0_B));

      f = -prefactor * K * zeta * phi_lin;

    }
    else
    {
      // HARMONIC
      en_term = 0.5 * K * delta_phi * delta_phi;
      dlam_term = 0.5 * (B_K - A_K) * delta_phi * delta_phi + K * delta_phi * (phi0_A - phi0_B);
      f = -prefactor * K * delta_phi;

    }

    energy = prefactor * en_term;
    (*ene_it) = energy;

    conf.current().energies.dihrest_energy[topo.atom_energy_group()
                                               [it->i]] += energy;

    const double ki = f * dkj / dmj2;
    const double kl = -f * dkj / dnk2;
    const double kj1 = dot(rij, rkj) / dkj2 - 1.0;
    const double kj2 = dot(rkl, rkj) / dkj2;

    fi = ki * rmj;
    fl = kl * rnk;
    fj = kj1 * fi - kj2 * fl;
    fk = -1.0 * (fi + fj + fl);

    force(it->i) += fi;
    force(it->j) += fj;
    force(it->k) += fk;
    force(it->l) += fl;

    // lambda derivative

    // divide by zero measure
    double dprefndl = 0.0, dprefmdl = 0.0; 
    if (it->n == 0)
      dprefndl = 0;
    else
      dprefndl = it->n * pow(l, it->n - 1) * pow(1.0 - l, it->m);

    if (it->m == 0)
      dprefmdl = 0;
    else
      dprefmdl = it->m * pow(l, it->n) * pow(1.0 - l, it->m - 1);

    double dprefdl = pow(2.0, it->m + it->n) *
                     (dprefndl - dprefmdl) * en_term;

    double dpotdl = prefactor * dlam_term;

    energy_derivative = l_deriv * (dprefdl + dpotdl);

    conf.current().perturbed_energy_derivatives.dihrest_energy
        [topo.atom_energy_group()[it->i]] += energy_derivative;

    /**
     * ext_TI code - Betty
     */

    if (sim.param().precalclam.nr_lambdas &&
        ((sim.steps() % sim.param().write.free_energy) == 0))
    {

      double lambda_step = (sim.param().precalclam.max_lam -
                            sim.param().precalclam.min_lam) /
                           (sim.param().precalclam.nr_lambdas - 1);

      //loop over nr_lambdas
      for (unsigned int lam_index = 0; lam_index < sim.param().precalclam.nr_lambdas; ++lam_index)
      {

        double lam = (lam_index * lambda_step) + sim.param().precalclam.min_lam;
        double prefactorlam = pow(2.0, it->n + it->m) * pow(lam, it->n) * pow(1.0 - lam, it->m);
        double phi0lam = (1 - lam) * phi0_A + lam * phi0_B;
        double delta_philam = phi - phi0lam;
        double delta_philam2 = delta_philam * delta_philam;
        double K_diff = B_K - A_K;
        double phi0_diff = phi0_A - phi0_B;

        double en_termlam = 0.0, dlam_termlam = 0.0;
        if (phi_lin >= 0.0 && fabs(delta_philam) > phi_lin)
        {
          double zeta = 1;
          if (delta_philam < 0)
            zeta = -1;
          en_termlam = K * (zeta * delta_philam - 0.5 * phi_lin) * phi_lin;
          dlam_termlam = phi_lin * (K_diff * (zeta * delta_philam - 0.5 * phi_lin) + K * zeta * phi0_diff);
        }
        else
        {
          en_termlam = 0.5 * K * delta_philam2;
          dlam_termlam = 0.5 * K_diff * delta_philam2 + K * delta_philam * phi0_diff;
        }
        double energylam = prefactorlam * en_termlam;

        double dprefndlam = 0.0, dprefmdlam = 0.0;
        if (it->n == 0)
          dprefndlam = 0;
        else
          dprefndlam = it->n * pow(lam, it->n - 1) * pow(1.0 - lam, it->m);
        if (it->m == 0)
          dprefmdl = 0;
        else
          dprefmdlam = it->m * pow(lam, it->n) * pow(1.0 - lam, it->m - 1);

        double dprefdlam = pow(2.0, it->m + it->n) * (dprefndlam - dprefmdlam) * en_termlam;
        double dpotdlam = prefactorlam * dlam_termlam;
        double energy_derivativlam = dprefdlam + dpotdlam;

        conf.current().energies.AB_dihres[lam_index] += energylam;
        conf.current().perturbed_energy_derivatives.AB_dihres[lam_index] += energy_derivativlam;
      }
    }
    // Betty
  }
  return 0;
}

int interaction::Perturbed_Dihedral_Restraint_Interaction ::calculate_interactions(topology::Topology &topo,
                                                                                   configuration::Configuration &conf,
                                                                                   simulation::Simulation &sim)
{

  SPLIT_VIRIAL_BOUNDARY(_calculate_perturbed_dihedral_restraint_interactions,
                        topo, conf, sim);

  return 0;
}

/**
 * initiate dihedral restraint interactions
 */
template <math::boundary_enum B>
static void _init_dihres_data(topology::Topology &topo,
                              configuration::Configuration &conf)
{
  math::Periodicity<B> periodicity(conf.current().box);
  math::VArray &pos = conf.current().pos;

  math::Vec rij, rkj, rkl, rmj, rnk;
  double dkj2 = 0.0, dkj = 0.0, dmj2 = 0.0, dmj = 0.0, dnk2 = 0.0, dnk = 0.0, ip = 0.0, phi = 0.0;

  for (std::vector<topology::perturbed_dihedral_restraint_struct>::const_iterator
           it = topo.perturbed_dihedral_restraints().begin(),
           to = topo.perturbed_dihedral_restraints().end();
       it != to; ++it)
  {

    DEBUG(9, "init: perturbed dihedral angle " << it->i << "-" << it->j << "-" << it->k << "-" << it->l);

    periodicity.nearest_image(pos(it->i), pos(it->j), rij);
    periodicity.nearest_image(pos(it->k), pos(it->j), rkj);
    periodicity.nearest_image(pos(it->k), pos(it->l), rkl);
      
    bool warn=false;
    for (int i=0; i<3;  i++) {
        if ((fabs(rij[i]) > conf.current().box(i)[i]*0.45 && fabs(rij[i]) < conf.current().box(i)[i]*0.55)
         || (fabs(rkj[i]) > conf.current().box(i)[i]*0.45 && fabs(rkj[i]) < conf.current().box(i)[i]*0.55) 
         || (fabs(rkl[i]) > conf.current().box(i)[i]*0.45 && fabs(rkl[i]) < conf.current().box(i)[i]*0.55)) {
          warn=true;
        }
    }
    if (warn)
    {
      std::ostringstream oss;
      oss << "one or more vectors of your dihedral angle restraint are\n"
          << "close to half a box length long in your initial structure, \n"
          << "the dihedral might be calculated from other periodic copies of the atoms than you intended!\n";
      io::messages.add(oss.str(), "dihedral_restraint", io::message::warning);
    }

    rmj = cross(rij, rkj);
    rnk = cross(rkj, rkl);

    dkj2 = abs2(rkj);
    dmj2 = abs2(rmj);
    dnk2 = abs2(rnk);
    dkj = sqrt(dkj2);
    dmj = sqrt(dmj2);
    dnk = sqrt(dnk2);

    DEBUG(15, "dkj=" << dkj << " dmj=" << dmj << " dnk=" << dnk);

    assert(dmj != 0.0);
    assert(dnk != 0.0);

    ip = dot(rmj, rnk);

    double acs = ip / (dmj * dnk);
    if (acs > 1.0)
    {
      if (acs < 1.0 + math::epsilon)
      {
        acs = 1.0;
      }
      else
      {
        io::messages.add("improper dihedral",
                         "acs > 1.0",
                         io::message::critical);
      }
    }

    phi = acos(acs);

    DEBUG(10, "raw phi=" << 180.0 * phi / math::Pi);

    conf.special().pertdihedralres.d.push_back(phi);
    conf.special().pertdihedralres.energy.push_back(0.0);
  }
}

int interaction::Perturbed_Dihedral_Restraint_Interaction::init(topology::Topology &topo,
                                                                configuration::Configuration &conf,
                                                                simulation::Simulation &sim,
                                                                std::ostream &os,
                                                                bool quiet)
{

  SPLIT_BOUNDARY(_init_dihres_data, topo, conf);

  if (!quiet)
  {
    os << "Perturbed dihedral restraint interaction";
    os << std::endl;
  }
  return 0;
}
