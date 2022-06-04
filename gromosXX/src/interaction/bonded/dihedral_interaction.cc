/**
 * @file dihedral_interaction.cc
 * template methods of Dihedral_interaction.
 */

#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"
#include "../../interaction/interaction.h"

#include "../../math/periodicity.h"

// interactions
#include "../../interaction/interaction_types.h"
#include "dihedral_interaction.h"

#include "../../util/template_split.h"
#include "../../util/debug.h"
#include "../../util/error.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE bonded

#include "../../util/debug.h"

static double _calculate_nearest_minimum(double phi, int m, double cospd);

/**
 * calculate dihedral forces and energies.
 */
template<math::boundary_enum B, math::virial_enum V>
static int _calculate_dihedral_interactions(topology::Topology & topo,
        configuration::Configuration & conf,
        simulation::Simulation & sim) {
  std::vector<interaction::dihedral_type_struct> const & param = topo.dihedral_types();
  // loop over the dihedrals
  std::vector<topology::four_body_term_struct>::iterator d_it =
          topo.solute().dihedrals().begin(),
          d_to = topo.solute().dihedrals().end();

  math::VArray &pos = conf.current().pos;
  math::VArray &force = conf.current().force;
  math::Vec rij, rkj, rkl, rlj, rim, rln, rmj, rnk, fi, fj, fk, fl;
  double dkj2 = 0.0, dim = 0.0, dln = 0.0, ip = 0.0;
  double energy = 0.0;

  math::Periodicity<B> periodicity(conf.current().box);

  for (int n = 0; d_it != d_to; ++d_it, ++n) {
    periodicity.nearest_image(pos(d_it->i), pos(d_it->j), rij);
    periodicity.nearest_image(pos(d_it->k), pos(d_it->j), rkj);
    periodicity.nearest_image(pos(d_it->k), pos(d_it->l), rkl);

    double cosdelta = param[d_it->type].cospd;


    rmj = cross(rij, rkj);
    rnk = cross(rkj, rkl);

    dkj2 = abs2(rkj);

    double frim = dot(rij, rkj) / dkj2;
    double frln = dot(rkl, rkj) / dkj2;

    rim = rij - frim * rkj;
    rln = frln * rkj - rkl;
    dim = sqrt(abs2(rim));
    dln = sqrt(abs2(rln));

    ip = dot(rim, rln);
    double cosphi = ip / (dim * dln);

    double cosphi2 = cosphi * cosphi;
    double cosphi3 = cosphi2 * cosphi;
    double cosphi4 = cosphi3 * cosphi;

    assert(unsigned(d_it->type) < param.size());

    double dcosmphi = 0;
    double cosmphi = 0;

    switch (param[d_it->type].m) {
      case 0:
        cosmphi = 0.0;
        dcosmphi = 0.0;
        break;
      case 1:
        cosmphi = cosphi;
        dcosmphi = 1;
        break;
      case 2:
        cosmphi = 2 * cosphi2 - 1;
        dcosmphi = 4 * cosphi;
        break;
      case 3:
        cosmphi = 4 * cosphi3 - 3 * cosphi;
        dcosmphi = 12 * cosphi2 - 3;
        break;
      case 4:
        cosmphi = 8 * cosphi4 - 8 * cosphi2 + 1;
        dcosmphi = 32 * cosphi3 - 16 * cosphi;
        break;
      case 5:
        cosmphi = 16 * cosphi4 * cosphi - 20 * cosphi3 + 5 * cosphi;
        dcosmphi = 80 * cosphi4 - 60 * cosphi2 + 5;
        break;
      case 6:
        cosmphi = 32 * cosphi4 * cosphi2 - 48 * cosphi4 + 18 * cosphi2 - 1;
        dcosmphi = 192 * cosphi4 * cosphi - 192 * cosphi3 + 36 * cosphi;
        break;
    }

    double K = param[d_it->type].K;

    DEBUG(10, "dihedral K=" << K << " cos(delta)=" << cosdelta << " dcos=" << dcosmphi);

    double ki = -K * cosdelta * dcosmphi / dim;
    double kl = -K * cosdelta * dcosmphi / dln;
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

    // if (V == math::atomic_virial){
    periodicity.nearest_image(pos(d_it->l), pos(d_it->j), rlj);

    for (int a = 0; a < 3; ++a)
      for (int bb = 0; bb < 3; ++bb)
        conf.current().virial_tensor(a, bb) +=
              rij(a) * fi(bb) +
        rkj(a) * fk(bb) +
        rlj(a) * fl(bb);

    DEBUG(11, "\tatomic virial done");
    // }

    energy = K * (1 + cosdelta * cosmphi);
    conf.current().energies.dihedral_energy
            [topo.atom_energy_group()[d_it->i]] += energy;

    // dihedral angle monitoring.
    if (sim.param().print.monitor_dihedrals) {
      DEBUG(8, "monitoring dihedrals");

      // cos_phi can be >1 or <-1 because of precision limits
      if (cosphi > 1) cosphi=1.0;
      if (cosphi < -1) cosphi=-1.0;

      double phi = acos(cosphi);

      ip = dot(rij, rnk);
      if (ip < 0) phi *= -1.0;
      DEBUG(11, "dihedral angle: " << phi
              << " previous minimum: " << conf.special().dihangle_trans.dihedral_angle_minimum[n]);

      if (fabs(conf.special().dihangle_trans.dihedral_angle_minimum[n] - phi) >
              2 * math::Pi / param[d_it->type].m) {
        double old_min = conf.special().dihangle_trans.dihedral_angle_minimum[n];
        conf.special().dihangle_trans.dihedral_angle_minimum[n] =
                _calculate_nearest_minimum(phi, param[d_it->type].m, cosdelta);
        // ugly check to see that it is not the first...
        if (old_min != 4 * math::Pi) {
          // could be written to a separate file or by a separate function
          // should at least be more descriptive.
          /*std::cout << "D-A-T: "
                << std::setw(4) << topo.solute().atom(d_it->i).residue_nr + 1
                << std::setw(4) << std::left
                << topo.residue_names()[topo.solute().atom(d_it->i).residue_nr]
                << std::setw(4)  << std::right<< topo.solute().atom(d_it->i).name
                << " -"
                << std::setw(4)  << std::right<< topo.solute().atom(d_it->j).name
                << " -"
                << std::setw(4)  << std::right<< topo.solute().atom(d_it->k).name
                << " -"
                << std::setw(4)  << std::right<< topo.solute().atom(d_it->l).name
                << std::setw(6) << d_it->i + 1 << " -"
                << std::setw(4) << d_it->j + 1 << " -"
                << std::setw(4) << d_it->k + 1 << " -"
                << std::setw(4) << d_it->l + 1
                << std::setw(10) << std::setprecision(1)
                << std::fixed << 180.0*old_min/math::Pi << " -> "
                << std::setw(8) << std::setprecision(1)
                << std::fixed
                << 180.0 * conf.special().dihedral_angle_minimum[n]/math::Pi << "\n";
          }*/
          conf.special().dihangle_trans.old_minimum[n] = old_min;
          conf.special().dihangle_trans.i[n] = d_it->i;
          conf.special().dihangle_trans.j[n] = d_it->j;
          conf.special().dihangle_trans.k[n] = d_it->k;
          conf.special().dihangle_trans.l[n] = d_it->l;
          conf.special().dihangle_trans.resid[n] = topo.solute().atom(d_it->i).residue_nr;
        } else {
          conf.special().dihangle_trans.old_minimum[n] = 0.0;
        }
      }
    }

  }
  return 0;

}

int interaction::Dihedral_Interaction
::calculate_interactions(topology::Topology &topo,
        configuration::Configuration &conf,
        simulation::Simulation &sim) {
  m_timer.start();

  SPLIT_VIRIAL_BOUNDARY(_calculate_dihedral_interactions,
          topo, conf, sim);

  m_timer.stop();

  return 0;
}

/**
 * calculate nearest minimum
 */
static inline double _calculate_nearest_minimum(double phi, int m, double cospd) {
  // copy from gromos++ nearest_minimum function
  double a_minimum = 0.5 * math::Pi * (3.0 - cospd) / m;
  double delta_phi = 2 * math::Pi / m;
  double nearest_min = a_minimum - int(rint((a_minimum - phi) / delta_phi)) * delta_phi;
  if (nearest_min >= 2 * math::Pi - math::epsilon) nearest_min -= 2 * math::Pi;

  return nearest_min;
}

int interaction::Dihedral_Interaction::init(topology::Topology &topo,
        configuration::Configuration &conf,
        simulation::Simulation &sim,
        std::ostream &os,
        bool quiet) {

  std::vector<interaction::dihedral_type_struct> const & param = topo.dihedral_types();
  std::vector<topology::four_body_term_struct>::iterator d_it =
          topo.solute().dihedrals().begin(),
          d_to = topo.solute().dihedrals().end();

  for (int n = 0; d_it != d_to; ++d_it, ++n) {
    double cosdelta = param[d_it->type].cospd;
    int m = param[d_it->type].m;

    if (((cosdelta > -1 - math::epsilon) && (cosdelta < -1 + math::epsilon)) || ((cosdelta > 1 - math::epsilon) && (cosdelta < 1 + math::epsilon))) {
      //
    } else {
      io::messages.add("dihedral function (NTBDN=0) not implemented for phase shifts not equal to 0 or 180 degrees. Please use NTBDN=1 in COVALENTFORM block", "dihedral interaction", io::message::error);
      return 1;
    }
    if (m > 6 || m < 0) {
      io::messages.add("dihedral function not implemented for m>6 or m<0. Please use NTBDN=1 in COVALENTFORM block", "dihedral_interaction", io::message::error);
      return 1;
    }
  }
  return 0;

};
