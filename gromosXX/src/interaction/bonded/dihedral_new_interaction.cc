/**
 * @file dihedral_new_interaction.cc
 * template methods of Dihedral_new_Interaction.
 * calculates the dihedral forces for any m and any shift angle
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
#include "dihedral_new_interaction.h"

#include "../../util/template_split.h"
#include "../../util/debug.h"
#include "../../util/error.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE bonded

static double _calculate_nearest_minimum(double phi, int m, double cospd);

/**
 * calculate dihedral forces and energies.
 */
template<math::boundary_enum B, math::virial_enum V>
static int _calculate_dihedral_new_interactions(topology::Topology & topo,
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
  double dkj2 = 0.0, dkj = 0.0, dmj2 = 0.0, dnk2 = 0.0, dim = 0.0, dln = 0.0, ip = 0.0;
  double energy = 0.0;

  math::Periodicity<B> periodicity(conf.current().box);

  for (int n = 0; d_it != d_to; ++d_it, ++n) {
    periodicity.nearest_image(pos(d_it->i), pos(d_it->j), rij);
    periodicity.nearest_image(pos(d_it->k), pos(d_it->j), rkj);
    periodicity.nearest_image(pos(d_it->k), pos(d_it->l), rkl);
    rmj = cross(rij, rkj);
    rnk = cross(rkj, rkl);
    dmj2 = abs2(rmj);
    dnk2 = abs2(rnk);
    dkj2 = abs2(rkj);
    dkj = abs(rkj);
    
	double frim = dot(rij, rkj) / dkj2;
    double frln = dot(rkl, rkj) / dkj2;

    rim = rij - frim * rkj;
    rln = frln * rkj - rkl;
    dim = sqrt(abs2(rim));
    dln = sqrt(abs2(rln));

    ip = dot(rim, rln);
    double cosphi = ip / (dim * dln);
    if (cosphi > 1) cosphi = 1;
    if (cosphi < -1) cosphi = -1;
    double phi = acos(cosphi);
    double sign = dot(rij, rnk);
    if (sign < 0) phi *= -1.0;

    //DEBUG(10, "dihedral angle phi= " << phi);

    assert(unsigned(d_it->type) < param.size());

    double K = param[d_it->type].K;
    double delta = param[d_it->type].pd;
    double m = param[d_it->type].m;
    double cosdelta = param[d_it->type].cospd;

    DEBUG(10, "dihedral K=" << K << " delta=" << delta);


    double kj1 = frim - 1.0;
    double kj2 = frln;
    double ki = K * m * sin(m * phi - delta);
    double kl = -ki;

    fi = ki * dkj * rmj;
	if (dmj2 < (1.0e-10 * dkj2)) {
		fi = 0;
	    io::messages.add("One bond angle is close to 180 degrees!","dihedral_new_interaction",io::message::warning);
	} else {
       fi = fi / dmj2;
	}
	fl = kl * dkj * rnk;
	if (dnk2 < (1.0e-10 * dkj2)) {
		fl = 0;
		io::messages.add("One bond angle is close to 180 degrees!","dihedral_new_interaction",io::message::warning);
	} else {
	   fl = fl / dnk2;
	}
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

    energy = K * (1 + cos(m * phi - delta));
    conf.current().energies.dihedral_energy[topo.atom_energy_group()[d_it->i]] += energy;
    
    // ORIOL_GAMD
    if(sim.param().gamd.gamd){
      unsigned int gamd_group = topo.gamd_accel_group(d_it->i);
      std::vector<unsigned int> key = {gamd_group, gamd_group};
      unsigned int igroup = topo.gamd_interaction_group(key);
      DEBUG(10, "\tGAMD interaction group is " << igroup);
      conf.special().gamd.dihe_force[igroup][d_it->i] += fi;
      conf.special().gamd.dihe_force[igroup][d_it->j] += fj;
      conf.special().gamd.dihe_force[igroup][d_it->k] += fk;
      conf.special().gamd.dihe_force[igroup][d_it->l] += fl;
      conf.current().energies.gamd_dihedral_total[igroup] += energy;
      //conf.current().energies.gamd_potential_total[igroup] += energy;
      // virial
      for(int a=0; a<3; ++a){
        for(int bb=0; bb < 3; ++bb){
          conf.special().gamd.virial_tensor_dihe[igroup](a, bb) += rij(a) * fi(bb) + rkj(a) * fk(bb) +  rlj(a) * fl(bb);
        }
      }

    } // end gamd

    // dihedral angle monitoring.
    if (sim.param().print.monitor_dihedrals) {
      DEBUG(8, "monitoring dihedrals");

      DEBUG(11, "dihedral angle: " << phi
              << " previous minimum: " << conf.special().dihangle_trans.dihedral_angle_minimum[n]);

      if (fabs(conf.special().dihangle_trans.dihedral_angle_minimum[n] - phi) >
              2 * math::Pi / param[d_it->type].m) {
        double old_min = conf.special().dihangle_trans.dihedral_angle_minimum[n];
        conf.special().dihangle_trans.dihedral_angle_minimum[n] =
                _calculate_nearest_minimum(phi, param[d_it->type].m, cosdelta);
        // ugly check to see that it is not the first...
        if (fabs(old_min - 4 * math::Pi) > math::epsilon) {
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

int interaction::Dihedral_new_Interaction
::calculate_interactions(topology::Topology &topo,
        configuration::Configuration &conf,
        simulation::Simulation & sim) {
  m_timer.start(sim);

  SPLIT_VIRIAL_BOUNDARY(_calculate_dihedral_new_interactions,
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

