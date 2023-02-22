/**
 * @file nonbonded_innerloop.cc
 * template methods of Nonbonded_Innerloop
 */

#include <vector>
#include <set>

#include "storage.h"
#include "nonbonded_innerloop.h"
#include "topology/topology.h"


#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE nonbonded

#include "solvent_innerloop.cc"

template<typename t_nonbonded_spec>
inline void
interaction::Nonbonded_Innerloop<t_nonbonded_spec>::lj_crf_innerloop_2
    (
     topology::Topology & topo, 
     unsigned int i,
     unsigned int j,
        const double dist2,
        double &f,
       double &e_lj, double &e_crf
            ) {

  switch (t_nonbonded_spec::interaction_func) {
    case simulation::lj_crf_func:
    {
      const lj_parameter_struct & lj =
              m_param->lj_parameter(topo.iac(i),
              topo.iac(j));

      DEBUG(11, "\tlj-parameter c6=" << lj.c6 << " c12=" << lj.c12);
      DEBUG(11, "\tcharge i=" << topo.charge()(i) << " j=" << topo.charge()(j));

      lj_crf_interaction_fast(dist2, lj.c6, lj.c12,
              charge_product(topo, i, j),
              f, e_lj, e_crf);
      DEBUG(12, "f: " << f);
      DEBUG(12, "e_lj: " << e_lj);
      DEBUG(12, "e_crf: " << e_crf);

      break;
    }
   
    default:
      io::messages.add("Nonbonded_Innerloop",
              "interaction function not implemented (new version)",
              io::message::critical);
  }
}

template<typename t_nonbonded_spec>
inline void
interaction::Nonbonded_Innerloop<t_nonbonded_spec>::lj_crf_innerloop
(
        topology::Topology & topo,
        configuration::Configuration & conf,
        unsigned int i,
        unsigned int j,
        Storage & storage,
        math::Periodicity<t_nonbonded_spec::boundary_type> const & periodicity,
        unsigned int gamd
        ) {
  DEBUG(8, "\tpair\t" << i << "\t" << j);

  math::Vec r, force;
  double f = 0.0;
  double e_lj = 0.0, e_crf = 0.0;
  //ORIOL_GAMD
  unsigned int igroup = 0;
  unsigned int gamd_i, gamd_j = 0;
  if (gamd){
    unsigned int gamdi = topo.gamd_accel_group(i);
    unsigned int gamdj = topo.gamd_accel_group(j);
    std::vector<unsigned int> key = {gamdi, gamdj};
    igroup = topo.gamd_interaction_group(key);
  }

  periodicity.nearest_image(conf.current().pos(i),
          conf.current().pos(j), r);
  DEBUG(10, "\tni i " << conf.current().pos(i)(0) << " / "
          << conf.current().pos(i)(1) << " / "
          << conf.current().pos(i)(2));
  DEBUG(10, "\tni j " << conf.current().pos(j)(0) << " / "
          << conf.current().pos(j)(1) << " / "
          << conf.current().pos(j)(2));
  DEBUG(10, "\tni r " << r(0) << " / " << r(1) << " / " << r(2));

  switch (t_nonbonded_spec::interaction_func) {
    case simulation::lj_crf_func:
    {
      const lj_parameter_struct & lj =
              m_param->lj_parameter(topo.iac(i),
              topo.iac(j));

      DEBUG(11, "\tlj-parameter c6=" << lj.c6 << " c12=" << lj.c12);
      DEBUG(11, "\tcharge i=" << topo.charge()(i) << " j=" << topo.charge()(j));

      lj_crf_interaction(r, lj.c6, lj.c12,
              charge_product(topo, i, j),
              f, e_lj, e_crf);
      DEBUG(12, "f: " << f);
      DEBUG(12, "e_lj: " << e_lj);
      DEBUG(12, "e_crf: " << e_crf);

      DEBUG(10, "\t\tatomic virial");
      //ORIOL_GAMD
      for (int a = 0; a < 3; ++a) {
        force(a) = f * r(a);
        storage.force(i)(a) += force(a);
        storage.force(j)(a) -= force(a);
        if (gamd){
          storage.force_gamd[igroup](i)(a) += force(a);
          storage.force_gamd[igroup](j)(a) -= force(a);
        }

        for (int b = 0; b < 3; ++b){
          storage.virial_tensor(b, a) += r(b) * force(a);
          if (gamd){
            storage.virial_tensor_gamd[igroup](b, a) += r(b) * force(a);
          }
        }
      }
      break;
    }

    case simulation::cgrain_func:
    {
      DEBUG(11, "\tiac(i) = " << topo.iac(i) << " iac(j) = " << topo.iac(j));
      DEBUG(11, "\tcg_parameter.size() = " << m_param->cg_parameter().size());
      DEBUG(11, "\tcg_parameter()[i].size() = " << m_param->cg_parameter()[i].size());

      const lj_parameter_struct & cg =
              m_param->cg_parameter(topo.iac(i),
              topo.iac(j));

      DEBUG(11, "\tcg-parameter c6=" << cg.c6 << " c12=" << cg.c12);
      DEBUG(11, "\tcharge i=" << topo.charge()(i) << " j=" << topo.charge()(j));

      cgrain_interaction(r, cg.c6, cg.c12,
              topo.charge(i) *
              topo.charge(j),
              f, e_lj, e_crf);

      DEBUG(10, "\t\tatomic virial");
      //ORIOL_GAMD
      for (int a = 0; a < 3; ++a) {
        force(a) = f * r(a);
        storage.force(i)(a) += force(a);
        storage.force(j)(a) -= force(a);
        if (gamd){
          storage.force_gamd[igroup](i)(a) += force(a);
          storage.force_gamd[igroup](j)(a) -= force(a);
        }
        for (int b = 0; b < 3; ++b){
          storage.virial_tensor(b, a) += r(b) * force(a);
          if (gamd){
            storage.virial_tensor_gamd[igroup](b, a) += r(b) * force(a);
          }
        }
      }
      break;
    }

    case simulation::cggromos_func:
    {
      DEBUG(11, "\tiac(i) = " << topo.iac(i) << " iac(j) = " << topo.iac(j));
      const lj_parameter_struct & lj = m_param->lj_parameter(topo.iac(i),
              topo.iac(j));

      DEBUG(11, "\tlj-parameter c6=" << lj.c6 << " c12=" << lj.c12);
      DEBUG(11, "\tcharge i=" << topo.charge()(i) << " j=" << topo.charge()(j));
      DEBUG(11, "\tCG-CG epsilon = " << cgrain_eps[0] << " , FG-CG epsilon = " << cgrain_eps[1]);
      DEBUG(11, "\ti = " << i << " , j = " << j);

      // check if...
      if (topo.is_coarse_grained(i) && topo.is_coarse_grained(j)) {  // CG-CG interaction
        DEBUG(11, "CG-CG pair");
        lj_crf_interaction(r, lj.c6, lj.c12,
              topo.charge(i) * topo.charge(j) / cgrain_eps[0],
              f, e_lj, e_crf, 0);
      } else if (topo.is_coarse_grained(i) || topo.is_coarse_grained(j)) {  // FG-CG interaction
        //DEBUG(1, "i = " << i << " , j = " << j);
        lj_crf_interaction(r, lj.c6, lj.c12,
              topo.charge(i) * topo.charge(j) / cgrain_eps[1],
              f, e_lj, e_crf, 1);
      } else {     // FG-FG interaction
        //DEBUG(1, "FG: i = " << i << " , j = " << j);
        lj_crf_interaction(r, lj.c6, lj.c12,
              topo.charge(i) * topo.charge(j),
              f, e_lj, e_crf, 2);
        DEBUG(10, "\t\tatomic virial");
      }
      //ORIOL_GAMD
      for (int a = 0; a < 3; ++a) {
        force(a) = f * r(a);
        storage.force(i)(a) += force(a);
        storage.force(j)(a) -= force(a);
        if (gamd){
          storage.force_gamd[igroup](i)(a) += force(a);
          storage.force_gamd[igroup](j)(a) -= force(a);
        }
        for (int b = 0; b < 3; ++b){
          storage.virial_tensor(b, a) += r(b) * force(a);
          if (gamd){
            storage.virial_tensor_gamd[igroup](b, a) += r(b) * force(a);
          }
        }
      }
      break;
    }

    case simulation::pol_lj_crf_func:
    {
      if (!topo.is_polarisable(i) && !topo.is_polarisable(j)) {
        const lj_parameter_struct & lj =
                m_param->lj_parameter(topo.iac(i),
                topo.iac(j));

        DEBUG(11, "\tlj-parameter c6=" << lj.c6 << " c12=" << lj.c12);
        DEBUG(11, "\tcharge i=" << topo.charge()(i) << " j=" << topo.charge()(j));

        lj_crf_interaction(r, lj.c6, lj.c12,
                topo.charge(i) *
                topo.charge(j),
                f, e_lj, e_crf);

        DEBUG(10, "\t\tatomic virial");
        //ORIOL_GAMD
        for (int a = 0; a < 3; ++a) {
          force(a) = f * r(a);
          storage.force(i)(a) += force(a);
          storage.force(j)(a) -= force(a);
          if (gamd){
            storage.force_gamd[igroup](i)(a) += force(a);
            storage.force_gamd[igroup](j)(a) -= force(a);
          }
          for (int b = 0; b < 3; ++b){
            storage.virial_tensor(b, a) += r(b) * force(a);
            if (gamd){
              storage.virial_tensor_gamd[igroup](b, a) += r(b) * force(a);
            }
          }
        }
      } else {
        double f_pol[4];

        const math::Vec rp1(r - conf.current().posV(j));
        const math::Vec rp2(r + conf.current().posV(i));
        const math::Vec rpp(rp2 - conf.current().posV(j));

        DEBUG(10, "\t rpp " << rpp(0) << " / " << rpp(1) << " / " << rpp(2));

        const lj_parameter_struct & lj = m_param->lj_parameter(topo.iac(i), topo.iac(j));

        DEBUG(11, "\tlj-parameter c6=" << lj.c6 << " c12=" << lj.c12);
        DEBUG(11, "\tcharge i=" << topo.charge()(i) << " j=" << topo.charge()(j));
        DEBUG(11, "\tcoscharge i=" << topo.coscharge()(i) << " j=" << topo.coscharge()(j));

        pol_lj_crf_interaction(r, rp1, rp2, rpp, lj.c6, lj.c12,
                topo.charge(i), topo.charge(j),
                topo.coscharge(i), topo.coscharge(j),
                f_pol, e_lj, e_crf);

        DEBUG(10, "\tatomic virial");
        //ORIOL_GAMD
        for (int a = 0; a < 3; ++a) {
          force(a) = f_pol[0] * r(a) + f_pol[1] * rp1(a) + f_pol[2] * rp2(a) + f_pol[3] * rpp(a);
          storage.force(i)(a) += force(a);
          storage.force(j)(a) -= force(a);
          if (gamd){
            storage.force_gamd[igroup](i)(a) += force(a);
            storage.force_gamd[igroup](j)(a) -= force(a);
          }
          for (int b = 0; b < 3; ++b){
            storage.virial_tensor(b, a) += r(b) * force(a);
            if (gamd){
            storage.virial_tensor_gamd[igroup](b, a) += r(b) * force(a);
            }
          }
        }
      }
      break;
    }
    case simulation::pol_off_lj_crf_func:
    {
      if (!topo.is_polarisable(i) && !topo.is_polarisable(j)) {
        const lj_parameter_struct & lj =
                m_param->lj_parameter(topo.iac(i),
                topo.iac(j));

        DEBUG(11, "\tlj-parameter c6=" << lj.c6 << " c12=" << lj.c12);
        DEBUG(11, "\tcharge i=" << topo.charge()(i) << " j=" << topo.charge()(j));

        lj_crf_interaction(r, lj.c6, lj.c12,
                topo.charge(i) *
                topo.charge(j),
                f, e_lj, e_crf);

        DEBUG(10, "\t\tatomic virial");
        //ORIOL_GAMD
        for (int a = 0; a < 3; ++a) {
          force(a) = f * r(a);
          storage.force(i)(a) += force(a);
          storage.force(j)(a) -= force(a);
          if (gamd){
            storage.force_gamd[igroup](i)(a) += force(a);
            storage.force_gamd[igroup](j)(a) -= force(a);
          }
          for (int b = 0; b < 3; ++b){
            storage.virial_tensor(b, a) += r(b) * force(a);
            if (gamd){
            storage.virial_tensor_gamd[igroup](b, a) += r(b) * force(a);
            }
          }
        }
      } else {
        math::Vec rm=r;
        DEBUG(10, "\t topo.gamma(i) " << topo.gamma(i));
        if(topo.gamma(i)!=0.0){
         math::Vec rij, rik, rim;
         periodicity.nearest_image(conf.current().pos(i),
                            conf.current().pos(topo.gamma_j(i)), rij);
         periodicity.nearest_image(conf.current().pos(i),
                            conf.current().pos(topo.gamma_k(i)), rik);
         rim=topo.gamma(i)*(rij+rik)/2;
         rm-=rim;
        }
        DEBUG(10, "\t topo.gamma(j) " << topo.gamma(j));
        if(topo.gamma(j)!=0.0){
          math::Vec rjj, rjk,rjm;
          periodicity.nearest_image(conf.current().pos(j),
                            conf.current().pos(topo.gamma_j(j)), rjj);
          periodicity.nearest_image(conf.current().pos(j),
                            conf.current().pos(topo.gamma_k(j)), rjk);
          rjm=topo.gamma(j)*(rjj+rjk)/2;
          rm+=rjm;
        }
        
        double f_pol[5];

        const math::Vec rp1(rm - conf.current().posV(j));
        const math::Vec rp2(rm + conf.current().posV(i));
        const math::Vec rpp(rp2 - conf.current().posV(j));

        DEBUG(10, "\t rpp " << rpp(0) << " / " << rpp(1) << " / " << rpp(2));

        const lj_parameter_struct & lj = m_param->lj_parameter(topo.iac(i), topo.iac(j));
        
        DEBUG(11, "\tlj-parameter c6=" << lj.c6 << " c12=" << lj.c12);
        DEBUG(11, "\tcharge i=" << topo.charge()(i) << " j=" << topo.charge()(j));
        DEBUG(11, "\tcoscharge i=" << topo.coscharge()(i) << " j=" << topo.coscharge()(j));

        pol_off_lj_crf_interaction(r, rm, rp1, rp2, rpp, lj.c6, lj.c12,
                topo.charge(i), topo.charge(j),
                topo.coscharge(i), topo.coscharge(j),
                f_pol, e_lj, e_crf);
        
        DEBUG(10, "\tatomic virial");
        for (int a = 0; a < 3; ++a) {

         
          const double term = f_pol[1] * rm(a) +
                  f_pol[2] * rp1(a) + f_pol[3] * rp2(a) + f_pol[4] * rpp(a);
          storage.force(i)(a) +=(1-topo.gamma(i))*term+f_pol[0]*r(a);
          storage.force(j)(a) -=(1-topo.gamma(j))*term+f_pol[0]*r(a);
          if (gamd){
            storage.force_gamd[igroup](i)(a) +=(1-topo.gamma(i))*term+f_pol[0]*r(a);
            storage.force_gamd[igroup](j)(a) -=(1-topo.gamma(j))*term+f_pol[0]*r(a);
          }


          if(topo.gamma(i)!=0.0){
           storage.force(topo.gamma_j(i))(a) +=term*topo.gamma(i)/2;
           storage.force(topo.gamma_k(i))(a) +=term*topo.gamma(i)/2;
           if (gamd){
            storage.force_gamd[igroup](topo.gamma_j(i))(a) +=term*topo.gamma(i)/2;
            storage.force_gamd[igroup](topo.gamma_k(i))(a) +=term*topo.gamma(i)/2;
           }
          }
          if(topo.gamma(j)!=0.0){
           storage.force(topo.gamma_j(j))(a) -=term*topo.gamma(j)/2;
           storage.force(topo.gamma_k(j))(a) -=term*topo.gamma(j)/2;
           if (gamd){
            storage.force_gamd[igroup](topo.gamma_j(j))(a) +=term*topo.gamma(j)/2;
            storage.force_gamd[igroup](topo.gamma_k(j))(a) +=term*topo.gamma(j)/2;
           }
          }
          //bugbugbug, debug...
          for (int b = 0; b < 3; ++b){
        //    storage.virial_tensor(b, a) += (term+f_pol[0]*r(a))*r(b);
            storage.virial_tensor(b, a) += (term*rm(b))+(f_pol[0]*r(a))*r(b);
            if (gamd){
              storage.virial_tensor_gamd[igroup](b, a) += (term*rm(b))+(f_pol[0]*r(a))*r(b);
            }
          }
         }
        }
      break;
    }
    case simulation::lj_func :
    case simulation::lj_ls_func : {
      const lj_parameter_struct & lj =
              m_param->lj_parameter(topo.iac(i),
              topo.iac(j));

      DEBUG(11, "\tlj-parameter c6=" << lj.c6 << " c12=" << lj.c12);

      lj_interaction(r, lj.c6, lj.c12, f, e_lj);
      e_crf = 0.0;

      DEBUG(10, "\t\tatomic virial");
      //ORIOL_GAMD
      for (int a = 0; a < 3; ++a) {
        force(a) = f * r(a);
        storage.force(i)(a) += force(a);
        storage.force(j)(a) -= force(a);
        if (gamd){
          storage.force_gamd[igroup](i)(a) += force(a);
          storage.force_gamd[igroup](j)(a) -= force(a);
        }
        for (int b = 0; b < 3; ++b){
          storage.virial_tensor(b, a) += r(b) * force(a);
          if (gamd){
            storage.virial_tensor_gamd[igroup](b, a) += r(b) * force(a);
          }
        }
      }
      break;
    }

    case simulation::lj_shifted_crf_corr_func:
    {
      const lj_parameter_struct & lj =
              m_param->lj_parameter(topo.iac(i),
              topo.iac(j));

      DEBUG(11, "\tlj-parameter c6=" << lj.c6 << " c12=" << lj.c12);
      DEBUG(11, "\tcharge i=" << topo.charge()(i) << " j=" << topo.charge()(j));

      double e_extra_orig;
      double e_extra_phys;
      lj_shifted_crf_corr_interaction(r, lj.c6, lj.c12,
              charge_product(topo, i, j),
              f, e_lj, e_crf, e_extra_orig, e_extra_phys);
      DEBUG(12, "f: " << f);
      DEBUG(12, "e_lj: " << e_lj);
      DEBUG(12, "e_crf: " << e_crf);

      DEBUG(10, "\t\tatomic virial");
      for (int a = 0; a < 3; ++a) {
        force(a) = f * r(a);
        storage.force(i)(a) += force(a);
        storage.force(j)(a) -= force(a);

        for (int b = 0; b < 3; ++b)
          storage.virial_tensor(b, a) += r(b) * force(a);
      }
      storage.energies.shift_extra_orig[topo.atom_energy_group(i)][topo.atom_energy_group(j)] += e_extra_orig;
      storage.energies.shift_extra_phys[topo.atom_energy_group(i)][topo.atom_energy_group(j)] += e_extra_phys;
      break;
    }

    default:
      io::messages.add("Nonbonded_Innerloop",
              "interaction function not implemented",
              io::message::critical);
  }

  //////////////////////////////////////////////////
  // molecular virial is calculated by a
  // correction to the atomic virial
  //////////////////////////////////////////////////

  /*
  if (fabs(f) > 10000){
    std::cerr << "pair " << i << " - " << j << " force = " << f << std::endl;
    std::cerr << "\tiac " << topo.iac(i) << " " << topo.iac(j)
          << "\tq " << topo.charge(i) << " " << topo.charge(j)
          << "\n\tr " << math::v2s(r)
          << "\n\tr2 " << abs2(r) << std::endl;
  }
   */

  // energy
  assert(storage.energies.lj_energy.size() >
          topo.atom_energy_group(i));
  assert(storage.energies.lj_energy.size() >
          topo.atom_energy_group(j));

  DEBUG(11, "\tenergy group i " << topo.atom_energy_group(i)
          << " j " << topo.atom_energy_group(j));

  const unsigned int eg_i = topo.atom_energy_group(i);
  const unsigned int eg_j = topo.atom_energy_group(j);
  storage.energies.lj_energy[eg_i][eg_j] += e_lj;
  storage.energies.crf_energy[eg_i][eg_j] += e_crf;

  //ORIOL_GAMD
  if (gamd){
    storage.energies.gamd_potential_total[igroup] += e_lj;
    storage.energies.gamd_potential_total[igroup] += e_crf;
  }

#ifdef XXFORCEGROUPS
  if (storage.force_groups.size()) {
    storage.force_groups[eg_i][eg_j][i] += force;
    storage.force_groups[eg_i][eg_j][j] -= force;
  }
#endif
}

template<typename t_nonbonded_spec>
inline void
interaction::Nonbonded_Innerloop<t_nonbonded_spec>::lj_innerloop
(
        topology::Topology & topo,
        configuration::Configuration & conf,
        unsigned int i,
        unsigned int j,
        Storage & storage,
        math::Periodicity<t_nonbonded_spec::boundary_type> const & periodicity
        ) {
  DEBUG(8, "\tpair\t" << i << "\t" << j);

  math::Vec r, force;
  double f = 0.0;
  double e_lj = 0.0;

  periodicity.nearest_image(conf.current().pos(i),
          conf.current().pos(j), r);
  DEBUG(10, "\tni i " << conf.current().pos(i)(0) << " / "
          << conf.current().pos(i)(1) << " / "
          << conf.current().pos(i)(2));
  DEBUG(10, "\tni j " << conf.current().pos(j)(0) << " / "
          << conf.current().pos(j)(1) << " / "
          << conf.current().pos(j)(2));
  DEBUG(10, "\tni r " << r(0) << " / " << r(1) << " / " << r(2));

  const lj_parameter_struct & lj =
          m_param->lj_parameter(topo.iac(i),
          topo.iac(j));

  DEBUG(11, "\tlj-parameter c6=" << lj.c6 << " c12=" << lj.c12);

  lj_interaction(r, lj.c6, lj.c12, f, e_lj);

  DEBUG(10, "\t\tatomic virial");
  for (int a = 0; a < 3; ++a) {
  force(a) = f * r(a);
  storage.force(i)(a) += force(a);
  storage.force(j)(a) -= force(a);

  for (int b = 0; b < 3; ++b)
          storage.virial_tensor(b, a) += r(b) * force(a);
  }

// energy
  assert(storage.energies.lj_energy.size() >
          topo.atom_energy_group(i));
  assert(storage.energies.lj_energy.size() >
          topo.atom_energy_group(j));

  DEBUG(11, "\tenergy group i " << topo.atom_energy_group(i)
          << " j " << topo.atom_energy_group(j));

  const unsigned int eg_i = topo.atom_energy_group(i);
  const unsigned int eg_j = topo.atom_energy_group(j);
  storage.energies.lj_energy[eg_i][eg_j] += e_lj;

#ifdef XXFORCEGROUPS
  if (storage.force_groups.size()) {
    storage.force_groups[eg_i][eg_j][i] += force;
    storage.force_groups[eg_i][eg_j][j] -= force;
  }
#endif
}


template<typename t_nonbonded_spec>
inline void
interaction::Nonbonded_Innerloop<t_nonbonded_spec>::sasa_calc_innerloop
(
        topology::Topology & topo,
        configuration::Configuration & conf,
        unsigned int i,
        simulation::Simulation & sim,
        math::Periodicity<t_nonbonded_spec::boundary_type> const & periodicity) {

  bool higher = false;
  const double p12 = sim.param().sasa.p_12;
  const double p13 = sim.param().sasa.p_13;
  const double p1x = sim.param().sasa.p_1x;
  std::set< unsigned int >::const_iterator it, to;

  // get sasa parameters for atom i
  const topology::sasa_parameter_struct & sasa_param_i = topo.sasa_parameter(i);
  DEBUG(15, "\tSASA parameters for atom(i) " << sasa_param_i.atom << ": p = " <<
          sasa_param_i.p << " r = " << sasa_param_i.r <<
          " sigma = " << sasa_param_i.sigma);

  // reduce total sasa by overlap from each type of neighbour in turn

  // first (1-2) neighbours
  it = topo.sasa_first_neighbour(i).begin();
  to = topo.sasa_first_neighbour(i).end();
  for (; it != to; ++it) {
    DEBUG(12, "\tFirst neighbours: sasa atoms " << i << " and " << *it);
    calculate_sasa(topo, conf, higher, p12, i, *it, sasa_param_i, periodicity);
  } // end first neighbours

  // second (1-3) neighbours
  it = topo.sasa_second_neighbour(i).begin();
  to = topo.sasa_second_neighbour(i).end();
  for (; it != to; ++it) {
    DEBUG(12, "\tSecond neighbours: sasa atoms " << i << " and " << *it);
    calculate_sasa(topo, conf, higher, p13, i, *it, sasa_param_i, periodicity);
  } // end second neighbours

  // reset higher
  higher = true;

  // third (1-4) neighbours
  it = topo.sasa_third_neighbour(i).begin();
  to = topo.sasa_third_neighbour(i).end();
  for (; it != to; ++it) {
    DEBUG(12, "\tThird neighbours: sasa atoms " << i << " and " << *it);
    calculate_sasa(topo, conf, higher, p1x, i, *it, sasa_param_i, periodicity);
  } // end third neighbours

  // higher (> 1-4) neighbours
  it = topo.sasa_higher_neighbour(i).begin();
  to = topo.sasa_higher_neighbour(i).end();
  for (; it != to; ++it) {
    DEBUG(12, "\tHigher neighbours: sasa atoms " << i << " and " << *it);
    calculate_sasa(topo, conf, higher, p1x, i, *it, sasa_param_i, periodicity);
  } // end 1-4 and higher neighbours

} // end sasa_calc_innerloop

template<typename t_nonbonded_spec>
inline void
interaction::Nonbonded_Innerloop<t_nonbonded_spec>::calculate_sasa
(
        topology::Topology & topo,
        configuration::Configuration & conf,
        bool higher,
        const double & pij,
        unsigned int i, unsigned int j,
        const topology::sasa_parameter_struct & sasa_param_i,
        math::Periodicity<t_nonbonded_spec::boundary_type> const & periodicity) {

  // get sasa parameters for atom j
  const topology::sasa_parameter_struct & sasa_param_j = topo.sasa_parameter(j);
  DEBUG(15, "\tSASA parameters for atom(j) " << sasa_param_j.atom << ": p = "
          << sasa_param_j.p << " r = " << sasa_param_j.r <<
          " sigma = " << sasa_param_j.sigma);

  // assign/compute sums of radii
  const double ri_rh2o = sasa_param_i.r_rh2o;
  const double rj_rh2o = sasa_param_j.r_rh2o;
  const double sum_of_radii = ri_rh2o +  rj_rh2o;

  // compute distance between atoms i and j
  math::Vec rij;
  math::VArray &pos = conf.current().pos;
  periodicity.nearest_image(pos(sasa_param_i.atom), pos(sasa_param_j.atom), rij);
  const double rdist = math::abs(rij);

  DEBUG(15, "\tRi + Rh2o = " << ri_rh2o << "\tRj + Rh2o = " << rj_rh2o <<
          "\tSum of radii = " << sum_of_radii << "\trdist = " << rdist);

  // if 1st or 2nd neighbours we don't need to check the inter-atomic distance
  if (!higher) {

    // compute components of area function
    const double c1 = (sum_of_radii - rdist) * math::Pi;
    const double c2 = (rj_rh2o - ri_rh2o) / rdist;
    // note that "bij", "bji" are actually pi*pij*bij/Si
    const double bij = (c1 * ri_rh2o * (1.0 + c2) * pij * sasa_param_i.p) /
            sasa_param_i.surface;
    const double bji = (c1 * rj_rh2o * (1.0 - c2) * pij * sasa_param_j.p) /
            sasa_param_j.surface;

    DEBUG(15, "\tc1 = " << c1 << "\tc2 = " << c2 << "\tbij = " << bij << "\tbji = "
            << bji);

    // modify areas
    conf.current().sasa_area[i] *= (1.0 - bij);
    conf.current().sasa_area[j] *= (1.0 - bji);
    DEBUG(15, "\tCurrent actual SASA of atom " << sasa_param_i.atom << " is " <<
            conf.current().sasa_area[i] << " and of atom " << sasa_param_j.atom
            << " is " << conf.current().sasa_area[j]);

  } // end not higher

  // for third and higher neighbours, we only adjust the SASA if the two atoms
  // are closer than two waters
  if (higher) {
    if (rdist <= sum_of_radii) {

      // compute components of area function
      const double c1 = (sum_of_radii - rdist) * math::Pi;
      const double c2 = (rj_rh2o - ri_rh2o) / rdist;
      // note that "bij", "bji" are are actually pi*pij*bij/Si
      const double bij = (c1 * ri_rh2o * (1.0 + c2) * pij * sasa_param_i.p) /
              sasa_param_i.surface;
      const double bji = (c1 * rj_rh2o * (1.0 - c2) * pij * sasa_param_j.p) /
              sasa_param_j.surface;

      DEBUG(15, "\tc1 = " << c1 << "\tc2 = " << c2 << "\tbij = " << bij << "\tbji = "
              << bji);

      // modify areas
      conf.current().sasa_area[i] *= (1.0 - bij);
      conf.current().sasa_area[j] *= (1.0 - bji);
      DEBUG(15, "\tCurrent actual SASA of atom " << sasa_param_i.atom << " is " <<
              conf.current().sasa_area[i] << " and of atom " << sasa_param_j.atom
              << " is " << conf.current().sasa_area[j]);

    } // end cutoff check
  } // end higher

} // end sasa_calculation

template<typename t_nonbonded_spec>
inline void
interaction::Nonbonded_Innerloop<t_nonbonded_spec>::sasa_force_innerloop
(
        topology::Topology & topo,
        configuration::Configuration & conf,
        unsigned int i,
        simulation::Simulation & sim,
        math::Periodicity<t_nonbonded_spec::boundary_type> const & periodicity) {

  double dg_i = 0.0;
  bool higher = false;
  math::Vec force_i(0.0);
  math::VArray & force = conf.current().force;
  std::set< unsigned int >::const_iterator it, to;
  const double p12 = sim.param().sasa.p_12;
  const double p13 = sim.param().sasa.p_13;
  const double p1x = sim.param().sasa.p_1x;


  // get sasa parameters for atom i and sum of radii
  const topology::sasa_parameter_struct & sasa_param_i = topo.sasa_parameter(i);
  const double ri_rh2o = sasa_param_i.r_rh2o;
  DEBUG(15, "\n\tSASA parameters for atom " << sasa_param_i.atom << " : p = "
          << sasa_param_i.p << " r = " << sasa_param_i.r << " sigma = "
          << sasa_param_i.sigma);

  // and if using volume term, get derivative switching function for atom i
  if (sim.param().sasa.switch_volume) {
    dg_i = conf.current().dgvol[i];
    DEBUG(15, "\tVolume switching function derivative for atom(i) " << sasa_param_i.atom <<
            " with SASA = " << conf.current().sasa_area[i] << " is = " << dg_i);
  }

  // to ignore forces for bonded atoms (as in POPS), comment this loop out
  // first neighbours
  it = topo.sasa_first_neighbour(i).begin();
  to = topo.sasa_first_neighbour(i).end();
  for (; it != to; ++it) {
    DEBUG(12, "\tFirst neighbours: sasa atoms " << i << " and " << *it);
    calculate_sasa_forces(topo, conf, higher, p12, i, *it, force_i,
            ri_rh2o, dg_i, sasa_param_i, sim, periodicity);
  } // end first neighbours
  

  // second neighbours
  it = topo.sasa_second_neighbour(i).begin();
  to = topo.sasa_second_neighbour(i).end();
  for (; it != to; ++it) {
    DEBUG(12, "\tSecond neighbours: sasa atoms " << i << " and " << *it);
    calculate_sasa_forces(topo, conf, higher, p13, i, *it, force_i, ri_rh2o,
            dg_i, sasa_param_i, sim, periodicity);
  } // end second neighbours of i

  // third (1-4) neighbours
  higher = true;
  it = topo.sasa_third_neighbour(i).begin();
  to = topo.sasa_third_neighbour(i).end();
  for (; it != to; ++it) {
    DEBUG(12, "\tThird neighbours: sasa atoms " << i << " and " << *it);
    calculate_sasa_forces(topo, conf, higher, p1x, i, *it, force_i, ri_rh2o,
            dg_i, sasa_param_i, sim, periodicity);
  } // end second neighbours of i

  // higher (>1-4) neighbours
  it = topo.sasa_higher_neighbour(i).begin();
  to = topo.sasa_higher_neighbour(i).end();
  for (; it != to; ++it) {
    DEBUG(12, "\tHigher neighbours: sasa atoms " << i << " and " << *it);
    calculate_sasa_forces(topo, conf, higher, p1x, i, *it, force_i, ri_rh2o,
            dg_i, sasa_param_i, sim, periodicity);
  } // higher (>1-4) neighbours of i

  DEBUG(12, "\tSASA/VOL forces for atom(i) " << sasa_param_i.atom << " :" << math::v2s(force_i));

  // then adjust overall forces for atom i
  force(sasa_param_i.atom) += force_i;
  
} // end sasa force

// calculate the sasa contribution to the forces
template<typename t_nonbonded_spec>
inline void
interaction::Nonbonded_Innerloop<t_nonbonded_spec>::calculate_sasa_forces
(
        topology::Topology & topo,
        configuration::Configuration & conf,
        bool higher,
        const double pij,
        unsigned int i,
        unsigned int j,
        math::Vec & force_i,
        const double ri_rh2o,
        double dg_i,
        const topology::sasa_parameter_struct & sasa_param_i,
        simulation::Simulation & sim,
        math::Periodicity<t_nonbonded_spec::boundary_type> const & periodicity) {

  math::VArray & pos = conf.current().pos;
  math::VArray & force = conf.current().force;
  math::Vec rij(0.0);
  double ddf_vol = 0.0;

  // get sasa parameters for atom j
  const topology::sasa_parameter_struct & sasa_param_j = topo.sasa_parameter(j);
  DEBUG(15, "\tSASA parameters for atom(j) " << sasa_param_j.atom << " : p = "
          << sasa_param_j.p << " r = " << sasa_param_j.r << " sigma = "
          << sasa_param_j.sigma);

  // assign/compute sums and diffs of radii
  const double rj_rh2o = sasa_param_j.r_rh2o;
  const double rdiff = sasa_param_j.r - sasa_param_i.r;
  const double sum_of_radii = ri_rh2o + rj_rh2o;

  // get diff in coords and compute dist_rij2 and dist_rij
  periodicity.nearest_image(pos(sasa_param_i.atom), pos(sasa_param_j.atom), rij);
  const double dist_rij2 = math::abs2(rij);
  const double dist_rij = sqrt(dist_rij2);

  DEBUG(15, "\tRi + Rh2o = " << ri_rh2o << "\tRj + Rh2o = " << rj_rh2o <<
          "\n\tsum of radii = " << sum_of_radii << "\tdiff of radii = " << rdiff <<
          "\trdist2 = " << dist_rij2 << "\trdist = " << dist_rij <<
          "\n\tSASA(j) = " << conf.current().sasa_area[j] <<
          "\tinitial force(j): " << math::v2s(force(sasa_param_j.atom)));

  if (!higher) {
    // compute components of force contribution
    const double c1 = math::Pi * (sum_of_radii - dist_rij);
    const double c2 = rdiff / dist_rij;
    const double bij = c1 * ri_rh2o * (1.0 + c2);
    const double bji = c1 * rj_rh2o * (1.0 - c2);
    const double cc1 = conf.current().sasa_area[i] / ((sasa_param_i.surface /
            (pij * sasa_param_i.p)) - bij);
    const double cc2 = conf.current().sasa_area[j] / ((sasa_param_j.surface /
            (pij * sasa_param_j.p)) - bji);
    const double termij = sum_of_radii * rdiff / dist_rij2;
    const double cc3 = math::Pi * ri_rh2o * (1.0 + termij);
    const double cc4 = math::Pi * rj_rh2o * (1.0 - termij);

    const double ddf = (sasa_param_i.sigma * cc1 * cc3 + sasa_param_j.sigma * cc2 * cc4)
            / dist_rij;

    // compute the force
    math::Vec f = ddf * rij;

    DEBUG(15, "\tContributions to SASA forces:\n\tc1 = " << c1 << "\tc2 = " << c2
            << "\tbij = " << bij << "\tbji = " << bji
            << "\tcc1 = " << cc1 << "\tcc2 = " << cc2 << "\n\ttermij = " << termij <<
            "\tcc3 = " << cc3 << "\tcc4 = " << cc4 << "\tterm_sasa = " << ddf);

    // volume contributions
    if (sim.param().sasa.switch_volume) {
      // switching function for atom j
      const double dg_j = conf.current().dgvol[j];
      DEBUG(15, "\tVolume switching function derivative for atom(j) " << sasa_param_j.atom <<
              " with SASA of " << conf.current().sasa_area[j] << " is " << dg_j);

      // extracting the dA/dr part in each direction
      // sasa_vol is zero if this atom doesn't contribute to the volume term
      const double part_i = sim.param().sasa.sigma_v * cc1 * cc3 * conf.current().sasa_buriedvol[i] * dg_i;
      const double part_j = sim.param().sasa.sigma_v * cc2 * cc4 * conf.current().sasa_buriedvol[j] * dg_j;
      ddf_vol = (part_i + part_j) / dist_rij;
      DEBUG(15, "\tVolume on, ddf_vol = " << ddf_vol);

      // update the force with volume contribution
      f += ddf_vol * rij;

    } // end volume

    // update the forces on j and store temporarily the forces on i
    force(sasa_param_j.atom) += f;
    force_i -= f;
    DEBUG(12, "\tUpdated force for atom(j) " << sasa_param_j.atom << ": " <<
            math::v2s(force(sasa_param_j.atom)));

  }// end not higher
  
  else {

    // check we are within cutoff
    if (dist_rij <= sum_of_radii) {

      // compute components of force contribution
      const double c1 = math::Pi * (sum_of_radii - dist_rij);
      const double c2 = rdiff / dist_rij;
      const double bij = c1 * ri_rh2o * (1.0 + c2);
      const double bji = c1 * rj_rh2o * (1.0 - c2);
      const double cc1 = conf.current().sasa_area[i] / ((sasa_param_i.surface /
              (pij * sasa_param_i.p)) - bij);
      const double cc2 = conf.current().sasa_area[j] / ((sasa_param_j.surface /
              (pij * sasa_param_j.p)) - bji);
      const double termij = sum_of_radii * rdiff / dist_rij2;
      const double cc3 = math::Pi * ri_rh2o * (1.0 + termij);
      const double cc4 = math::Pi * rj_rh2o * (1.0 - termij);

      const double ddf = (sasa_param_i.sigma * cc1 * cc3 + sasa_param_j.sigma * cc2 * cc4)
              / dist_rij;

      // compute the force
      math::Vec f = ddf * rij;

      DEBUG(15, "\tContributions to SASA forces:\n\tc1 = " << c1 << "\tc2 = " << c2
              << "\tbij = " << bij << "\tbji = " << bji
              << "\tcc1 = " << cc1 << "\tcc2 = " << cc2 << "\n\ttermij = " << termij <<
              "\tcc3 = " << cc3 << "\tcc4 = " << cc4 << "\tterm_sasa = " << ddf);

      // volume contributions
      if (sim.param().sasa.switch_volume) {
        // switching function for atom j
        const double dg_j = conf.current().dgvol[j];
        DEBUG(15, "\tVolume switching function deriv. for atom(j) " << sasa_param_j.atom <<
                " with SASA of " << conf.current().sasa_area[j] << " is " << dg_j);

        // extracting the dA/dr part in each direction
        const double part_i = sim.param().sasa.sigma_v * cc1 * cc3 * conf.current().sasa_buriedvol[i] * dg_i;
        const double part_j = sim.param().sasa.sigma_v * cc2 * cc4 * conf.current().sasa_buriedvol[j] * dg_j;
        ddf_vol = (part_i + part_j) / dist_rij;
        DEBUG(15, "\tVolume on, ddf_vol = " << ddf_vol);

        // update the force with volume contribution
        f += ddf_vol * rij;

      } // end volume

      // update forces on atoms
      force(sasa_param_j.atom) += f;
      force_i -= f;
      DEBUG(12, "\tUpdated force for atom(j) " << sasa_param_j.atom << ": " <<
              math::v2s(force(sasa_param_j.atom))); 
    } // end cutoff
  } // end higher

} // end calculate_sasa_forces

template<typename t_nonbonded_spec>
inline void
interaction::Nonbonded_Innerloop<t_nonbonded_spec>::sasa_volume_innerloop
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 unsigned int i,
 simulation::Simulation & sim)
{
  // get sasa parameters for atom i
  const topology::sasa_parameter_struct & sasa_param_i = topo.sasa_parameter(i);
  // get switching function for atom i
  const double g = conf.current().gvol[i];

  // only store the volume of buried atoms
  // (those that will contribute to the VOL, not SASA, energy term)
  if (g != 0.0) {
    conf.current().sasa_buriedvol[i] = sasa_param_i.vol;
    conf.current().sasa_buriedvol_tot += sasa_param_i.vol;
    // compute volume energy term for atom i
    const double evolume = sim.param().sasa.sigma_v * g * conf.current().sasa_buriedvol[i];
    DEBUG(15, "\tvol. energy of atom " << sasa_param_i.atom << " is " << evolume);
    // add volume energy to energy group of atom i (this is zeroed at start)
    conf.current().energies.sasa_volume_energy[topo.atom_energy_group(sasa_param_i.atom)] += evolume;
  } else {
    // make sure vol is zero
    conf.current().sasa_buriedvol[i] = 0.0;
  }

  DEBUG(15, "\tVolume of atom " << sasa_param_i.atom << " is " << conf.current().sasa_buriedvol[i]
          << " and current sum of volume energy for this energy group is " <<
          conf.current().energies.sasa_volume_energy[topo.atom_energy_group(sasa_param_i.atom)]);
}

template<typename t_nonbonded_spec>
void interaction::Nonbonded_Innerloop<t_nonbonded_spec>::sasa_switching_fct
(
        configuration::Configuration & conf,
        unsigned int i,
        const double amin,
        const double amax,
        const double adiff) {

  // get sasa of atom i
  const double a = conf.current().sasa_area[i];

  // sasa is less than lowest cutoff: it contributes to the volume term
  if (a <= amin) { //return 1.0;
    conf.current().gvol[i] = 1.0;
    conf.current().dgvol[i] = 0.0; // not strictly necessary
  }
  // sasa is between the lower and upper cutoffs: use switching function
  else if (a <= amax) {
    const double da = a - amin;
    const double da2 = da * da;
    const double adiff2 = adiff * adiff;
    const double adiff3 = adiff2 * adiff;
    const double g = ((2.0 * da * da2) / (adiff3)) - ((3.0 * da2) / (adiff2)) + 1.0;
    DEBUG(15, "\tSwitching function is " << g);
    conf.current().gvol[i] = g;
    const double dg = ((6.0 * da2) / (adiff3)) - ((6.0 * da) / (adiff2));
    DEBUG(15, "\tDerivative of switching function is " << dg);
    conf.current().dgvol[i] = dg;
  }
  else {
  // sasa is greater than highest cutoff: it contributes to the sasa term instead
  conf.current().gvol[i] = 0.0; // already zero
  conf.current().dgvol[i] = 0.0; // already zero
  }
}

template<typename t_nonbonded_spec>
void interaction::Nonbonded_Innerloop<t_nonbonded_spec>::one_four_interaction_innerloop
(
        topology::Topology & topo,
        configuration::Configuration & conf,
        int i,
        int j,
        Storage & storage,
        math::Periodicity<t_nonbonded_spec::boundary_type> const & periodicity,
        unsigned int gamd) {
  DEBUG(8, "\t1,4-pair\t" << i << "\t" << j);

  math::Vec r;
  double f = 0.0, e_lj = 0.0, e_crf = 0.0, e_ls = 0.0;
  //ORIOL_GAMD
  unsigned int igroup = 0;
  unsigned int gamd_i, gamd_j = 0;
  if (gamd){
    unsigned int gamdi = topo.gamd_accel_group(i);
    unsigned int gamdj = topo.gamd_accel_group(j);
    std::vector<unsigned int> key = {gamdi, gamdj};
    igroup = topo.gamd_interaction_group(key);
  }

  periodicity.nearest_image(conf.current().pos(i),
          conf.current().pos(j), r);
  DEBUG(10, "\tni i " << conf.current().pos(i)(0) << " / "
          << conf.current().pos(i)(1) << " / "
          << conf.current().pos(i)(2));
  DEBUG(10, "\tni j " << conf.current().pos(j)(0) << " / "
          << conf.current().pos(j)(1) << " / "
          << conf.current().pos(j)(2));
  DEBUG(10, "\tni r " << r(0) << " / " << r(1) << " / " << r(2));
  switch (t_nonbonded_spec::interaction_func) {
    case simulation::lj_crf_func:
    {
      const lj_parameter_struct & lj =
              m_param->lj_parameter(topo.iac(i),
              topo.iac(j));

      DEBUG(11, "\tlj-parameter cs6=" << lj.cs6 << " cs12=" << lj.cs12);
      DEBUG(11, "\tcharge i=" << topo.charge()(i) << " j=" << topo.charge()(j));
      DEBUG(11, "\tcoulomb scaling = " << m_param->get_coulomb_scaling());

      lj_crf_interaction(r, lj.cs6, lj.cs12,
              charge_product(topo, i, j),
              f, e_lj, e_crf, 0,
              m_param->get_coulomb_scaling());

      DEBUG(10, "\t\tatomic virial");
      //ORIOL_GAMD
      for (int a = 0; a < 3; ++a) {
        const double term = f * r(a);
        storage.force(i)(a) += term;
        storage.force(j)(a) -= term;
        if (gamd){
          storage.force_gamd[igroup](i)(a) += term;
          storage.force_gamd[igroup](j)(a) -= term;
        }
        for (int b = 0; b < 3; ++b){
          storage.virial_tensor(b, a) += r(b) * term;
          if (gamd){
            storage.virial_tensor_gamd[igroup](b, a) += r(b) * term;
          }
        }
      }
      break;
    }
    case simulation::cgrain_func:
    {
      io::messages.add("Nonbonded_Innerloop",
              "no 1,4 interactions for Martini coarse-grained simulations!",
              io::message::critical);
      break;
    }
    case simulation::cggromos_func :
    {
      // check if...
      if (topo.is_coarse_grained(i)) { // CG-CG or CG-FG interaction
        io::messages.add("Nonbonded_Innerloop",
                "no 1,4 interactions for Gromos coarse-grained simulations!",
                io::message::critical);
      } else { // FG-FG interaction
        const lj_parameter_struct & lj =
                m_param->lj_parameter(topo.iac(i),
                topo.iac(j));

        DEBUG(11, "\tlj-parameter cs6=" << lj.cs6 << " cs12=" << lj.cs12);
        DEBUG(11, "\tcharge i=" << topo.charge()(i) << " j=" << topo.charge()(j));
        
        lj_crf_interaction(r, lj.cs6, lj.cs12,
                topo.charge()(i) *
                topo.charge()(j),
                f, e_lj, e_crf, 2);

        DEBUG(10, "\t\tatomic virial");
      //ORIOL_GAMD
      for (int a = 0; a < 3; ++a) {
        const double term = f * r(a);
        storage.force(i)(a) += term;
        storage.force(j)(a) -= term;
        if (gamd){
          storage.force_gamd[igroup](i)(a) += term;
          storage.force_gamd[igroup](j)(a) -= term;
        }

        for (int b = 0; b < 3; ++b){
          storage.virial_tensor(b, a) += r(b) * term;
          if (gamd){
            storage.virial_tensor_gamd[igroup](b, a) += r(b) * term;
          }
        }
       }
      }
      break;
    }
    case simulation::pol_lj_crf_func :
      {
        double f_pol[4];
        
        const math::Vec rp1(r - conf.current().posV(j));
        const math::Vec rp2(r + conf.current().posV(i));
        const math::Vec rpp(rp2 - conf.current().posV(j));
        
        DEBUG(10, "\t rpp " << rpp(0) << " / " << rpp(1) << " / " << rpp(2));
        
        const lj_parameter_struct & lj = m_param->lj_parameter(topo.iac(i), topo.iac(j));
	
      	DEBUG(11, "\tlj-parameter cs6=" << lj.cs6 << " cs12=" << lj.cs12);
	      DEBUG(11, "\tcharge i=" << topo.charge()(i) << " j=" << topo.charge()(j));
        DEBUG(11, "\tcoscharge i=" << topo.coscharge()(i) << " j=" << topo.coscharge()(j));
        
        pol_lj_crf_interaction(r, rp1, rp2, rpp, lj.cs6, lj.cs12,
                               topo.charge(i), topo.charge(j),
                               topo.coscharge(i), topo.coscharge(j), 
                               f_pol, e_lj, e_crf);
                
        DEBUG(10, "\tatomic virial");
        for (int a=0; a<3; ++a){
          const double term = f_pol[0]*r(a) + f_pol[1]*rp1(a) + f_pol[2]*rp2(a) + f_pol[3]*rpp(a);
          storage.force(i)(a) += term;
          storage.force(j)(a) -= term;
          if (gamd){
            storage.force_gamd[igroup](i)(a) += term;
            storage.force_gamd[igroup](j)(a) -= term;
          }
          for(int b=0; b<3; ++b){
            storage.virial_tensor(b, a) += r(b)*term;
            if (gamd){
              storage.virial_tensor_gamd[igroup](b, a) += r(b)*term;
            }
          }
        }
      break;
    }
    case simulation::pol_off_lj_crf_func :
      {
         math::Vec rm=r;
        if(topo.gamma(i)!=0.0){
        math::Vec rij, rik, rim;
        periodicity.nearest_image(conf.current().pos(i),
                            conf.current().pos(topo.gamma_j(i)), rij);
        periodicity.nearest_image(conf.current().pos(i),
                            conf.current().pos(topo.gamma_k(i)), rik);
        rim=topo.gamma(i)*(rij+rik)/2;
        rm-=rim;
        }
        if(topo.gamma(j)!=0.0){
          math::Vec rjj, rjk,rjm;
          periodicity.nearest_image(conf.current().pos(j),
                            conf.current().pos(topo.gamma_j(j)), rjj);
          periodicity.nearest_image(conf.current().pos(j),
                            conf.current().pos(topo.gamma_k(j)), rjk);
          rjm=topo.gamma(j)*(rjj+rjk)/2;
          rm+=rjm;
        }
        
        double f_pol[5];

        const math::Vec rp1(rm - conf.current().posV(j));
        const math::Vec rp2(rm + conf.current().posV(i));
        const math::Vec rpp(rp2 - conf.current().posV(j));

        DEBUG(10, "\t rpp " << rpp(0) << " / " << rpp(1) << " / " << rpp(2));

        const lj_parameter_struct & lj = m_param->lj_parameter(topo.iac(i), topo.iac(j));

        DEBUG(11, "\tlj-parameter cs6=" << lj.cs6 << " cs12=" << lj.cs12);
        DEBUG(11, "\tcharge i=" << topo.charge()(i) << " j=" << topo.charge()(j));
        DEBUG(11, "\tcoscharge i=" << topo.coscharge()(i) << " j=" << topo.coscharge()(j));

        pol_off_lj_crf_interaction(r,rm, rp1, rp2, rpp, lj.cs6, lj.cs12,
                               topo.charge(i), topo.charge(j),
                               topo.coscharge(i), topo.coscharge(j),
                               f_pol, e_lj, e_crf);

        DEBUG(10, "\tatomic virial");
        for (int a=0; a<3; ++a){
          const double term = f_pol[1]*r(a) + f_pol[2]*rp1(a) + f_pol[3]*rp2(a) + f_pol[4]*rpp(a);
          storage.force(i)(a) +=(1-topo.gamma(i))*term+f_pol[0]*r(a);
          storage.force(j)(a) -=(1-topo.gamma(j))*term+f_pol[0]*r(a);
          if (gamd){
            storage.force_gamd[igroup](i)(a) +=(1-topo.gamma(i))*term+f_pol[0]*r(a);
            storage.force_gamd[igroup](j)(a) -=(1-topo.gamma(j))*term+f_pol[0]*r(a);
          }
          if(topo.gamma(i)!=0.0){
           storage.force(topo.gamma_j(i))(a) +=topo.gamma(i)/2*term;
           storage.force(topo.gamma_k(i))(a) -=topo.gamma(i)/2*term;
           if (gamd){
            storage.force_gamd[igroup](topo.gamma_j(i))(a) +=topo.gamma(i)/2*term;
            storage.force_gamd[igroup](topo.gamma_k(i))(a) -=topo.gamma(i)/2*term;
           }
          }
          if(topo.gamma(j)!=0.0){
           storage.force(topo.gamma_j(j))(a) +=topo.gamma(j)/2*term;
           storage.force(topo.gamma_k(j))(a) -=topo.gamma(j)/2*term;
           if (gamd){
            storage.force_gamd[igroup](topo.gamma_j(j))(a) +=topo.gamma(j)/2*term;
            storage.force_gamd[igroup](topo.gamma_k(j))(a) -=topo.gamma(j)/2*term;
           }
          }
          for(int b=0; b<3; ++b){
           // storage.virial_tensor(b,a ) += (term+f_pol[0]*r(a))*r(b);
            storage.virial_tensor(b, a) += (term*rm(b))+(f_pol[0]*r(a))*r(b);
            if (gamd){
              storage.virial_tensor_gamd[igroup](b, a) += (term*rm(b))+(f_pol[0]*r(a))*r(b);
            }
          }
        }
        break;
    }
    case simulation::lj_func : {
      const lj_parameter_struct & lj =
              m_param->lj_parameter(topo.iac(i),
              topo.iac(j));

      DEBUG(11, "\tlj-parameter cs6=" << lj.cs6 << " cs12=" << lj.cs12);
      lj_interaction(r, lj.cs6, lj.cs12, f, e_lj);

      e_crf = 0.0;

      DEBUG(10, "\t\tatomic virial");
      for (int a = 0; a < 3; ++a) {
        const double term = f * r(a);
        storage.force(i)(a) += term;
        storage.force(j)(a) -= term;

        for (int b = 0; b < 3; ++b)
          storage.virial_tensor(b, a) += r(b) * term;
      }
      break;
    }
    case simulation::lj_ls_func : {
      const lj_parameter_struct & lj =
              m_param->lj_parameter(topo.iac(i),
              topo.iac(j));

      DEBUG(11, "\tlj-parameter cs6=" << lj.cs6 << " cs12=" << lj.cs12);
 //     lj_interaction(r, lj.cs6, lj.cs12, f, e_lj);

//  call lj_ls!!!
       lj_ls_interaction(r, lj.cs6, lj.cs12,
              topo.charge(i) * topo.charge(j),
              f, e_lj, e_ls);



      e_crf = 0.0;

      DEBUG(10, "\t\tatomic virial");
      //ORIOL_GAMD
      for (int a = 0; a < 3; ++a) {
        const double term = f * r(a);
        storage.force(i)(a) += term;
        storage.force(j)(a) -= term;
        if (gamd){
          storage.force_gamd[igroup](i)(a) += term;
          storage.force_gamd[igroup](j)(a) -= term;
        }
        for (int b = 0; b < 3; ++b){
          storage.virial_tensor(b, a) += r(b) * term;
          if (gamd){
            storage.virial_tensor_gamd[igroup](b, a) += r(b) * term;
          }
        }
      }
      break;
    }
    case simulation::lj_shifted_crf_corr_func:
    {
      const lj_parameter_struct & lj =
              m_param->lj_parameter(topo.iac(i),
              topo.iac(j));

      DEBUG(11, "\tlj-parameter cs6=" << lj.cs6 << " cs12=" << lj.cs12);
      DEBUG(11, "\tcharge i=" << topo.charge()(i) << " j=" << topo.charge()(j));
      DEBUG(11, "\tcoulomb scaling = " << m_param->get_coulomb_scaling());

      double e_extra_orig;
      double e_extra_phys;

      lj_shifted_crf_corr_interaction(r, lj.cs6, lj.cs12,
              charge_product(topo, i, j),
              f, e_lj, e_crf, e_extra_orig, e_extra_phys,
              0, m_param->get_coulomb_scaling());

      DEBUG(10, "\t\tatomic virial");
      for (int a = 0; a < 3; ++a) {
        const double term = f * r(a);
        storage.force(i)(a) += term;
        storage.force(j)(a) -= term;

        for (int b = 0; b < 3; ++b)
          storage.virial_tensor(b, a) += r(b) * term;
      }
      storage.energies.shift_extra_orig[topo.atom_energy_group(i)][topo.atom_energy_group(j)] += e_extra_orig;
      storage.energies.shift_extra_phys[topo.atom_energy_group(i)][topo.atom_energy_group(j)] += e_extra_phys;
      break;
    }
    default:
      io::messages.add("Nonbonded_Innerloop",
              "interaction function not implemented",
              io::message::critical);
  }

  // energy
  storage.energies.lj_energy[topo.atom_energy_group(i)]
          [topo.atom_energy_group(j)] += e_lj;

  storage.energies.crf_energy[topo.atom_energy_group(i)]
          [topo.atom_energy_group(j)] += e_crf;

  storage.energies.ls_real_energy[topo.atom_energy_group(i)]
          [topo.atom_energy_group(j)] += e_ls;

  // ORIOL_GAMD
  if (gamd){
    storage.energies.gamd_potential_total[igroup] += e_lj;
    storage.energies.gamd_potential_total[igroup] += e_crf;
  }
  DEBUG(11, "\tenergy group i " << topo.atom_energy_group(i)
          << " j " << topo.atom_energy_group(j));

}

template<typename t_nonbonded_spec>
void interaction::Nonbonded_Innerloop<t_nonbonded_spec>::lj_exception_innerloop
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 topology::lj_exception_struct const & ljex,
 Storage & storage,
 math::Periodicity<t_nonbonded_spec::boundary_type> const & periodicity,
 unsigned int gamd)
{
  DEBUG(8, "\tLJ exception\t" << ljex.i << "\t" << ljex.j);

  math::Vec r;
  double f = 0.0, e_lj = 0.0, e_crf = 0.0, e_ls = 0.0;

  unsigned int i = ljex.i, j = ljex.j;
  //ORIOL_GAMD
  unsigned int igroup = 0;
  if (gamd){
    unsigned int gamdi = topo.gamd_accel_group(i);
    unsigned int gamdj = topo.gamd_accel_group(j);
    std::vector<unsigned int> key = {gamdi, gamdj};
    igroup = topo.gamd_interaction_group(key);
  }

  periodicity.nearest_image(conf.current().pos(i),
			    conf.current().pos(j), r);
  DEBUG(10, "\tni i " << conf.current().pos(i)(0) << " / "
          << conf.current().pos(i)(1) << " / "
          << conf.current().pos(i)(2));
  DEBUG(10, "\tni j " << conf.current().pos(j)(0) << " / "
          << conf.current().pos(j)(1) << " / "
          << conf.current().pos(j)(2));
  DEBUG(10, "\tni r " << r(0) << " / " << r(1) << " / " << r(2));

  // do cutoff check
  if (math::abs2(r) > m_cut2)
    return;

  switch(t_nonbonded_spec::interaction_func){
    case simulation::lj_crf_func :
      {
	DEBUG(11, "\tlj-parameter cs6=" << ljex.c6 << " cs12=" << ljex.c12);
	DEBUG(11, "\tcharge i=" << topo.charge()(i) << " j=" << topo.charge()(j));

        lj_crf_interaction(r, ljex.c6, ljex.c12,
			 topo.charge()(i) *
			 topo.charge()(j),
			 f, e_lj, e_crf);

        DEBUG(10, "\t\tatomic virial");
      //ORIOL_GAMD
      for (int a = 0; a < 3; ++a) {
        const double term = f * r(a);
        storage.force(i)(a) += term;
        storage.force(j)(a) -= term;
        if (gamd){
          storage.force_gamd[igroup](i)(a) += term;
          storage.force_gamd[igroup](j)(a) -= term;
        }
        for (int b = 0; b < 3; ++b){
          storage.virial_tensor(b, a) += r(b) * term;
          if (gamd){
            storage.virial_tensor_gamd[igroup](b, a) += r(b) * term;
          }
        }
      }
      break;
    }

    case simulation::cggromos_func :
    {
      if (topo.is_coarse_grained(i) || topo.is_coarse_grained(j)) { // CG-CG or CG-FG
        io::messages.add("Nonbonded_Innerloop",
                         "No LJ exceptions for coarse-grained simulations!",
		         io::message::critical);
      } else {
        lj_crf_interaction(r, ljex.c6, ljex.c12,
			 topo.charge()(i) *
			 topo.charge()(j),
			 f, e_lj, e_crf, 2);

        DEBUG(10, "\t\tatomic virial");
      //ORIOL_GAMD
      for (int a = 0; a < 3; ++a) {
        const double term = f * r(a);
        storage.force(i)(a) += term;
        storage.force(j)(a) -= term;
        if (gamd){
          storage.force_gamd[igroup](i)(a) += term;
          storage.force_gamd[igroup](j)(a) -= term;
        }
        for (int b = 0; b < 3; ++b){
          storage.virial_tensor(b, a) += r(b) * term;
          if (gamd){
            storage.virial_tensor_gamd[igroup](b, a) += r(b) * term;
          }
        }
       }
      }
      break;
    }
    case simulation::cgrain_func :
      {
        io::messages.add("Nonbonded_Innerloop",
                         "No LJ exceptions for coarse-grained simulations!",
		         io::message::critical);
        break;
      }
    case simulation::pol_lj_crf_func :
      {

        double f_pol[4];
        const math::Vec rp1(r - conf.current().posV(j));
        const math::Vec rp2(r + conf.current().posV(i));
        const math::Vec rpp(rp2 - conf.current().posV(j));

        DEBUG(10, "\t rpp " << rpp(0) << " / " << rpp(1) << " / " << rpp(2));
	      DEBUG(11, "\tlj-parameter cs6=" << ljex.c6 << " cs12=" << ljex.c12);
      	DEBUG(11, "\tcharge i=" << topo.charge()(i) << " j=" << topo.charge()(j));
        DEBUG(11, "\tcoscharge i=" << topo.coscharge()(i) << " j=" << topo.coscharge()(j));

        pol_lj_crf_interaction(r, rp1, rp2, rpp, ljex.c6, ljex.c12,
                               topo.charge(i), topo.charge(j),
                               topo.coscharge(i), topo.coscharge(j),
                               f_pol, e_lj, e_crf);

        DEBUG(10, "\tatomic virial");
        //ORIOL_GAMD
        for (int a=0; a<3; ++a){
          const double term = f_pol[0]*r(a) + f_pol[1]*rp1(a) + f_pol[2]*rp2(a) + f_pol[3]*rpp(a);
          storage.force(i)(a) += term;
          storage.force(j)(a) -= term;
          if (gamd){
            storage.force_gamd[igroup](i)(a) += term;
            storage.force_gamd[igroup](j)(a) -= term;
          }
          for(int b=0; b<3; ++b){
            storage.virial_tensor(b, a) += r(b)*term;
            if (gamd){
            storage.virial_tensor_gamd[igroup](b, a) += r(b)*term;
            }
          }
        }

	break;
    }
   case simulation::pol_off_lj_crf_func :
      {
        math::Vec rm=r;
        if(topo.gamma(i)!=0.0){
        math::Vec rij, rik, rim;
        periodicity.nearest_image(conf.current().pos(i),
                            conf.current().pos(topo.gamma_j(i)), rij);
        periodicity.nearest_image(conf.current().pos(i),
                            conf.current().pos(topo.gamma_k(i)), rik);
        rim=topo.gamma(i)*(rij+rik)/2;
        rm-=rim;
        }
        if(topo.gamma(j)!=0.0){
          math::Vec rjj, rjk,rjm;
          periodicity.nearest_image(conf.current().pos(j),
                            conf.current().pos(topo.gamma_j(j)), rjj);
          periodicity.nearest_image(conf.current().pos(j),
                            conf.current().pos(topo.gamma_k(j)), rjk);
          rjm=topo.gamma(j)*(rjj+rjk)/2;
          rm+=rjm;
        }

        double f_pol[5];
        const math::Vec rp1(rm - conf.current().posV(j));
        const math::Vec rp2(rm + conf.current().posV(i));
        const math::Vec rpp(rp2 - conf.current().posV(j));

        DEBUG(10, "\t rpp " << rpp(0) << " / " << rpp(1) << " / " << rpp(2));
        DEBUG(11, "\tlj-parameter cs6=" << ljex.c6 << " cs12=" << ljex.c12);
        DEBUG(11, "\tcharge i=" << topo.charge()(i) << " j=" << topo.charge()(j));
        DEBUG(11, "\tcoscharge i=" << topo.coscharge()(i) << " j=" << topo.coscharge()(j));

        pol_off_lj_crf_interaction(r, rm, rp1, rp2, rpp, ljex.c6, ljex.c12,
                               topo.charge(i), topo.charge(j),
                               topo.coscharge(i), topo.coscharge(j),
                               f_pol, e_lj, e_crf);

        DEBUG(10, "\tatomic virial");
        for (int a=0; a<3; ++a){
          const double term = f_pol[1]*r(a) + f_pol[2]*rp1(a) + f_pol[3]*rp2(a) + f_pol[4]*rpp(a);
          storage.force(i)(a) +=(1-topo.gamma(i))*term+f_pol[0]*r(a);
          storage.force(j)(a) -=(1-topo.gamma(j))*term+f_pol[0]*r(a);
          if (gamd){
            storage.force_gamd[igroup](i)(a) +=(1-topo.gamma(i))*term+f_pol[0]*r(a);
            storage.force_gamd[igroup](j)(a) -=(1-topo.gamma(j))*term+f_pol[0]*r(a);
          }
          if(topo.gamma(i)!=0.0){
           storage.force(topo.gamma_j(i))(a) +=topo.gamma(i)/2*term;
           storage.force(topo.gamma_k(i))(a) -=topo.gamma(i)/2*term;
           if (gamd){
            storage.force_gamd[igroup](topo.gamma_j(i))(a) +=topo.gamma(i)/2*term;
            storage.force_gamd[igroup](topo.gamma_k(i))(a) -=topo.gamma(i)/2*term;
           }
          }
          if(topo.gamma(j)!=0.0){
           storage.force(topo.gamma_j(j))(a) +=topo.gamma(j)/2*term;
           storage.force(topo.gamma_k(j))(a) -=topo.gamma(j)/2*term;
           if (gamd){
            storage.force_gamd[igroup](topo.gamma_j(j))(a) +=topo.gamma(j)/2*term;
            storage.force_gamd[igroup](topo.gamma_k(j))(a) -=topo.gamma(j)/2*term;
           }
          }
          for(int b=0; b<3; ++b){
          // storage.virial_tensor(b,a ) += (term+f_pol[0]*r(a))*r(b);
           storage.virial_tensor(b, a) += (term*rm(b))+(f_pol[0]*r(a))*r(b);
           if (gamd){
            storage.virial_tensor_gamd[igroup](b, a) += (term*rm(b))+(f_pol[0]*r(a))*r(b);
           }
          }
        }
        break;
    }
    case simulation::lj_func : 
    case simulation::lj_ls_func : {
      DEBUG(11, "\tlj-parameter c6=" << ljex.c6 << " c12=" << ljex.c12);
      lj_interaction(r, ljex.c6, ljex.c12, f, e_lj);
      e_crf = 0.0;

      DEBUG(10, "\t\tatomic virial");
      //ORIOL_GAMD
      for (int a = 0; a < 3; ++a) {
        const double term = f * r(a);
        storage.force(i)(a) += term;
        storage.force(j)(a) -= term;
        if (gamd){
          storage.force_gamd[igroup](i)(a) += term;
          storage.force_gamd[igroup](j)(a) -= term;
        }
        for (int b = 0; b < 3; ++b){
          storage.virial_tensor(b, a) += r(b) * term;
          if (gamd){
            storage.virial_tensor_gamd[igroup](b, a) += r(b) * term;
          }
        }
      }
      break;
    }
    default:
      io::messages.add("Nonbonded_Innerloop",
		       "interaction function not implemented",
		       io::message::critical);
  }

  // energy
  storage.energies.lj_energy[topo.atom_energy_group(i)]
    [topo.atom_energy_group(j)] += e_lj;

  storage.energies.crf_energy[topo.atom_energy_group(i)]
    [topo.atom_energy_group(j)] += e_crf;

  storage.energies.ls_real_energy[topo.atom_energy_group(i)]
    [topo.atom_energy_group(j)] += e_ls;

  //ORIOL_GAMD
  if (gamd){
    storage.energies.gamd_potential_total[igroup] +=  e_lj;
    storage.energies.gamd_potential_total[igroup] +=  e_crf;
  }
  DEBUG(11, "\tenergy group i " << topo.atom_energy_group(i)
	<< " j " << topo.atom_energy_group(j));

}

template<typename t_nonbonded_spec>
inline void 
interaction::Nonbonded_Innerloop<t_nonbonded_spec>::RF_excluded_interaction_innerloop
(
        topology::Topology & topo,
        configuration::Configuration & conf,
        int i,
        Storage & storage,
        math::Periodicity<t_nonbonded_spec::boundary_type> const & periodicity,
        unsigned int gamd) {
  math::Vec r, f;
  double e_crf = 0.0;
  unsigned int gamd_i, gamd_j = 0;

  math::VArray &pos = conf.current().pos;
  math::VArray &force = storage.force;

  topology::excl_cont_t::value_type::const_iterator it, to;
  it = topo.exclusion(i).begin();
  to = topo.exclusion(i).end();

  DEBUG(8, "\tself-term " << i);
  r = 0;

  switch (t_nonbonded_spec::interaction_func) {

    case simulation::lj_crf_func:
    {
      // this will only contribute in the energy, the force should be zero.
      double q = topo.charge()(i);
      if (t_nonbonded_spec::charge_type == simulation::qm_buffer_charge) {
        q += topo.qm_delta_charge(i);
      }
      rf_interaction(r, q * q, f, e_crf);
      storage.energies.crf_energy[topo.atom_energy_group(i)]
              [topo.atom_energy_group(i)] += 0.5 * e_crf;
      // ORIOL_GAMD
      unsigned int igroup = 0;
      if (gamd){
        unsigned int gamd_i = topo.gamd_accel_group(i);
        std::vector<unsigned int> key = {gamd_i, gamd_i};
        igroup = topo.gamd_interaction_group(key);
        storage.energies.gamd_potential_total[igroup] +=  0.5 * e_crf;
      }
      DEBUG(11, "\tcontribution " << 0.5 * e_crf);

      for (; it != to; ++it) {

        //ORIOL_GAMD
        if (gamd){
          unsigned int gamd_j = topo.gamd_accel_group(*it);
          std::vector<unsigned int> key = {gamd_i, gamd_j};
          igroup = topo.gamd_interaction_group(key);
        }

        DEBUG(11, "\texcluded pair " << i << " - " << *it);

        periodicity.nearest_image(pos(i), pos(*it), r);
        DEBUG(10, "\tni i " << pos(i)(0) << " / "
                << pos(i)(1) << " / "
                << pos(i)(2));
        DEBUG(10, "\tni j " << pos(*it)(0) << " / "
                << pos(*it)(1) << " / "
                << pos(*it)(2));
        DEBUG(10, "\tni r " << r(0) << " / " << r(1) << " / " << r(2));
        rf_interaction(r, charge_product(topo, i, *it), f, e_crf);

        force(i) += f;
        force(*it) -= f;

        //ORIOL_GAMD
        if (gamd){
          storage.force_gamd[igroup](i) += f;
          storage.force_gamd[igroup](*it) -= f;
        }
        for (int a = 0; a < 3; ++a){
          for (int b = 0; b < 3; ++b){
            storage.virial_tensor(a, b) += r(a) * f(b);
            if (gamd){
              storage.virial_tensor_gamd[igroup](a,b) += r(a) * f(b);
            }
          }
        }
        DEBUG(11, "\tatomic virial done");

        // energy
        storage.energies.crf_energy[topo.atom_energy_group(i)]
                [topo.atom_energy_group(*it)] += e_crf;

        // ORIOL_GAMD
        if (gamd){
          storage.energies.gamd_potential_total[igroup] +=  e_crf;
        }

        DEBUG(11, "\tcontribution " << e_crf);

      } // loop over excluded pairs
      break;
    }
    case simulation::cgrain_func:
    {
      io::messages.add("Nonbonded_Innerloop",
              "no RF excluded interactions for Martini coarse-grained simulations!",
              io::message::critical);
      break;
    }
    case simulation::cggromos_func:
    {
      // this will only contribute in the energy, the force should be zero.
      // check if...
      if (topo.is_coarse_grained(i)) {  // CG particle
        rf_interaction(r, topo.charge()(i) * topo.charge()(i) / cgrain_eps[0],
                f, e_crf, 0);
      } else {     // FG particle
        rf_interaction(r, topo.charge()(i) * topo.charge()(i), f, e_crf, 2);
      }
      storage.energies.crf_energy[topo.atom_energy_group(i)]
              [topo.atom_energy_group(i)] += 0.5 * e_crf;

      // ORIOL_GAMD
      unsigned int igroup = 0;
      if (gamd){
        unsigned int gamd_i = topo.gamd_accel_group(i);
        std::vector<unsigned int> key = {gamd_i, gamd_i};
        igroup = topo.gamd_interaction_group(key);
        storage.energies.gamd_potential_total[igroup] +=  0.5 * e_crf;
      }
      DEBUG(11, "\tcontribution " << 0.5 * e_crf);

      for (; it != to; ++it) {

        DEBUG(11, "\texcluded pair " << i << " - " << *it);

        periodicity.nearest_image(pos(i), pos(*it), r);
        DEBUG(10, "\tni i " << pos(i)(0) << " / "
                << pos(i)(1) << " / "
                << pos(i)(2));
        DEBUG(10, "\tni j " << pos(*it)(0) << " / "
                << pos(*it)(1) << " / "
                << pos(*it)(2));
        DEBUG(10, "\tni r " << r(0) << " / " << r(1) << " / " << r(2));

        // check if...
        if (topo.is_coarse_grained(i) && topo.is_coarse_grained(*it)) { // CG-CG interaction
          rf_interaction(r, topo.charge()(i) * topo.charge()(*it) / cgrain_eps[0],
                  f, e_crf, 0);
        } else if (topo.is_coarse_grained(i) || topo.is_coarse_grained(*it)) { // FG-CG interaction
          rf_interaction(r, topo.charge()(i) * topo.charge()(*it) / cgrain_eps[1],
                  f, e_crf, 1);
        } else { // FG-FG interaction
          rf_interaction(r, topo.charge()(i) * topo.charge()(*it), f, e_crf, 2);
        }
        //ORIOL_GAMD
        if (gamd){
          unsigned int gamd_j = topo.gamd_accel_group(*it);
          std::vector<unsigned int> key = {gamd_i, gamd_j};
          igroup = topo.gamd_interaction_group(key);
          storage.force_gamd[igroup](i) += f;
          storage.force_gamd[igroup](*it) -= f;
        }
        force(i) += f;
        force(*it) -= f;

        for (int a = 0; a < 3; ++a){
          for (int b = 0; b < 3; ++b){
            storage.virial_tensor(a, b) += r(a) * f(b);
            if (gamd){
              storage.virial_tensor_gamd[igroup](a, b) += r(a) * f(b);
            }
          }
        }
        DEBUG(11, "\tatomic virial done");

        // energy
        storage.energies.crf_energy[topo.atom_energy_group(i)]
                [topo.atom_energy_group(*it)] += e_crf;

        // ORIOL_GAMD
        if (gamd){
          storage.energies.gamd_potential_total[igroup] +=  e_crf;
        }
        DEBUG(11, "\tcontribution " << e_crf);

      } // loop over excluded pairs
      break;
    }
    case simulation::pol_lj_crf_func:
    {
      math::Vec rp1, rp2, rpp;
      math::VArray f_pol(4);
      f_pol = 0.0;
      //ORIOL_GAMD
      unsigned int igroup = 0;
      if (gamd){
        unsigned int gamd_i = topo.gamd_accel_group(i);
        std::vector<unsigned int> key = {gamd_i, gamd_i};
        igroup = topo.gamd_interaction_group(key);
      }
      DEBUG(8, "\tself-term " << i);
      rp1 = -conf.current().posV(i);
      rp2 = conf.current().posV(i);
      rpp = 0.0;

      // this will only contribute in the energy, the force should be zero.
      pol_rf_interaction(r, rp1, rp2, rpp, topo.charge(i), topo.charge(i),
              topo.coscharge(i), topo.coscharge(i), f_pol, e_crf);
      storage.energies.crf_energy[topo.atom_energy_group(i)]
              [topo.atom_energy_group(i)] += 0.5 * e_crf;

      // ORIOL_GAMD
      if (gamd){
        storage.energies.gamd_potential_total[igroup] +=  0.5 * e_crf; 
      storage.energies.gamd_potential_total[igroup] +=  0.5 * e_crf; 
        storage.energies.gamd_potential_total[igroup] +=  0.5 * e_crf; 
      }
      DEBUG(11, "\tcontribution " << 0.5 * e_crf);

      for (; it != to; ++it) {

        DEBUG(11, "\texcluded pair " << i << " - " << *it);

        periodicity.nearest_image(pos(i), pos(*it), r);
        rp1 = r - conf.current().posV(*it);
        rp2 = r + conf.current().posV(i);
        rpp = rp2 - conf.current().posV(*it);

        pol_rf_interaction(r, rp1, rp2, rpp, topo.charge()(i),
                topo.charge()(*it), topo.coscharge(i),
                topo.coscharge(*it), f_pol, e_crf);
        
        //ORIOL_GAMD
        if (gamd){
          unsigned int gamd_j = topo.gamd_accel_group(*it);
          std::vector<unsigned int> key = {gamd_i, gamd_j};
          igroup = topo.gamd_interaction_group(key);
        }
        for (int a = 0; a < 3; ++a) {
          const double term = f_pol(0)(a) + f_pol(1)(a) + f_pol(2)(a) + f_pol(3)(a);
          force(i)(a) += term;
          force(*it)(a) -= term;
          storage.force_gamd[igroup](i)(a) += term;
          storage.force_gamd[igroup](*it)(a) -= term;
          for (int b = 0; b < 3; ++b){
            storage.virial_tensor(b, a) += r(b) * term;
            if (gamd){
              storage.virial_tensor_gamd[igroup](b, a) += r(b) * term;
            }
          }
        }
        DEBUG(11, "\tatomic virial done");

        // energy
        storage.energies.crf_energy[topo.atom_energy_group(i)]
                [topo.atom_energy_group(*it)] += e_crf;
        
        //ORIOL_GAMD
        if (gamd){
          storage.energies.gamd_potential_total[igroup] += e_crf;
        }
        DEBUG(11, "\tcontribution " << e_crf);
      } // loop over excluded pairs
      break;
    }
    case simulation::pol_off_lj_crf_func : {
      math::Vec rp1, rp2, rpp;
      math::VArray f_pol(4);
      f_pol = 0.0;
      //ORIOL_GAMD
      unsigned int igroup = 0;
      if (gamd){
        unsigned int gamd_i = topo.gamd_accel_group(i);
        std::vector<unsigned int> key = {gamd_i, gamd_i};
        igroup = topo.gamd_interaction_group(key);
      }
      DEBUG(8, "\tself-term " << i );
      rp1 = -conf.current().posV(i);
      //rp1 = 0.0;
      rp2 = conf.current().posV(i);
      rpp = 0.0;
      math::Vec rm=r;
      math::Vec rim= math::Vec(0.0);
      if(topo.gamma(i)!=0.0){
        math::Vec rij, rik;
        periodicity.nearest_image(conf.current().pos(i),
                            conf.current().pos(topo.gamma_j(i)), rij);
        periodicity.nearest_image(conf.current().pos(i),
                            conf.current().pos(topo.gamma_k(i)), rik);
        rim=topo.gamma(i)*(rij+rik)/2;
       }
      
      // this will only contribute in the energy, the force should be zero.
      pol_rf_interaction(r, rp1, rp2, rpp, topo.charge(i), topo.charge(i),
              topo.coscharge(i), topo.coscharge(i), f_pol, e_crf);
      storage.energies.crf_energy[topo.atom_energy_group(i)]
              [topo.atom_energy_group(i)] += 0.5 * e_crf;

      //ORIOL_GAMD
      if (gamd){
        storage.energies.gamd_potential_total[igroup] += 0.5 * e_crf;
      }
      DEBUG(11, "\tcontribution " << 0.5*e_crf);

      for( ; it != to; ++it){
        //ORIOL_GAMD
        if (gamd){
          unsigned int gamd_j = topo.gamd_accel_group(*it);
          std::vector<unsigned int> key = {gamd_i, gamd_j};
          igroup = topo.gamd_interaction_group(key);
        }
        periodicity.nearest_image(pos(i), pos(*it), r);
	math::Vec rjm= math::Vec(0.0);
        if(topo.gamma(*it)!=0.0){
          math::Vec rjj, rjk;
          periodicity.nearest_image(conf.current().pos(*it),
                            conf.current().pos(topo.gamma_j(*it)), rjj);
          periodicity.nearest_image(conf.current().pos(*it),
                            conf.current().pos(topo.gamma_k(*it)), rjk);
          rjm+=topo.gamma(*it)*(rjj+rjk)/2;
        }
        rm=r-rim+rjm; 
        DEBUG(11, "\texcluded pair " << i << " - " << *it);

        rp1 = rm - conf.current().posV(*it);
        rp2 = rm + conf.current().posV(i);
        rpp = rp2 - conf.current().posV(*it);

        pol_rf_interaction(rm, rp1, rp2, rpp, topo.charge()(i),
                topo.charge()(*it), topo.coscharge(i),
                topo.coscharge(*it), f_pol, e_crf);

        for (int a = 0; a < 3; ++a) {
          const double term = f_pol(0)(a) + f_pol(1)(a) + f_pol(2)(a) + f_pol(3)(a);
          force(i)(a)   +=(1-topo.gamma(i))*term;
          force(*it)(a) -=(1-topo.gamma(*it))*term;
          if (gamd){
            storage.force_gamd[igroup](i)(a)   +=(1-topo.gamma(i))*term;
            storage.force_gamd[igroup](*it)(a) -=(1-topo.gamma(*it))*term;
          }
          if(topo.gamma(i)!=0.0){
           force(topo.gamma_j(i))(a) +=topo.gamma(i)/2*term;
           force(topo.gamma_k(i))(a) +=topo.gamma(i)/2*term;
           if (gamd){
            storage.force_gamd[igroup](topo.gamma_j(i))(a) +=topo.gamma(i)/2*term;
            storage.force_gamd[igroup](topo.gamma_k(i))(a) +=topo.gamma(i)/2*term;
           }
          }
          if(topo.gamma(*it)!=0.0){
           force(topo.gamma_j(*it))(a) -=topo.gamma(*it)/2*term;
           force(topo.gamma_k(*it))(a) -=topo.gamma(*it)/2*term;
           if (gamd){
            storage.force_gamd[igroup](topo.gamma_j(*it))(a) -=topo.gamma(*it)/2*term;
            storage.force_gamd[igroup](topo.gamma_k(*it))(a) -=topo.gamma(*it)/2*term;
           }
          }
          for(int b=0; b<3; ++b){
            storage.virial_tensor(b, a) += rm(b)*term;
            if (gamd){ 
              storage.virial_tensor_gamd[igroup](b, a) += rm(b)*term;
            }
          }
        }
        DEBUG(11, "\tatomic virial done");

        // energy
        storage.energies.crf_energy[topo.atom_energy_group(i)]
                [topo.atom_energy_group(*it)] += e_crf;

        //ORIOL_GAMD
        if (gamd){
          storage.energies.gamd_potential_total[igroup] += 0.5 * e_crf;
        }
        DEBUG(11, "\tcontribution " << e_crf);
      } // loop over excluded pairs
      break;
    }
    case simulation::lj_shifted_crf_corr_func:
    {
      // this will only contribute in the energy, the force should be zero.
      double q = topo.charge()(i);
      if (t_nonbonded_spec::charge_type == simulation::qm_buffer_charge) {
        q += topo.qm_delta_charge(i);
      }
      double e_extra_orig;
      double e_extra_phys;

      shifted_rf_corr_interaction(r, q * q, f, e_crf, e_extra_orig, e_extra_phys);
      storage.energies.crf_energy[topo.atom_energy_group(i)]
              [topo.atom_energy_group(i)] += 0.5 * e_crf;
      storage.energies.shift_extra_orig[topo.atom_energy_group(i)][topo.atom_energy_group(i)] += 0.5 * e_extra_orig;
      storage.energies.shift_extra_phys[topo.atom_energy_group(i)][topo.atom_energy_group(i)] += 0.5 * e_extra_phys;
      DEBUG(11, "\tcontribution " << 0.5 * e_crf);

      for (; it != to; ++it) {

        DEBUG(11, "\texcluded pair " << i << " - " << *it);

        periodicity.nearest_image(pos(i), pos(*it), r);
        DEBUG(10, "\tni i " << pos(i)(0) << " / "
                << pos(i)(1) << " / "
                << pos(i)(2));
        DEBUG(10, "\tni j " << pos(*it)(0) << " / "
                << pos(*it)(1) << " / "
                << pos(*it)(2));
        DEBUG(10, "\tni r " << r(0) << " / " << r(1) << " / " << r(2));
        shifted_rf_corr_interaction(r, charge_product(topo, i, *it), f, e_crf, e_extra_orig, e_extra_phys);

        force(i) += f;
        force(*it) -= f;

        for (int a = 0; a < 3; ++a)
          for (int b = 0; b < 3; ++b)
            storage.virial_tensor(a, b) += r(a) * f(b);
        DEBUG(11, "\tatomic virial done");

        // energy
        storage.energies.crf_energy[topo.atom_energy_group(i)]
                [topo.atom_energy_group(*it)] += e_crf;
        storage.energies.shift_extra_orig[topo.atom_energy_group(i)][topo.atom_energy_group(*it)] += e_extra_orig;
        storage.energies.shift_extra_phys[topo.atom_energy_group(i)][topo.atom_energy_group(*it)] += e_extra_phys;

        DEBUG(11, "\tcontribution " << e_crf);

      } // loop over excluded pairs
      break;
    }
    default:
      io::messages.add("Nonbonded_Innerloop",
              "interaction function not implemented",
              io::message::critical);
  }


}

template<typename t_nonbonded_spec>
inline void
interaction::Nonbonded_Innerloop<t_nonbonded_spec>::RF_solvent_interaction_innerloop
(
        topology::Topology & topo,
        configuration::Configuration & conf,
        topology::Chargegroup_Iterator const & cg_it,
        Storage & storage,
        math::Periodicity<t_nonbonded_spec::boundary_type> const & periodicity,
        unsigned int gamd
        ) {
  math::Vec r;
  double e_crf = 0.0;

  math::VArray &pos = conf.current().pos;

  // loop over the atoms
  topology::Atom_Iterator at_it = cg_it.begin(),
          at_to = cg_it.end();

  switch (t_nonbonded_spec::interaction_func) {
    case simulation::lj_crf_func:
    {
      for (; at_it != at_to; ++at_it) {
        DEBUG(11, "\tsolvent self term " << *at_it);
        // no solvent self term. The distance dependent part and the forces
        // are zero. The distance independent part should add up to zero
        // for the energies and is left out.

        for (topology::Atom_Iterator at2_it = at_it + 1; at2_it != at_to; ++at2_it) {
          //ORIOL_GAMD
          unsigned int igroup = 0;
          if (gamd){
            unsigned int gamdi = topo.gamd_accel_group(*at_it);
            unsigned int gamdj = topo.gamd_accel_group(*at2_it);
            std::vector<unsigned int> key = {gamdi, gamdj};
            igroup = topo.gamd_interaction_group(key);
          }
          DEBUG(11, "\tsolvent " << *at_it << " - " << *at2_it);
          periodicity.nearest_image(pos(*at_it),
                  pos(*at2_it), r);

          // for solvent, we don't calculate internal forces (rigid molecules)
          // and the distance independent parts should go to zero
          e_crf = -charge_product(topo, *at_it, *at2_it) *
                  math::four_pi_eps_i * crf_2cut3i() * abs2(r);
          DEBUG(15, "\tqi = " << topo.charge()(*at_it) << ", qj = " << topo.charge()(*at2_it));
          DEBUG(15, "\tcrf_2cut3i = " << crf_2cut3i() << ", abs2(r) = " << abs2(r));
          // energy
          storage.energies.crf_energy
                  [topo.atom_energy_group(*at_it) ]
                  [topo.atom_energy_group(*at2_it)] += e_crf;
          // ORIOL_GAMD
          if (gamd){
            storage.energies.gamd_potential_total[igroup] += e_crf;
          }
          DEBUG(11, "\tsolvent rf excluded contribution: " << e_crf);
        } // loop over at2_it
      } // loop over at_it
      break;
    }
    case simulation::cggromos_func :
    {
      if (topo.is_coarse_grained(*at_it)) { // CG-CG or CG-FG
        io::messages.add("Nonbonded_Innerloop",
              "no RF solvent interaction innerloop for coarse-grained simulations!",
              io::message::critical);
      } else { // FG-FG
        for (; at_it != at_to; ++at_it) {
          DEBUG(11, "\tsolvent self term " << *at_it);
          // no solvent self term. The distance dependent part and the forces
          // are zero. The distance independent part should add up to zero
          // for the energies and is left out.

          for (topology::Atom_Iterator at2_it = at_it + 1; at2_it != at_to; ++at2_it) {
            //ORIOL_GAMD
            unsigned int igroup = 0;
            if (gamd){
              unsigned int gamdi = topo.gamd_accel_group(*at_it);
              unsigned int gamdj = topo.gamd_accel_group(*at2_it);
              std::vector<unsigned int> key = {gamdi, gamdj};
              igroup = topo.gamd_interaction_group(key);
            }
            DEBUG(11, "\tsolvent " << *at_it << " - " << *at2_it);
            periodicity.nearest_image(pos(*at_it),
                    pos(*at2_it), r);

            // for solvent, we don't calculate internal forces (rigid molecules)
            // and the distance independent parts should go to zero
            e_crf = -topo.charge()(*at_it) * topo.charge()(*at2_it) *
                    math::four_pi_eps_i * crf_2cut3i() * abs2(r);
            DEBUG(15, "\tqi = " << topo.charge()(*at_it) << ", qj = " << topo.charge()(*at2_it));
            DEBUG(15, "\tcrf_2cut3i = " << crf_2cut3i() << ", abs2(r) = " << abs2(r));
            // energy
            storage.energies.crf_energy
                    [topo.atom_energy_group(*at_it) ]
                    [topo.atom_energy_group(*at2_it)] += e_crf;

            // ORIOL_GAMD
            if (gamd){
              storage.energies.gamd_potential_total[igroup] += e_crf;
            }
            DEBUG(11, "\tsolvent rf excluded contribution: " << e_crf);
          } // loop over at2_it
        } // loop over at_it
      }
      break;
    }
    case simulation::cgrain_func :
    {
      io::messages.add("Nonbonded_Innerloop",
              "no RF solvent interaction innerloop for coarse-grained simulations!",
              io::message::critical);
      break;
    }
    case simulation::pol_lj_crf_func:
    {
      math::Vec rp1, rp2, rpp;

      for (; at_it != at_to; ++at_it) {
        DEBUG(11, "\tsolvent self term " << *at_it);
        // no solvent self term. The distance dependent part and the forces
        // are zero. The distance independent part should add up to zero
        // for the energies and is left out.

        for (topology::Atom_Iterator at2_it = at_it + 1; at2_it != at_to; ++at2_it) {

          //ORIOL_GAMD
          unsigned int igroup = 0;
          if (gamd){
            unsigned int gamdi = topo.gamd_accel_group(*at_it);
            unsigned int gamdj = topo.gamd_accel_group(*at2_it);
            std::vector<unsigned int> key = {gamdi, gamdj};
            igroup = topo.gamd_interaction_group(key);
          }
          DEBUG(11, "\tsolvent " << *at_it << " - " << *at2_it);
          periodicity.nearest_image(pos(*at_it),
                  pos(*at2_it), r);
          rp1 = r - conf.current().posV(*at2_it);
          rp2 = r + conf.current().posV(*at_it);
          rpp = rp2 - conf.current().posV(*at2_it);

          // for solvent, we don't calculate internal forces (rigid molecules)
          // and the distance independent parts should go to zero
          const double cqi = topo.coscharge(*at_it);
          const double cqj = topo.coscharge(*at2_it);
          //const double qeps = -topo.charge()(*at_it) * topo.charge()(*at2_it) * math::four_pi_eps_i;
          //const double qepsp1 = -topo.charge()(*at_it)*(topo.charge()(*at2_it) - cqj) * math::four_pi_eps_i;
          //const double qepsp2 = (-topo.charge()(*at_it) - cqi) * topo.charge()(*at2_it) * math::four_pi_eps_i;
          //const double qepspp = cqi * cqj * math::four_pi_eps_i;
	  const double qeps = -(topo.charge()(*at_it)-cqi) * (topo.charge()(*at2_it)-cqj) * math::four_pi_eps_i;
	  const double qepsp1 = (-topo.charge()(*at_it) - cqi) * cqj * math::four_pi_eps_i;
	  const double qepsp2 = -cqi*(topo.charge()(*at2_it) - cqj) * math::four_pi_eps_i;
	  const double qepspp = cqi * cqj * math::four_pi_eps_i;

          e_crf = qeps * crf_2cut3i() * abs2(r) + qepsp1 * crf_2cut3i() * abs2(rp1)
                  + qepsp2 * crf_2cut3i() * abs2(rp2) + qepspp * crf_2cut3i() * abs2(rpp);

          // energy
          storage.energies.crf_energy
                  [topo.atom_energy_group(*at_it) ]
                  [topo.atom_energy_group(*at2_it)] += e_crf;

          // ORIOL_GAMD
          if (gamd){
            storage.energies.gamd_potential_total[igroup] += e_crf;
          }
        } // loop over at2_it
      } // loop over at_it

      break;
    }
    case simulation::pol_off_lj_crf_func: {
      math::Vec rp1, rp2, rpp;

      for ( ; at_it != at_to; ++at_it){
        DEBUG(11, "\tsolvent self term " << *at_it);
        // no solvent self term. The distance dependent part and the forces
        // are zero. The distance independent part should add up to zero
        // for the energies and is left out.
        math::Vec rim = math::Vec(0.0);
        math::Vec rm;
        if(topo.gamma(*at_it)!=0.0){
         math::Vec rij, rik;
         periodicity.nearest_image(conf.current().pos(*at_it),
                            conf.current().pos(topo.gamma_j(*at_it)), rij);
         periodicity.nearest_image(conf.current().pos(*at_it),
                            conf.current().pos(topo.gamma_k(*at_it)), rik);

         rim=topo.gamma(*at_it)*(rij+rik)/2;
        }
        for(topology::Atom_Iterator at2_it=at_it+1; at2_it!=at_to; ++at2_it){
          math::Vec rjm = math::Vec(0.0);
          //ORIOL_GAMD
          unsigned int igroup = 0;
          if (gamd){
            unsigned int gamdi = topo.gamd_accel_group(*at_it);
            unsigned int gamdj = topo.gamd_accel_group(*at2_it);
            std::vector<unsigned int> key = {gamdi, gamdj};
            igroup = topo.gamd_interaction_group(key);
          }
          if(topo.gamma(*at2_it)!=0.0){
           math::Vec rjj, rjk;
           periodicity.nearest_image(conf.current().pos(*at2_it),
                            conf.current().pos(topo.gamma_j(*at2_it)), rjj);
           periodicity.nearest_image(conf.current().pos(*at2_it),
                            conf.current().pos(topo.gamma_k(*at2_it)), rjk);
           rjm=topo.gamma(*at2_it)*(rjj+rjk)/2;
          }

          DEBUG(11, "\tsolvent " << *at_it << " - " << *at2_it);
          periodicity.nearest_image(pos(*at_it),
                  pos(*at2_it), r);
          rm=r+rjm-rim;
          rp1 = rm - conf.current().posV(*at2_it);
          rp2 = rm + conf.current().posV(*at_it);
          rpp = rp2 - conf.current().posV(*at2_it);

          // for solvent, we don't calculate internal forces (rigid molecules)
          // and the distance independent parts should go to zero
          const double cqi = topo.coscharge(*at_it);
          const double cqj = topo.coscharge(*at2_it);
          //const double qeps = -topo.charge()(*at_it) * topo.charge()(*at2_it) * math::four_pi_eps_i;
          //const double qepsp1 = -topo.charge()(*at_it)*(topo.charge()(*at2_it) - cqj) * math::four_pi_eps_i;
          //const double qepsp2 = (-topo.charge()(*at_it) - cqi) * topo.charge()(*at2_it) * math::four_pi_eps_i;
          //const double qepspp = cqi * cqj * math::four_pi_eps_i;
	  const double qeps = -(topo.charge()(*at_it)-cqi) * (topo.charge()(*at2_it)-cqj) * math::four_pi_eps_i;
	  const double qepsp1 = -(topo.charge()(*at_it) - cqi) * cqj * math::four_pi_eps_i;
	  const double qepsp2 = -cqi*(topo.charge()(*at2_it) - cqj) * math::four_pi_eps_i;
	  const double qepspp = cqi * cqj * math::four_pi_eps_i;
          e_crf = qeps * crf_2cut3i() * abs2(rm) + qepsp1 * crf_2cut3i() * abs2(rp1)
                  + qepsp2 * crf_2cut3i() * abs2(rp2) + qepspp * crf_2cut3i() * abs2(rpp);

          // energy
          storage.energies.crf_energy
                  [topo.atom_energy_group(*at_it) ]
                  [topo.atom_energy_group(*at2_it)] += e_crf;
          // ORIOL_GAMD
          if (gamd){
            storage.energies.gamd_potential_total[igroup] += e_crf;
          }
        } // loop over at2_it
      } // loop over at_it

      break;
    }
    case simulation::lj_shifted_crf_corr_func:
    {
      for (; at_it != at_to; ++at_it) {
        DEBUG(11, "\tsolvent self term " << *at_it);
        // no solvent self term. The distance dependent part and the forces
        // are zero. The distance independent part should add up to zero
        // for the energies and is left out.

        for (topology::Atom_Iterator at2_it = at_it + 1; at2_it != at_to; ++at2_it) {

          DEBUG(11, "\tsolvent " << *at_it << " - " << *at2_it);
          periodicity.nearest_image(pos(*at_it),
                  pos(*at2_it), r);

          const double dist2 = abs2(r);
          const double dist4 = dist2*dist2;
          const double dist6 = dist2*dist4;

          // for solvent, we don't calculate internal forces (rigid molecules)
          // and the distance independent parts should go to zero
          e_crf = -charge_product(topo, *at_it, *at2_it) *
                  math::four_pi_eps_i * (crf_2cut3i() * dist2 - a_RFm * dist4 - a_RFn * dist6);
          double e_extra_orig = -charge_product(topo, *at_it, *at2_it) * math::four_pi_eps_i * (a_RFm * dist4 + a_RFn * dist6);
          double e_extra_phys = -charge_product(topo, *at_it, *at2_it) * math::four_pi_eps_i * (a_RFm * dist4 + a_RFn * dist6);

          DEBUG(15, "\tqi = " << topo.charge()(*at_it) << ", qj = " << topo.charge()(*at2_it));
          DEBUG(15, "\tcrf_2cut3i = " << crf_2cut3i() << ", abs2(r) = " << abs2(r));
          DEBUG(15, "\a_RFm = " << a_RFm << " a_RFn = " << a_RFn);
          // energy
          storage.energies.crf_energy
                  [topo.atom_energy_group(*at_it) ]
                  [topo.atom_energy_group(*at2_it)] += e_crf;
          storage.energies.shift_extra_orig[topo.atom_energy_group(*at_it) ]
                  [topo.atom_energy_group(*at2_it)] += e_extra_orig;
          storage.energies.shift_extra_phys[topo.atom_energy_group(*at_it) ]
                  [topo.atom_energy_group(*at2_it)] += e_extra_phys;
          DEBUG(11, "\tsolvent rf excluded contribution: " << e_crf);
        } // loop over at2_it
      } // loop over at_it
      break;
    }
    default:
      io::messages.add("Nonbonded_Innerloop",
              "interaction function not implemented",
              io::message::critical);
  }

}

template<typename t_nonbonded_spec>
inline void
interaction::Nonbonded_Innerloop<t_nonbonded_spec>::electric_field_innerloop
(
        topology::Topology & topo,
        configuration::Configuration & conf,
        unsigned int i, unsigned int j, math::Vec &e_eli, math::Vec &e_elj,
        math::Periodicity<t_nonbonded_spec::boundary_type> const & periodicity
        ) {

  math::Vec r, rp1, rp2, e_el1, e_el2;

  // energy field term at position i and j
  DEBUG(11, "\tenergy field calculation i: " << i << " j: " << j);

  periodicity.nearest_image(conf.current().pos(i),
          conf.current().pos(j), r);
  DEBUG(10, "\tni i " << conf.current().pos(i)(0) << " / "
          << conf.current().pos(i)(1) << " / "
          << conf.current().pos(i)(2));
  DEBUG(10, "\tni j " << conf.current().pos(j)(0) << " / "
          << conf.current().pos(j)(1) << " / "
          << conf.current().pos(j)(2));
  DEBUG(10, "\tni r " << r(0) << " / " << r(1) << " / " << r(2));
  math::Vec rm = r;
  DEBUG(10, "\t topo.gamma(i) " << topo.gamma(i));
  if (topo.gamma(i) != 0.0 && simulation::pol_off_lj_crf_func != 0)  {
     math::Vec rij, rik, rim;
     periodicity.nearest_image(conf.current().pos(i),
        conf.current().pos(topo.gamma_j(i)), rij);
     periodicity.nearest_image(conf.current().pos(i),
        conf.current().pos(topo.gamma_k(i)), rik);
     rim = topo.gamma(i)*(rij + rik) / 2;
     rm -= rim;
  }
  DEBUG(10, "\t topo.gamma(j) " << topo.gamma(j));
  if (topo.gamma(j) != 0.0 && simulation::pol_off_lj_crf_func != 0) {
     math::Vec rjj, rjk, rjm;
     periodicity.nearest_image(conf.current().pos(j),
         conf.current().pos(topo.gamma_j(j)), rjj);
     periodicity.nearest_image(conf.current().pos(j),
         conf.current().pos(topo.gamma_k(j)), rjk);
     rjm = topo.gamma(j)*(rjj + rjk) / 2;
     rm += rjm;
  }
 

  switch(t_nonbonded_spec::efield_site) {
    case simulation::ef_atom : {
      rp1 = rm- conf.current().posV(j);
      rp2 = -rm- conf.current().posV(i);

      electric_field_interaction(rm, rp1, topo.charge(j),
              topo.coscharge(j), e_el1);
      electric_field_interaction(-rm, rp2, topo.charge(i),
              topo.coscharge(i), e_el2);
      break;
    }
    case simulation::ef_cos:
    {
      const math::Vec r_cos1 = conf.current().posV(i) + rm,
                      r_cos2 = conf.current().posV(j) - rm;
      rp1 = r_cos1 - conf.current().posV(j);
      rp2 = r_cos2 - conf.current().posV(i);

      electric_field_interaction(r_cos1, rp1, topo.charge(j),
              topo.coscharge(j), e_el1);
      electric_field_interaction(r_cos2, rp2, topo.charge(i),
              topo.coscharge(i), e_el2);
      break;
    }
    default:
      io::messages.add("electric_field_innerloop",
              "this site for electric field calculation was not implemented.",
              io::message::critical);

  }

  e_eli = e_el1;
  e_elj = e_el2;

}

template<typename t_nonbonded_spec>
inline void
interaction::Nonbonded_Innerloop<t_nonbonded_spec>::self_energy_innerloop
(
        topology::Topology & topo,
        configuration::Configuration & conf,
        unsigned int i,
        Storage & storage,
        math::Periodicity<t_nonbonded_spec::boundary_type> const & periodicity
        ) {
  DEBUG(8, "\tself energy of molecule i " << i);

  double self_e = 0.0;
  const double e_i2 = math::abs2(storage.electric_field(i));

  if (t_nonbonded_spec::pol_damping)
    self_energy_interaction(topo.polarisability(i), e_i2,
          topo.damping_level(i), topo.damping_power(i), self_e);
  else
    self_energy_interaction(topo.polarisability(i), e_i2, self_e);



  // self energy term
  DEBUG(10, "\tself energy storage:\t" << self_e);
  storage.energies.self_energy[topo.atom_energy_group(i)] += self_e;

}

template<typename t_nonbonded_spec>
inline void
interaction::Nonbonded_Innerloop<t_nonbonded_spec>::lj_ls_real_innerloop
(
        topology::Topology & topo,
        configuration::Configuration & conf,
        unsigned int i,
        unsigned int j,
        Storage & storage,
        math::Periodicity<t_nonbonded_spec::boundary_type> const & periodicity,
        unsigned int gamd
        ) {
  DEBUG(8, "\tpair\t" << i << "\t" << j);

  math::Vec r;
  double f = 0.0;
  double e_lj = 0.0, e_ls = 0.0;
  unsigned int igroup = 0;
  //ORIOL_GAMD
  if (gamd){
    unsigned int gamdi = topo.gamd_accel_group(i);
    unsigned int gamdj = topo.gamd_accel_group(j);
    std::vector<unsigned int> key = {gamdi, gamdj};
    igroup = topo.gamd_interaction_group(key);
  }
  periodicity.nearest_image(conf.current().pos(i),
          conf.current().pos(j), r);
  DEBUG(10, "\tni r " << r(0) << " / " << r(1) << " / " << r(2));

  switch (t_nonbonded_spec::interaction_func) {
    case simulation::lj_ls_func:
    {
      const lj_parameter_struct & lj =
              m_param->lj_parameter(topo.iac(i),
              topo.iac(j));

      DEBUG(11, "\tlj-parameter c6=" << lj.c6 << " c12=" << lj.c12);
      DEBUG(11, "\tcharge i=" << topo.charge()(i) << " j=" << topo.charge()(j));

      lj_ls_interaction(r, lj.c6, lj.c12,
              topo.charge(i) * topo.charge(j),
              f, e_lj, e_ls);
      DEBUG(12, "\t\t e_ls = " << e_ls << " f = " << f);
      DEBUG(10, "\t\tatomic virial");
      //ORIOL_GAMD
      for (int a = 0; a < 3; ++a) {
        const double term = f * r(a);
        storage.force(i)(a) += term;
        storage.force(j)(a) -= term;
        if (gamd){
          storage.force_gamd[igroup](i)(a) += term;
          storage.force_gamd[igroup](j)(a) -= term;
        }
        for (int b = 0; b < 3; ++b){
          storage.virial_tensor(b, a) += r(b) * term;
          if (gamd){
            storage.virial_tensor_gamd[igroup](b, a) += r(b) * term;
          }
        }
      }
      break;
    }

    default:
      io::messages.add("Nonbonded_Innerloop",
              "interaction function not implemented",
              io::message::critical);
  }

  // energy
  assert(storage.energies.lj_energy.size() >
          topo.atom_energy_group(i));
  assert(storage.energies.lj_energy.size() >
          topo.atom_energy_group(j));

  DEBUG(11, "\tenergy group i " << topo.atom_energy_group(i)
          << " j " << topo.atom_energy_group(j));

  storage.energies.lj_energy[topo.atom_energy_group(i)]
          [topo.atom_energy_group(j)] += e_lj;

  storage.energies.ls_real_energy[topo.atom_energy_group(i)]
          [topo.atom_energy_group(j)] += e_ls;
  if (gamd){
    storage.energies.gamd_potential_total[igroup] += e_lj +  e_ls;
  }
}

template<typename t_nonbonded_spec>
inline void
interaction::Nonbonded_Innerloop<t_nonbonded_spec>::ls_real_excluded_innerloop
(
        topology::Topology & topo,
        configuration::Configuration & conf,
        unsigned int i,
        unsigned int j,
        Storage & storage,
        math::Periodicity<t_nonbonded_spec::boundary_type> const & periodicity,
        unsigned int gamd
        ) {
  DEBUG(8, "\tpair\t" << i << "\t" << j);

  math::Vec r;
  double f = 0.0;
  double e_ls = 0.0;
  unsigned int igroup = 0;
  //ORIOL_GAMD
  if (gamd){
    unsigned int gamdi = topo.gamd_accel_group(i);
    unsigned int gamdj = topo.gamd_accel_group(j);
    std::vector<unsigned int> key = {gamdi, gamdj};
    igroup = topo.gamd_interaction_group(key);
  }
  periodicity.nearest_image(conf.current().pos(i),
          conf.current().pos(j), r);
  DEBUG(10, "\tni r " << r(0) << " / " << r(1) << " / " << r(2));

  switch (t_nonbonded_spec::interaction_func) {
    case simulation::lj_ls_func:
    {
      DEBUG(11, "\tcharge i=" << topo.charge()(i) << " j=" << topo.charge()(j));

      ls_excluded_interaction(r, topo.charge(i) * topo.charge(j),
              f, e_ls);
      DEBUG(10, "\t\t e_ls (exl. atoms) = " << e_ls);
      DEBUG(10, "\t\tatomic virial");
      //ORIOL_GAMD
      for (int a = 0; a < 3; ++a) {
        const double term = f * r(a);
        storage.force(i)(a) += term;
        storage.force(j)(a) -= term;
        if (gamd){
          storage.force_gamd[igroup](i)(a) += term;
          storage.force_gamd[igroup](j)(a) -= term; 
        storage.force_gamd[igroup](j)(a) -= term;
          storage.force_gamd[igroup](j)(a) -= term; 
        }

        for (int b = 0; b < 3; ++b){
          storage.virial_tensor(b, a) += r(b) * term;
          if (gamd){
            storage.virial_tensor_gamd[igroup](b, a) += r(b) * term;
          }
        }
      }
      break;
    }

    default:
      io::messages.add("Nonbonded_Innerloop",
              "interaction function not implemented",
              io::message::critical);
  }

  // energy
  assert(storage.energies.lj_energy.size() >
          topo.atom_energy_group(j));

  DEBUG(11, "\tenergy group i " << topo.atom_energy_group(i)
          << " j " << topo.atom_energy_group(j));

  storage.energies.ls_real_energy[topo.atom_energy_group(i)]
          [topo.atom_energy_group(j)] += e_ls;
  if (gamd){
    storage.energies.gamd_potential_total[igroup] += e_ls;
  }
}

/**
 * calculate the product of charges based on the charge type
 * this function implements variable charges for the QM buffer region
 */
template <typename t_nonbonded_spec>
double interaction::Nonbonded_Innerloop<t_nonbonded_spec>::charge_product(
        topology::Topology const & topo, 
        unsigned i, unsigned j) {
  double q = topo.charge(i) * topo.charge(j);
  switch (t_nonbonded_spec::charge_type) {
    case simulation::mm_charge : break;
    case simulation::qm_buffer_charge : {
      if (topo.is_adaptive_qm_buffer(i) != topo.is_adaptive_qm_buffer(j)) {
        DEBUG(11, "\tqm_delta_charge i=" << topo.qm_delta_charge(i) << " j=" << topo.qm_delta_charge(j));
        q +=  topo.charge(i) * topo.qm_delta_charge(j)
            + topo.charge(j) * topo.qm_delta_charge(i);
      }
      break;
    }
    default : io::messages.add("Charge type not implemented.", "nonbonded_innerloop", io::message::warning);
  }
  return q;
}
