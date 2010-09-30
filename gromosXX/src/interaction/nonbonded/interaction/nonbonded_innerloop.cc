/**
 * @file nonbonded_innerloop.cc
 * template methods of Nonbonded_Innerloop
 */

#include <vector>
#include <set>

#include "storage.h"
#include "nonbonded_innerloop.h"


#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE nonbonded

#include "solvent_innerloop.cc"

template<typename t_nonbonded_spec>
inline void 
interaction::Nonbonded_Innerloop<t_nonbonded_spec>::lj_crf_innerloop
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
  DEBUG(10, "\tni i " << conf.current().pos(i)(0) << " / " 
          <<  conf.current().pos(i)(1) << " / " 
          <<  conf.current().pos(i)(2));
    DEBUG(10, "\tni j " << conf.current().pos(j)(0) << " / " 
          <<  conf.current().pos(j)(1) << " / " 
          <<  conf.current().pos(j)(2));
  DEBUG(10, "\tni r " << r(0) << " / " << r(1) << " / " << r(2));
  
  switch(t_nonbonded_spec::interaction_func){
    case simulation::lj_crf_func :
      {
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
        for (int a=0; a<3; ++a){
          const double term = f * r(a);
          storage.force(i)(a) += term;
          storage.force(j)(a) -= term;
            
          for(int b=0; b<3; ++b)
            storage.virial_tensor(b, a) += r(b) * term;
        }
	break;
      }
    
    case simulation::cgrain_func :
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
        for (int a=0; a<3; ++a){
          const double term = f * r(a);
          storage.force(i)(a) += term;
          storage.force(j)(a) -= term;
            
          for(int b=0; b<3; ++b)
            storage.virial_tensor(b, a) += r(b) * term;
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
        for (int a = 0; a < 3; ++a) {
          const double term = f * r(a);
          storage.force(i)(a) += term;
          storage.force(j)(a) -= term;

          for (int b = 0; b < 3; ++b)
            storage.virial_tensor(b, a) += r(b) * term;
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
        for (int a = 0; a < 3; ++a) {
          const double term = f_pol[0] * r(a) + f_pol[1] * rp1(a) + f_pol[2] * rp2(a) + f_pol[3] * rpp(a);
          storage.force(i)(a) += term;
          storage.force(j)(a) -= term;

          for (int b = 0; b < 3; ++b)
            storage.virial_tensor(b, a) += r(b) * term;
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
        for (int a = 0; a < 3; ++a) {
          const double term = f * r(a);
          storage.force(i)(a) += term;
          storage.force(j)(a) -= term;

          for (int b = 0; b < 3; ++b)
            storage.virial_tensor(b, a) += r(b) * term;
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


          if(topo.gamma(i)!=0.0){
           storage.force(topo.gamma_j(i))(a) +=term*topo.gamma(i)/2;
           storage.force(topo.gamma_k(i))(a) +=term*topo.gamma(i)/2;
          }
          if(topo.gamma(j)!=0.0){
           storage.force(topo.gamma_j(j))(a) -=term*topo.gamma(j)/2;
           storage.force(topo.gamma_k(j))(a) -=term*topo.gamma(j)/2;
          }
          //bugbugbug, debug...
          for (int b = 0; b < 3; ++b)
        //    storage.virial_tensor(b, a) += (term+f_pol[0]*r(a))*r(b);
            storage.virial_tensor(b, a) += (term*rm(b))+(f_pol[0]*r(a))*r(b);
          }
        
      }
      break;
    }
    case simulation::lj_ls_func : {
      const lj_parameter_struct & lj =
      m_param->lj_parameter(topo.iac(i),
              topo.iac(j));
      
      DEBUG(11, "\tlj-parameter c6=" << lj.c6 << " c12=" << lj.c12);
      
      lj_interaction(r, lj.c6, lj.c12, f, e_lj);
      e_crf = 0.0;
      
      DEBUG(10, "\t\tatomic virial");
      for (int a=0; a<3; ++a){
        const double term = f * r(a);
        storage.force(i)(a) += term;
        storage.force(j)(a) -= term;
        
        for(int b=0; b<3; ++b)
          storage.virial_tensor(b, a) += r(b) * term;
      }
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
  
  storage.energies.lj_energy[topo.atom_energy_group(i)]
    [topo.atom_energy_group(j)] += e_lj;
  
  storage.energies.crf_energy[topo.atom_energy_group(i)]
    [topo.atom_energy_group(j)] += e_crf;
}

template<typename t_nonbonded_spec>
inline void
interaction::Nonbonded_Innerloop<t_nonbonded_spec>::sasa_innerloop
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 unsigned int i,
 Storage & storage,
 simulation::Simulation & sim,
 math::Periodicity<t_nonbonded_spec::boundary_type> const & periodicity)
{
  math::Vec rij;
  double e_sasa, e_sasa_temp;
  math::VArray &pos = conf.current().pos;

  const sasa_parameter_struct & sasa_i = m_param->sasa_parameter(i);

  DEBUG(1, "\tsasa-parameters: p_i=" << sasa_i.p_i << " r_i=" << sasa_i.r_i << " sigma_i=" << sasa_i.sigma_i);

  double surface = 4 * math::Pi * (sasa_i.r_i + sim.param().sasa.r_solv) *
      (sasa_i.r_i + sim.param().sasa.r_solv);

  // set initial energy to surface of atom i to avoid multiplication by zero
  e_sasa = surface;

  if (sasa_i.p_i != 0.0 && surface != 0.0) {
    // direct neighbours
    std::set< int >::const_iterator it, to;
    it = topo.sasa_first_neighbour(i).begin();
    to = topo.sasa_first_neighbour(i).end();
    for (; it != to; ++it) {
      DEBUG(11, "\tdirect neighbours " << i << " - " << *it);

      periodicity.nearest_image(pos(i), pos(*it), rij);

      DEBUG(11, "\tni r " << rij(0) << " / " << rij(1) << " / " << rij(2));

      double bij = sasa_overlap(topo, conf, i, *it, sim, periodicity);
      double pij = sim.param().sasa.p_12;

      if (bij != 0.0) {
        sasa_interaction(rij, bij, pij, sasa_i.p_i, surface, e_sasa_temp);
        e_sasa = e_sasa * e_sasa_temp;
      }
    }

    // second neighbours, 1,3
    it = topo.sasa_second_neighbour(i).begin();
    to = topo.sasa_second_neighbour(i).end();
    for (; it != to; ++it) {
      DEBUG(11, "\tsecond neighbours " << i << " - " << *it);

      periodicity.nearest_image(pos(i), pos(*it), rij);

      DEBUG(11, "\tni r " << rij(0) << " / " << rij(1) << " / " << rij(2));

      double bij = sasa_overlap(topo, conf, i, *it, sim, periodicity);
      double pij = sim.param().sasa.p_13;

      if (bij != 0.0) {
        sasa_interaction(rij, bij, pij, sasa_i.p_i, surface, e_sasa_temp);
        e_sasa = e_sasa * e_sasa_temp;
      }
    }

    // third neighbours, 1,4
    it = topo.sasa_third_neighbour(i).begin();
    to = topo.sasa_third_neighbour(i).end();
    for (; it != to; ++it) {
      DEBUG(11, "\tthird neighbours " << i << " - " << *it);

      periodicity.nearest_image(pos(i), pos(*it), rij);

      DEBUG(11, "\tni r " << rij(0) << " / " << rij(1) << " / " << rij(2));

      double bij = sasa_overlap(topo, conf, i, *it, sim, periodicity);
      double pij = sim.param().sasa.p_13; // same as for second neighbours

      if (bij != 0.0) {
        sasa_interaction(rij, bij, pij, sasa_i.p_i, surface, e_sasa_temp);
        e_sasa = e_sasa * e_sasa_temp;
      }
    }

    // higher neighbours, > 1,4, this will be excluded...
    //if (sim.param().sasa.switch_1x) {
      it = topo.sasa_higher_neighbour(i).begin();
      to = topo.sasa_higher_neighbour(i).end();
      for (; it != to; ++it) {
        DEBUG(11, "\thigher neighbours " << i << " - " << *it);

        periodicity.nearest_image(pos(i), pos(*it), rij);

        DEBUG(10, "\tni r " << rij(0) << " / " << rij(1) << " / " << rij(2));

        double bij = sasa_overlap(topo, conf, i, *it, sim, periodicity);
        double pij = sim.param().sasa.p_1x;

        if (bij != 0.0) {
          sasa_interaction(rij, bij, pij, sasa_i.p_i, surface, e_sasa_temp);
          e_sasa = e_sasa * e_sasa_temp;
        }
      }
    //}
    // sasa
    conf.current().sasa_area[i] = e_sasa;
    conf.current().sasa_tot += e_sasa;
    DEBUG(1, "\tcurrent total sasa: " << conf.current().sasa_tot <<
        "\tsasa of atom i: " << conf.current().sasa_area[i] << " - step:" << i);

    // now calculate the sasa energy
    e_sasa *= sasa_i.sigma_i;

    conf.current().energies.sasa_energy[topo.atom_energy_group(i)] += e_sasa;
  }
}

template<typename t_nonbonded_spec>
inline void
interaction::Nonbonded_Innerloop<t_nonbonded_spec>::sasa_innerloop_force
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 unsigned int i, double amax,
 Storage & storage,
 simulation::Simulation & sim,
 math::Periodicity<t_nonbonded_spec::boundary_type> const & periodicity)
{
  math::Vec rij;
  double fsasa, fvolume, ai;
  math::VArray &pos = conf.current().pos;
  math::VArray &force = conf.current().force;
  //math::VArray &fs = conf.current().fsasa; // only needed for testing
  //math::VArray &volume = conf.current().fvolume;// only needed for testing

  double const a_max = amax * sim.param().sasa.max_cut;
  double const a_min = a_max - sim.param().sasa.min_cut;
  ai = conf.current().sasa_area[i];

  const unsigned int end = topo.num_solute_atoms();

  for (unsigned int j = 0; j < end; ++j) {
    if (i != j) {
      const sasa_parameter_struct & sasa_j = m_param->sasa_parameter(j);
      DEBUG(11, "\tsasa-parameter p_j=" << sasa_j.p_i);

      double surface = 4 * math::Pi * (sasa_j.r_i + sim.param().sasa.r_solv) *
          (sasa_j.r_i + sim.param().sasa.r_solv);
      double a, aj;
      a = conf.current().sasa_area[j];

      double pji;
      double bji = sasa_overlap(topo, conf, j, i, sim, periodicity);
      double dbji = sasa_overlap_der(topo, conf, j, i, sim, periodicity);

      if (dbji != 0.0 && sasa_j.p_i != 0.0 && surface != 0.0) {

        DEBUG(11, "\tatom pair " << i << " - " << j);

        periodicity.nearest_image(pos(j), pos(i), rij);

        if (topo.sasa_first_neighbour(j).count(i)) pji = sim.param().sasa.p_12;
        else if (topo.sasa_second_neighbour(j).count(i)) pji = sim.param().sasa.p_13;
        else if (topo.sasa_third_neighbour(j).count(i)) pji = sim.param().sasa.p_13;
        else if (topo.sasa_higher_neighbour(j).count(i)) pji = sim.param().sasa.p_1x;

        aj = a / (1 - sasa_j.p_i * pji * bji / surface);

        // force for sasa interaction
        fsasa = -sasa_j.p_i * pji * dbji * aj / surface;

        // force for volume interaction
        if (sim.param().sasa.switch_volume) {
          double dg = sasa_switching_fct_der(topo, conf, a, a_max, a_min, storage, sim, periodicity);
          fvolume = 4 / 3 * math::Pi * dg * fsasa * sasa_j.r_i * sasa_j.r_i * sasa_j.r_i;
        }

        for (int a = 0; a < 3; ++a) {
          // stores each component of sasa contribution
          const double term_s = sasa_j.sigma_i * fsasa * rij(a) / math::abs(rij);
          force(i)(a) += term_s;
          //fs(i)(a) += term_s;

          // stores each component of volume contribution
          if (sim.param().sasa.switch_volume) {
            const double term_v = sim.param().sasa.sigma_v * fvolume * rij(a) / math::abs(rij);
            force(i)(a) += term_v;
            //volume(i)(a) += term_v;
          }
        }
      } // store forces
    } // for i != j
  } // loop over all atoms

  // forces for j = i
  const sasa_parameter_struct & sasa_i = m_param->sasa_parameter(i);
  DEBUG(11, "\tsasa-parameter p_i=" << sasa_i.p_i);

  double surface_i = 4 * math::Pi * (sasa_i.r_i + sim.param().sasa.r_solv) *
      (sasa_i.r_i + sim.param().sasa.r_solv);

  // first neighbours
  std::set< int >::const_iterator it, to;
  it = topo.sasa_first_neighbour(i).begin();
  to = topo.sasa_first_neighbour(i).end();
  for (; it != to; ++it) {
    DEBUG(11, "\tdirect neighbours " << i << " - " << *it);

    double bij = sasa_overlap(topo, conf, i, *it, sim, periodicity);
    double dbij = sasa_overlap_der(topo, conf, i, *it, sim, periodicity);
    double pij = sim.param().sasa.p_12;

    if (dbij != 0.0 && sasa_i.p_i != 0.0 && surface_i != 0.0) {

      DEBUG(11, "\tatom pair " << i << " - " << *it);

      periodicity.nearest_image(pos(i), pos(*it), rij);

      double aij = ai / (1 - sasa_i.p_i * pij * bij / surface_i);

      // sasa contribution for i=j
      fsasa = -(sasa_i.p_i * pij * dbij * aij / surface_i);
      if (sim.param().sasa.switch_volume){
        double dg = sasa_switching_fct_der(topo, conf, ai, a_max, a_min, storage, sim, periodicity);
        fvolume = 4 / 3 * math::Pi * dg * sasa_i.r_i * sasa_i.r_i * sasa_i.r_i * fsasa;
      }

      for (int a = 0; a < 3; ++a) {
        // sasa contribution for i=j
        const double term_s = sasa_i.sigma_i * fsasa * rij(a) / math::abs(rij);
        force(i)(a) -= term_s;
        //fs(i)(a) -= term_s;

        // volume contribution for i=j
        if (sim.param().sasa.switch_volume) {
          const double term_v = sim.param().sasa.sigma_v * fvolume * rij(a) / math::abs(rij);
          force(i)(a) -= term_v;
          //volume(i)(a) -= term_v;
        }
      }
    }
  }

  // second neighbours, 1,3
  it = topo.sasa_second_neighbour(i).begin();
  to = topo.sasa_second_neighbour(i).end();
  for (; it != to; ++it) {
    DEBUG(11, "\tsecond neighbours " << i << " - " << *it);

    double bij = sasa_overlap(topo, conf, i, *it, sim, periodicity);
    double dbij = sasa_overlap_der(topo, conf, i, *it, sim, periodicity);
    double pij = sim.param().sasa.p_13;

    if (dbij != 0.0 && sasa_i.p_i != 0.0 && surface_i != 0.0) {

      DEBUG(11, "\tatom pair " << i << " - " << *it);

      periodicity.nearest_image(pos(i), pos(*it), rij);

      double aij = ai / (1 - sasa_i.p_i * pij * bij / surface_i);

      // sasa contribution for i=j
      fsasa = -(sasa_i.p_i * pij * dbij * aij / surface_i);
      if (sim.param().sasa.switch_volume){
        double dg = sasa_switching_fct_der(topo, conf, ai, a_max, a_min, storage, sim, periodicity);
        fvolume = 4 / 3 * math::Pi * dg * sasa_i.r_i * sasa_i.r_i * sasa_i.r_i * fsasa;
      }

      for (int a = 0; a < 3; ++a) {
        // sasa contribution for i=j
        const double term_s = sasa_i.sigma_i * fsasa * rij(a) / math::abs(rij);
        force(i)(a) -= term_s;
        //fs(i)(a) -= term_s;

        // volume contribution for i=j
        if (sim.param().sasa.switch_volume){
          const double term_v = sim.param().sasa.sigma_v * fvolume * rij(a) / math::abs(rij);
          force(i)(a) -= term_v;
          //volume(i)(a) -= term_v;
        }
      }
    }
  }

  // third neighbours, 1,4
  it = topo.sasa_third_neighbour(i).begin();
  to = topo.sasa_third_neighbour(i).end();
  for (; it != to; ++it) {
    DEBUG(11, "\tthird neighbours " << i << " - " << *it);

    double bij = sasa_overlap(topo, conf, i, *it, sim, periodicity);
    double dbij = sasa_overlap_der(topo, conf, i, *it, sim, periodicity);
    double pij = sim.param().sasa.p_13;

    if (dbij != 0.0 && sasa_i.p_i != 0.0 && surface_i != 0.0) {

      DEBUG(11, "\tatom pair " << i << " - " << *it);

      periodicity.nearest_image(pos(i), pos(*it), rij);

      double aij = ai / (1 - sasa_i.p_i * pij * bij / surface_i);

      // sasa contribution for i=j
      fsasa = -(sasa_i.p_i * pij * dbij * aij / surface_i);

      if (sim.param().sasa.switch_volume){
        double dg = sasa_switching_fct_der(topo, conf, ai, a_max, a_min, storage, sim, periodicity);
        fvolume = 4 / 3 * math::Pi * dg * sasa_i.r_i * sasa_i.r_i * sasa_i.r_i * fsasa;
      }

      for (int a = 0; a < 3; ++a) {
        // sasa contribution for i=j
        const double term_s = sasa_i.sigma_i * fsasa * rij(a) / math::abs(rij);
        force(i)(a) -= term_s;
        //fs(i)(a) -= term_s;

        // volume contribution for i=j
        if (sim.param().sasa.switch_volume) {
          const double term_v = sim.param().sasa.sigma_v * fvolume * rij(a) / math::abs(rij);
          force(i)(a) -= term_v;
          //volume(i)(a) -= term_v;
        }
      }
    }
  }

  // higher neighbours, > 1,4
  //if (sim.param().sasa.switch_1x) {
    it = topo.sasa_higher_neighbour(i).begin();
    to = topo.sasa_higher_neighbour(i).end();
    for (; it != to; ++it) {
      DEBUG(11, "\thigher neighbours " << i << " - " << *it);

      double bij = sasa_overlap(topo, conf, i, *it, sim, periodicity);
      double dbij = sasa_overlap_der(topo, conf, i, *it, sim, periodicity);
      double pij = sim.param().sasa.p_1x;

      if (dbij != 0.0 && sasa_i.p_i != 0.0 && surface_i != 0.0) {

        DEBUG(11, "\tatom pair " << i << " - " << *it);

        periodicity.nearest_image(pos(i), pos(*it), rij);

        double aij = ai / (1 - sasa_i.p_i * pij * bij / surface_i);
        fsasa = -(sasa_i.p_i * pij * dbij * aij / surface_i);

        // sasa contribution for i=j
        if (sim.param().sasa.switch_volume) {
          double dg = sasa_switching_fct_der(topo, conf, ai, a_max, a_min, storage, sim, periodicity);
          fvolume = 4 / 3 * math::Pi * dg * sasa_i.r_i * sasa_i.r_i * sasa_i.r_i * fsasa;
        }

        for (int a = 0; a < 3; ++a) {
          // sasa contribution for i=j
          const double term_s = sasa_i.sigma_i * fsasa * rij(a) / math::abs(rij);
          force(i)(a) -= term_s;
          //fs(i)(a) -= term_s;

          // volume contribution for i=j
          if (sim.param().sasa.switch_volume) {
            const double term_v = sim.param().sasa.sigma_v * fvolume * rij(a) / math::abs(rij);
            force(i)(a) -= term_v;
            //volume(i)(a) -= term_v;
          }
        }
      }
    }
  //}
}

template<typename t_nonbonded_spec>
inline void
interaction::Nonbonded_Innerloop<t_nonbonded_spec>::sasa_volume_innerloop
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 unsigned int i, double amax,
 Storage & storage,
 simulation::Simulation & sim,
 math::Periodicity<t_nonbonded_spec::boundary_type> const & periodicity)
{
  const sasa_parameter_struct & sasa_i =
      m_param->sasa_parameter(i);
  double a_min, a_max, a, r, evolume, volume;

  a_max = amax * sim.param().sasa.max_cut;
  a_min = a_max - sim.param().sasa.min_cut;
  assert(a_min <= a_max);

  a = conf.current().sasa_area[i];
  r = sasa_i.r_i;

  double g = sasa_switching_fct(topo, conf, a, a_max, a_min, storage, sim, periodicity);
  // radius of hydrogen atoms is set to -r_solv, hydrogen must be excluded to avoid negative values
  // for the ineterior of a protein...
  if (r > 0) {
    volume = 4 / 3 * g * math::Pi * r * r * r;
    evolume = sim.param().sasa.sigma_v * volume;

    conf.current().sasa_vol[i] = volume;
    conf.current().sasavol_tot += volume;
    conf.current().energies.sasa_volume_energy[topo.atom_energy_group(i)] += evolume;
  }
}

template<typename t_nonbonded_spec>
inline double
interaction::Nonbonded_Innerloop<t_nonbonded_spec>::sasa_overlap
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 unsigned int i, unsigned int j,
 simulation::Simulation & sim,
 math::Periodicity<t_nonbonded_spec::boundary_type> const & periodicity)
{
  double b;
  math::Vec r;
  const sasa_parameter_struct & sasa_i = m_param->sasa_parameter(i);
  const sasa_parameter_struct & sasa_j = m_param->sasa_parameter(j);

  periodicity.nearest_image(conf.current().pos(i),
      conf.current().pos(j), r);
  DEBUG(10, "\tni r " << r(0) << " / " << r(1) << " / " << r(2));

  double absr = math::abs(r);
  double term = sasa_i.r_i + sasa_j.r_i + 2 * sim.param().sasa.r_solv;
  double diff = sasa_j.r_i - sasa_i.r_i;

  if (sasa_j.r_i == (-sim.param().sasa.r_solv)) {
    DEBUG(1, "\texcluded atom " << j << " - radius " << sasa_j.r_i);
    return 0.0;
  }
  else if (absr < term) {
    b = math::Pi * (sasa_i.r_i + sim.param().sasa.r_solv) * (term - absr) * (1 + diff/absr);
    return b;
  }
  return 0.0;
}

template<typename t_nonbonded_spec>
inline double
interaction::Nonbonded_Innerloop<t_nonbonded_spec>::sasa_overlap_der
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 unsigned int i, unsigned int j,
 simulation::Simulation & sim,
 math::Periodicity<t_nonbonded_spec::boundary_type> const & periodicity)
{
  double b;
  math::Vec r;
  const sasa_parameter_struct & sasa_i = m_param->sasa_parameter(i);
  const sasa_parameter_struct & sasa_j = m_param->sasa_parameter(j);

  periodicity.nearest_image(conf.current().pos(i),
      conf.current().pos(j), r);
  DEBUG(10, "\tni r " << r(0) << " / " << r(1) << " / " << r(2));

  double absr = math::abs(r);
  double absr2 = absr * absr;
  double term = sasa_i.r_i + sasa_j.r_i + 2 * sim.param().sasa.r_solv;
  double diff = sasa_j.r_i - sasa_i.r_i;

  if (sasa_j.r_i == (-sim.param().sasa.r_solv)) {
    DEBUG(1, "\texcluded atom " << j << " - radius " << sasa_j.r_i);
    return 0.0;
  }
  else if (absr < term && absr > 0.0) {
    b = (-1) * math::Pi * (sasa_i.r_i + sim.param().sasa.r_solv) * (1 + diff/absr) +
        math::Pi * (sasa_i.r_i + sim.param().sasa.r_solv) * (term - absr) * (- diff/absr2);
    return b;
  }
  return 0.0;
}


template<typename t_nonbonded_spec>
inline double
interaction::Nonbonded_Innerloop<t_nonbonded_spec>::sasa_switching_fct
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 double a, double amax, double amin,
 Storage & storage,
 simulation::Simulation & sim,
 math::Periodicity<t_nonbonded_spec::boundary_type> const & periodicity)
{
  DEBUG(1, "sasa surface of atom i: " << a << "\tmax sasa: " << amax << "\tmin sasa: " << amin);
  double g, diff, da;
  assert(a >= 0);
  assert(amax >= 0);
  assert(amin >= 0);
  if (a <= amin && a >= 0) return 1.0;
  else if (a <= amax && a > amin) {
    diff = amax - amin;
    da = a - amin;
    g = 2*(da * da * da)/(diff * diff * diff) - 3*(da * da)/(diff * diff) + 1;
    return g;
  }
  else return 0.0;
}

template<typename t_nonbonded_spec>
inline double
interaction::Nonbonded_Innerloop<t_nonbonded_spec>::sasa_switching_fct_der
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 double a, double amax, double amin,
 Storage & storage,
 simulation::Simulation & sim,
 math::Periodicity<t_nonbonded_spec::boundary_type> const & periodicity)
{
  double g, diff, da;
  assert(a >= 0);
  assert(amax >= 0);
  assert(amin >= 0);
  if (a <= amin && a >= 0) return 0.0;
  else if (a <= amax && a > amin) {
    diff = amax - amin;
    da = a - amin;
    g = 6*(da * da)/(diff * diff * diff) - 6*(da)/(diff * diff);
    return g;
  }
  else return 0.0;
}

template<typename t_nonbonded_spec>
void interaction::Nonbonded_Innerloop<t_nonbonded_spec>::one_four_interaction_innerloop
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 int i,
 int j,
 Storage & storage,
 math::Periodicity<t_nonbonded_spec::boundary_type> const & periodicity)
{
  DEBUG(8, "\t1,4-pair\t" << i << "\t" << j);
  
  math::Vec r;
  double f, e_lj, e_crf = 0.0, e_ls = 0.0;
  
  periodicity.nearest_image(conf.current().pos(i), 
			    conf.current().pos(j), r);
  DEBUG(10, "\tni i " << conf.current().pos(i)(0) << " / " 
          <<  conf.current().pos(i)(1) << " / " 
          <<  conf.current().pos(i)(2));
    DEBUG(10, "\tni j " << conf.current().pos(j)(0) << " / " 
          <<  conf.current().pos(j)(1) << " / " 
          <<  conf.current().pos(j)(2));
  DEBUG(10, "\tni r " << r(0) << " / " << r(1) << " / " << r(2));
  switch(t_nonbonded_spec::interaction_func){
    case simulation::lj_crf_func :
      {
	const lj_parameter_struct & lj = 
	  m_param->lj_parameter(topo.iac(i),
				topo.iac(j));
	
	DEBUG(11, "\tlj-parameter cs6=" << lj.cs6 << " cs12=" << lj.cs12);
	DEBUG(11, "\tcharge i=" << topo.charge()(i) << " j=" << topo.charge()(j));
	
        lj_crf_interaction(r, lj.cs6, lj.cs12,
			 topo.charge()(i) * 
			 topo.charge()(j),
			 f, e_lj, e_crf);
        
        DEBUG(10, "\t\tatomic virial");
        for (int a=0; a<3; ++a){
          const double term = f * r(a);
          storage.force(i)(a) += term;
          storage.force(j)(a) -= term;
            
          for(int b=0; b<3; ++b)
            storage.virial_tensor(b, a) += r(b) * term;
        }
	break;
      }
    case simulation::cgrain_func :
      {
        io::messages.add("Nonbonded_Innerloop",
                         "no 1,4 interactions for coarse-grained simulations!",
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
          
          for(int b=0; b<3; ++b)
            storage.virial_tensor(b, a) += r(b)*term;
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

          if(topo.gamma(i)!=0.0){
           storage.force(topo.gamma_j(i))(a) +=topo.gamma(i)/2*term;
           storage.force(topo.gamma_k(i))(a) -=topo.gamma(j)/2*term;
          }
          if(topo.gamma(j)!=0.0){
           storage.force(topo.gamma_j(j))(a) +=topo.gamma(i)/2*term;
           storage.force(topo.gamma_k(j))(a) -=topo.gamma(j)/2*term;
          }
          for(int b=0; b<3; ++b)
           // storage.virial_tensor(b,a ) += (term+f_pol[0]*r(a))*r(b);
            storage.virial_tensor(b, a) += (term*rm(b))+(f_pol[0]*r(a))*r(b);
        }
    std::cout << "just a bit of debugging 2\n";
        break;
    }
    case simulation::lj_ls_func : {
      const lj_parameter_struct & lj =
      m_param->lj_parameter(topo.iac(i),
              topo.iac(j));
      
      DEBUG(11, "\tlj-parameter c6=" << lj.c6 << " c12=" << lj.c12);
      lj_interaction(r, lj.cs6, lj.cs12, f, e_lj);
      e_crf = 0.0;
           
      DEBUG(10, "\t\tatomic virial");
      for (int a=0; a<3; ++a){
        const double term = f * r(a);
        storage.force(i)(a) += term;
        storage.force(j)(a) -= term;
        
        for(int b=0; b<3; ++b)
          storage.virial_tensor(b, a) += r(b) * term;
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
 math::Periodicity<t_nonbonded_spec::boundary_type> const & periodicity)
{
  DEBUG(8, "\tLJ exception\t" << ljex.i << "\t" << ljex.j);

  math::Vec r;
  double f, e_lj, e_crf = 0.0, e_ls = 0.0;

  unsigned int i = ljex.i, j = ljex.j;

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
        for (int a=0; a<3; ++a){
          const double term = f * r(a);
          storage.force(i)(a) += term;
          storage.force(j)(a) -= term;

          for(int b=0; b<3; ++b)
            storage.virial_tensor(b, a) += r(b) * term;
        }
	break;
      }

    case simulation::cgrain_func :
      {
        io::messages.add("Nonbonded_Innerloop",
                         "LJ exceptions for coarse-grained simulations!",
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
        for (int a=0; a<3; ++a){
          const double term = f_pol[0]*r(a) + f_pol[1]*rp1(a) + f_pol[2]*rp2(a) + f_pol[3]*rpp(a);
          storage.force(i)(a) += term;
          storage.force(j)(a) -= term;

          for(int b=0; b<3; ++b)
            storage.virial_tensor(b, a) += r(b)*term;
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

          if(topo.gamma(i)!=0.0){
           storage.force(topo.gamma_j(i))(a) +=topo.gamma(i)/2*term;
           storage.force(topo.gamma_k(i))(a) -=topo.gamma(j)/2*term;
          }
          if(topo.gamma(j)!=0.0){
           storage.force(topo.gamma_j(j))(a) +=topo.gamma(i)/2*term;
           storage.force(topo.gamma_k(j))(a) -=topo.gamma(j)/2*term;
          }
          for(int b=0; b<3; ++b)
          // storage.virial_tensor(b,a ) += (term+f_pol[0]*r(a))*r(b);
           storage.virial_tensor(b, a) += (term*rm(b))+(f_pol[0]*r(a))*r(b);
        }
   std::cout << "just a bit of debugging 3\n";
        break;
    }
    case simulation::lj_ls_func : {
      DEBUG(11, "\tlj-parameter c6=" << ljex.c6 << " c12=" << ljex.c12);
      lj_interaction(r, ljex.c6, ljex.c12, f, e_lj);
      e_crf = 0.0;

      DEBUG(10, "\t\tatomic virial");
      for (int a=0; a<3; ++a){
        const double term = f * r(a);
        storage.force(i)(a) += term;
        storage.force(j)(a) -= term;

        for(int b=0; b<3; ++b)
          storage.virial_tensor(b, a) += r(b) * term;
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
 math::Periodicity<t_nonbonded_spec::boundary_type> const & periodicity)
{
  math::Vec r, f;
  double e_crf;
  
  math::VArray &pos   = conf.current().pos;
  math::VArray &force = storage.force;

  std::set<int>::const_iterator it, to;
  it = topo.exclusion(i).begin();
  to = topo.exclusion(i).end();
  
  DEBUG(8, "\tself-term " << i );
  r=0;
  
  switch(t_nonbonded_spec::interaction_func){
    case simulation::lj_crf_func :
    {
      // this will only contribute in the energy, the force should be zero.
      rf_interaction(r, topo.charge()(i) * topo.charge()(i), f, e_crf);
      storage.energies.crf_energy[topo.atom_energy_group(i)]
              [topo.atom_energy_group(i)] += 0.5 * e_crf;
      DEBUG(11, "\tcontribution " << 0.5*e_crf);
      
      for( ; it != to; ++it){
        
        DEBUG(11, "\texcluded pair " << i << " - " << *it);
        
        periodicity.nearest_image(pos(i), pos(*it), r);
          DEBUG(10, "\tni i " << pos(i)(0) << " / " 
          <<  pos(i)(1) << " / " 
          <<  pos(i)(2));
    DEBUG(10, "\tni j " << pos(*it)(0) << " / " 
          <<  pos(*it)(1) << " / " 
          <<  pos(*it)(2));
  DEBUG(10, "\tni r " << r(0) << " / " << r(1) << " / " << r(2));
        rf_interaction(r, topo.charge()(i) * topo.charge()(*it), f, e_crf);
        
        force(i) += f;
        force(*it) -= f;
        
        for(int a=0; a<3; ++a)
          for(int b=0; b<3; ++b)
            storage.virial_tensor(a, b) += r(a) * f(b);
        DEBUG(11, "\tatomic virial done");
             
        // energy
        storage.energies.crf_energy[topo.atom_energy_group(i)]
                [topo.atom_energy_group(*it)] += e_crf;
        
        DEBUG(11, "\tcontribution " << e_crf);
        
      } // loop over excluded pairs
      break;
    }
    case simulation::cgrain_func :{
      io::messages.add("Nonbonded_Innerloop",
                       "no RF excluded interactions for coarse-grained simulations!",
		       io::message::critical);
      break;
    }
    case simulation::pol_lj_crf_func : {
      math::Vec rp1, rp2, rpp;
      math::VArray f_pol(4);
      f_pol = 0.0;
  
      DEBUG(8, "\tself-term " << i );
      rp1 = -conf.current().posV(*it);
      rp2 = conf.current().posV(i);
      rpp = 0.0;
      
      // this will only contribute in the energy, the force should be zero.
      pol_rf_interaction(r, rp1, rp2, rpp, topo.charge(i), topo.charge(i),
              topo.coscharge(i), topo.coscharge(i), f_pol, e_crf);
      storage.energies.crf_energy[topo.atom_energy_group(i)]
              [topo.atom_energy_group(i)] += 0.5 * e_crf;
      DEBUG(11, "\tcontribution " << 0.5*e_crf);
      
      for( ; it != to; ++it){
        
        DEBUG(11, "\texcluded pair " << i << " - " << *it);
        
        periodicity.nearest_image(pos(i), pos(*it), r);
        rp1 = r - conf.current().posV(*it);
        rp2 = r + conf.current().posV(i);
        rpp = rp2 - conf.current().posV(*it);
        
        pol_rf_interaction(r, rp1, rp2, rpp, topo.charge()(i),
                topo.charge()(*it), topo.coscharge(i),
                topo.coscharge(*it), f_pol, e_crf);

        for (int a = 0; a < 3; ++a) {
          const double term = f_pol(0)(a) + f_pol(1)(a) + f_pol(2)(a) + f_pol(3)(a);
          force(i)(a) += term;
          force(*it)(a) -= term;
          for (int b = 0; b < 3; ++b)
            storage.virial_tensor(b, a) += r(b) * term;
        }
        DEBUG(11, "\tatomic virial done");

        // energy
        storage.energies.crf_energy[topo.atom_energy_group(i)]
                [topo.atom_energy_group(*it)] += e_crf;
        
        DEBUG(11, "\tcontribution " << e_crf);
        
      } // loop over excluded pairs
      break;
    }
    case simulation::pol_off_lj_crf_func : {
      math::Vec rp1, rp2, rpp;
      math::VArray f_pol(4);
      f_pol = 0.0;

      DEBUG(8, "\tself-term " << i );
      rp1 = -conf.current().posV(*it);
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
      DEBUG(11, "\tcontribution " << 0.5*e_crf);

      for( ; it != to; ++it){
        periodicity.nearest_image(pos(i), pos(*it), r);
        if(topo.gamma(*it)!=0.0){
          math::Vec rjj, rjk;
          periodicity.nearest_image(conf.current().pos(*it),
                            conf.current().pos(topo.gamma_j(*it)), rjj);
          periodicity.nearest_image(conf.current().pos(*it),
                            conf.current().pos(topo.gamma_k(*it)), rjk);
          rm+=topo.gamma(*it)*(rjj+rjk)/2;
        }
        rm-=rim;
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
          if(topo.gamma(i)!=0.0){
           force(topo.gamma_j(i))(a) +=topo.gamma(i)/2*term;
           force(topo.gamma_k(i))(a) +=topo.gamma(i)/2*term;
          }
          if(topo.gamma(*it)!=0.0){
           force(topo.gamma_j(*it))(a) -=topo.gamma(*it)/2*term;
           force(topo.gamma_k(*it))(a) -=topo.gamma(*it)/2*term;
          }
          for(int b=0; b<3; ++b)
            storage.virial_tensor(b, a) += rm(b)*term; 
        }
        DEBUG(11, "\tatomic virial done");

        // energy
        storage.energies.crf_energy[topo.atom_energy_group(i)]
                [topo.atom_energy_group(*it)] += e_crf;

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
 math::Periodicity<t_nonbonded_spec::boundary_type> const & periodicity
 )
{
  math::Vec r;
  double e_crf;
  
  math::VArray &pos   = conf.current().pos;

  // loop over the atoms
  topology::Atom_Iterator at_it = cg_it.begin(),
    at_to = cg_it.end();
  
  switch(t_nonbonded_spec::interaction_func){
    case simulation::lj_crf_func :
    {
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
          DEBUG(15,"\tqi = " << topo.charge()(*at_it) << ", qj = " << topo.charge()(*at2_it));
          DEBUG(15,"\tcrf_2cut3i = " << crf_2cut3i() << ", abs2(r) = " << abs2(r));
          // energy
          storage.energies.crf_energy
                  [topo.atom_energy_group(*at_it) ]
                  [topo.atom_energy_group(*at2_it)] += e_crf;
          DEBUG(11,"\tsolvent rf excluded contribution: " << e_crf);
        } // loop over at2_it
      } // loop over at_it
      break;
    }
    case simulation::cgrain_func :{
      io::messages.add("Nonbonded_Innerloop",
                       "no RF solvent interaction innerloop for coarse-grained simulations!",
		       io::message::critical);
      break;
    }
    case simulation::pol_lj_crf_func: {
      math::Vec rp1, rp2, rpp;
      
      for ( ; at_it != at_to; ++at_it){
        DEBUG(11, "\tsolvent self term " << *at_it);
        // no solvent self term. The distance dependent part and the forces
        // are zero. The distance independent part should add up to zero
        // for the energies and is left out.
        
        for(topology::Atom_Iterator at2_it=at_it+1; at2_it!=at_to; ++at2_it){
          
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
          const double qeps = -topo.charge()(*at_it) * topo.charge()(*at2_it) * math::four_pi_eps_i;
          const double qepsp1 = -topo.charge()(*at_it)*(topo.charge()(*at2_it) - cqj) * math::four_pi_eps_i;
          const double qepsp2 = (-topo.charge()(*at_it) - cqi) * topo.charge()(*at2_it) * math::four_pi_eps_i;
          const double qepspp = cqi * cqj * math::four_pi_eps_i;
          e_crf = qeps * crf_2cut3i() * abs2(r) + qepsp1 * crf_2cut3i() * abs2(rp1)
                  + qepsp2 * crf_2cut3i() * abs2(rp2) + qepspp * crf_2cut3i() * abs2(rpp);

          // energy
          storage.energies.crf_energy
                  [topo.atom_energy_group(*at_it) ]
                  [topo.atom_energy_group(*at2_it)] += e_crf;
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
          const double qeps = -topo.charge()(*at_it) * topo.charge()(*at2_it) * math::four_pi_eps_i;
          const double qepsp1 = -topo.charge()(*at_it)*(topo.charge()(*at2_it) - cqj) * math::four_pi_eps_i;
          const double qepsp2 = (-topo.charge()(*at_it) - cqi) * topo.charge()(*at2_it) * math::four_pi_eps_i;
          const double qepspp = cqi * cqj * math::four_pi_eps_i;
          e_crf = qeps * crf_2cut3i() * abs2(r) + qepsp1 * crf_2cut3i() * abs2(rp1)
                  + qepsp2 * crf_2cut3i() * abs2(rp2) + qepspp * crf_2cut3i() * abs2(rpp);

          // energy
          storage.energies.crf_energy
                  [topo.atom_energy_group(*at_it) ]
                  [topo.atom_energy_group(*at2_it)] += e_crf;
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
 )
{

  math::Vec r, rp1, rp2, e_el1, e_el2;

  // energy field term at position i and j
  DEBUG(11, "\tenergy field calculation i: "<<i<<" j: "<<j);

  periodicity.nearest_image(conf.current().pos(i),
                            conf.current().pos(j), r);
    DEBUG(10, "\tni i " << conf.current().pos(i)(0) << " / " 
          <<  conf.current().pos(i)(1) << " / " 
          <<  conf.current().pos(i)(2));
    DEBUG(10, "\tni j " << conf.current().pos(j)(0) << " / " 
          <<  conf.current().pos(j)(1) << " / " 
          <<  conf.current().pos(j)(2));
  DEBUG(10, "\tni r " << r(0) << " / " << r(1) << " / " << r(2));
  math::Vec rm = r;
  DEBUG(10, "\t topo.gamma(i) " << topo.gamma(i));
  if (topo.gamma(i) != 0.0 && simulation::pol_off_lj_crf_func)  {
     math::Vec rij, rik, rim;
     periodicity.nearest_image(conf.current().pos(i),
        conf.current().pos(topo.gamma_j(i)), rij);
     periodicity.nearest_image(conf.current().pos(i),
        conf.current().pos(topo.gamma_k(i)), rik);
     rim = topo.gamma(i)*(rij + rik) / 2;
     rm -= rim;
  }
  DEBUG(10, "\t topo.gamma(j) " << topo.gamma(j));
  if (topo.gamma(j) != 0.0 && simulation::pol_off_lj_crf_func) {
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
    case simulation::ef_cos : {
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
    default :
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
 )
{
  DEBUG(8, "\tself energy of molecule i " << i);

  double self_e;
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
 math::Periodicity<t_nonbonded_spec::boundary_type> const & periodicity
 )
{
  DEBUG(8, "\tpair\t" << i << "\t" << j);
  
  math::Vec r;
  double f;
  double e_lj, e_ls;
  
  periodicity.nearest_image(conf.current().pos(i), 
			    conf.current().pos(j), r);
  DEBUG(10, "\tni r " << r(0) << " / " << r(1) << " / " << r(2));
  
  switch(t_nonbonded_spec::interaction_func){
    case simulation::lj_ls_func : {
      const lj_parameter_struct & lj =
      m_param->lj_parameter(topo.iac(i),
              topo.iac(j));
      
      DEBUG(11, "\tlj-parameter c6=" << lj.c6 << " c12=" << lj.c12);
      DEBUG(11, "\tcharge i=" << topo.charge()(i) << " j=" << topo.charge()(j));
      
      lj_ls_interaction(r, lj.c6, lj.c12,
                        topo.charge(i) * topo.charge(j), 
                        f, e_lj, e_ls);
      DEBUG(12,"\t\t e_ls = " << e_ls << " f = " << f);
      DEBUG(10, "\t\tatomic virial");
      for (int a=0; a<3; ++a){
        const double term = f * r(a);
        storage.force(i)(a) += term;
        storage.force(j)(a) -= term;
        
        for(int b=0; b<3; ++b)
          storage.virial_tensor(b, a) += r(b) * term;
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
 math::Periodicity<t_nonbonded_spec::boundary_type> const & periodicity
 )
{
  DEBUG(8, "\tpair\t" << i << "\t" << j);
  
  math::Vec r;
  double f;
  double e_ls;
  
  periodicity.nearest_image(conf.current().pos(i), 
			    conf.current().pos(j), r);
  DEBUG(10, "\tni r " << r(0) << " / " << r(1) << " / " << r(2));
  
  switch(t_nonbonded_spec::interaction_func){
    case simulation::lj_ls_func : {
      DEBUG(11, "\tcharge i=" << topo.charge()(i) << " j=" << topo.charge()(j));
      
      ls_excluded_interaction(r, topo.charge(i) * topo.charge(j), 
                              f, e_ls);
      DEBUG(10,"\t\t e_ls (exl. atoms) = " << e_ls);
      DEBUG(10, "\t\tatomic virial");
      for (int a=0; a<3; ++a){
        const double term = f * r(a);
        storage.force(i)(a) += term;
        storage.force(j)(a) -= term;
        
        for(int b=0; b<3; ++b)
          storage.virial_tensor(b, a) += r(b) * term;
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
}
