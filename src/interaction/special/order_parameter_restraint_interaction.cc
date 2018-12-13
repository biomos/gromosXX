/**
 * @file order_parameter_restraint_interaction.cc
 * order parameter restraint interaction implementation
 */
#include <stdheader.h>

#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <simulation/simulation.h>
#include <configuration/configuration.h>
#include <interaction/interaction.h>
#include <interaction/forcefield/forcefield.h>

#include <math/periodicity.h>

// special interactions
#include <interaction/interaction_types.h>

#include <interaction/special/order_parameter_restraint_interaction.h>

#include <util/template_split.h>
#include <util/debug.h>

#include <io/topology/in_order.h>
#include <list>

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE special


/**
 * calculate order parameter restraint interactions
 */
template<math::boundary_enum B, math::virial_enum V>
int _calculate_order_parameter_restraint_interactions
(topology::Topology & topo,
        configuration::Configuration & conf,
        simulation::Simulation & sim) {
  
  // loop over the order parameter restraints
  std::vector<topology::order_parameter_restraint_struct>::iterator
  it = topo.order_parameter_restraints().begin(),
          to = topo.order_parameter_restraints().end();

  math::Periodicity<B> periodicity(conf.current().box);

  bool weighted =
          sim.param().orderparamrest.orderparamrest == simulation::oparam_restr_av_weighted;
  
  double exptau = 0.0;
  if (sim.param().orderparamrest.tau > 0.0)
    exptau = exp(-sim.time_step_size() / sim.param().orderparamrest.tau);
  

  for(unsigned int l = 0; it != to; ++it, ++l) {
    // get the two positions and their connections
    // they should be gathered!
    const math::Vec & r_i = it->v1.pos(conf, topo);
    math::Vec r_ij;
    periodicity.nearest_image(r_i, it->v2.pos(conf, topo), r_ij);
    //const math::Vec r_j(r_i + r_ij);
    const math::Vec r_j(r_i - r_ij);  // vector r_ij points from j to i !!!

    DEBUG(9, "r_i  :" << math::v2s(r_i));
    DEBUG(9, "r_j  :" << math::v2s(r_j));
    DEBUG(9, "r_ij :" << math::v2s(r_ij));

    const double d_r_ij = math::abs(r_ij);
    //const double d_r_ij_2 = d_r_ij_2 * d_r_ij_2; // THIS WAS THE ACTUAL BUG IN THE CODE
    const double d_r_ij_2 = d_r_ij * d_r_ij;   
    DEBUG(9, "d_r_ij : " << d_r_ij);
    const double inv_r_ij = 1.0 / d_r_ij;
    const double inv_r_ij_2 = inv_r_ij * inv_r_ij;
    const double inv_r_ij_3 = inv_r_ij_2 * inv_r_ij;
    const double inv_r_ij_5 = inv_r_ij_2 * inv_r_ij_3;
    const double inv_r_ij_7 = inv_r_ij_2 * inv_r_ij_5;  
    
    const double r_eff_6 = pow(it->normalisation_distance, 6.0);
    
    // get Q, dQ/dr and dD/dr
    math::Matrix Q;
    math::GenericMatrix<math::Vec> dQdr;
    math::Vec dDdr;
    for(unsigned int a = 0; a < 3; ++a) {
      const double ria_m_rja = r_i(a) - r_j(a);
      dDdr(a) = -3.0 * ria_m_rja * inv_r_ij_5;
      for(unsigned int b = 0; b < 3; ++b) {
        const double rib_m_rjb = r_i(b) - r_j(b);
        Q(a,b) = ria_m_rja * rib_m_rjb * inv_r_ij_5;
        for(unsigned int g = 0; g < 3; ++g) {
          const double rig_m_rjg = r_i(g) - r_j(g);
          double term = 0.0;
          if (g == a) term += rib_m_rjb;
          if (g == b) term += ria_m_rja;
          dQdr(a,b)(g) = (d_r_ij_2 * term - 5.0 * ria_m_rja * rib_m_rjb * rig_m_rjg) * inv_r_ij_7;
        }
        DEBUG(12, "dQdr(" << a << "," << b << "): " << math::v2s(dQdr(a,b)));
      }
    }
    DEBUG(10, "Q:\n" << math::m2s(Q));
    math::Matrix & Q_avg = conf.special().orderparamres.Q_avg[l];
    
    // get D
    const double D = inv_r_ij_3;
    double & D_avg = conf.special().orderparamres.D_avg[l];
    DEBUG(10, "D: " << D);
    DEBUG(10, "dDdr: " << math::v2s(dDdr));

    // time-averaging
    if (sim.param().orderparamrest.orderparamrest == simulation::oparam_restr_av ||
        sim.param().orderparamrest.orderparamrest == simulation::oparam_restr_av_weighted) {
      // initalise average?
      if (sim.steps() == 0 && !sim.param().orderparamrest.read) {
        Q_avg = Q;
        D_avg = D;
      }
      // apply time averaging
      Q_avg = (1.0 - exptau) * Q + exptau * Q_avg;
      D_avg = (1.0 - exptau) * D + exptau * D_avg;
    } else if (sim.param().orderparamrest.orderparamrest == simulation::oparam_restr_winav ||
        sim.param().orderparamrest.orderparamrest == simulation::oparam_restr_winav_weighted) {
      unsigned int window_size = int(sim.param().orderparamrest.tau / sim.time_step_size()) /
              sim.param().orderparamrest.update_step;
      DEBUG(8, "window size: " << window_size);
      if (sim.steps() == 0 && !sim.param().orderparamrest.read) {
        for(unsigned int i = 0; i < window_size; ++i) {
          conf.special().orderparamres.Q_winavg[l].push_back(Q);
          conf.special().orderparamres.D_winavg[l].push_back(D);
        }
      }

      if (sim.steps() == 0 || (sim.steps() && sim.steps() % sim.param().orderparamrest.update_step == 0)) {
        Q_avg = 0.0;
        conf.special().orderparamres.Q_winavg[l].pop_front();
        conf.special().orderparamres.Q_winavg[l].push_back(Q);
        for (std::list<math::Matrix>::const_iterator it = conf.special().orderparamres.Q_winavg[l].begin(),
                to = conf.special().orderparamres.Q_winavg[l].end(); it != to; ++it) {
          Q_avg += *it;
        }
        Q_avg *= (1.0 / window_size);

        D_avg = 0.0;
        conf.special().orderparamres.D_winavg[l].pop_front();
        conf.special().orderparamres.D_winavg[l].push_back(D);
        for (std::list<double>::const_iterator it = conf.special().orderparamres.D_winavg[l].begin(),
                to = conf.special().orderparamres.D_winavg[l].end(); it != to; ++it) {
          D_avg += *it;
        }
        D_avg *= (1.0 / window_size);
      }
    }
    

    // compute order parameter
    double & S2_avg = conf.special().orderparamres.S2_avg[l];

    double sum = 0.0;
    math::Vec sum_force(0.0, 0.0, 0.0);
    for(unsigned int a = 0; a < 3; ++a) {
      for(unsigned int b = 0; b < 3; ++b) {
        sum += Q_avg(a,b) * Q_avg(a,b);
        sum_force += Q_avg(a,b) * dQdr(a,b);
      }
    }
    S2_avg = 0.5 * (3.0 * sum - D_avg * D_avg) * r_eff_6;
    DEBUG(8, " S2_avg: " << S2_avg);
    
    DEBUG(10, "sum_force: " << math::v2s(sum_force));

    // compute the energy
    double & energy = conf.special().orderparamres.energy[l];
    const double & K = sim.param().orderparamrest.K;
    math::Vec f_i(0.0, 0.0, 0.0), f_j(0.0, 0.0, 0.0);
    switch (sim.param().orderparamrest.orderparamrest) {
      case simulation::oparam_restr_off:
        break;
      case simulation::oparam_restr_av:
      case simulation::oparam_restr_av_weighted:
      case simulation::oparam_restr_winav:
      case simulation::oparam_restr_winav_weighted:
      {
        double term = 0.0;
        if (S2_avg > it->S0 + it->dS0) {
          term = S2_avg - it->S0 - it->dS0; 
        } else if (S2_avg < it->S0 - it->dS0) {
          term = S2_avg - it->S0 + it->dS0;
        }
        
        energy = 0.5 * K * term * term;
        //double force_term = 0.5 * K * term; // why not -K * term
        double force_term = -K * term;  // for the minus sign we have to use the correct definition of r_ij, the 0.5 was wrong, see eq. (55)
        f_i = (force_term * r_eff_6) * (3.0 * sum_force - D_avg * dDdr);
        f_j = -f_i;
        break;
      }
      default:
        io::messages.add("Restraining method not implemented.",
                "Order_Parameter_Restraint_Interaction", io::message::error);
        return 1;
    }

    // weight if necessary
    if (weighted) {
      energy *= it->w;
      f_i *= it->w;
      f_j *= it->w;
    }
    std::cout.precision(15); // useful for debugging
    DEBUG(7, "energy: " << energy);
    DEBUG(7, "f_i: " << math::v2s(f_i));
    DEBUG(7, "f_j: " << math::v2s(f_j));
    
    // apply energy and force
    conf.current().energies.oparam_total += energy;
    it->v1.force(conf, topo, f_i); 
    it->v2.force(conf, topo, f_j); 
  }

  return 0;
}

int interaction::Order_Parameter_Restraint_Interaction::calculate_interactions
(
        topology::Topology & topo,
        configuration::Configuration & conf,
        simulation::Simulation & sim) {
  m_timer.start();
  SPLIT_VIRIAL_BOUNDARY(_calculate_order_parameter_restraint_interactions,
          topo, conf, sim);
  m_timer.stop();
  return 0;
}

/**
 * init
 */
int interaction::Order_Parameter_Restraint_Interaction::init
(
        topology::Topology &topo,
        configuration::Configuration &conf,
        simulation::Simulation &sim,
        std::ostream &os,
        bool quiet) {

  // resizing the containers in the configuration
  const unsigned int & num_res = topo.order_parameter_restraints().size();
  conf.special().orderparamres.Q_avg.resize(num_res);
  conf.special().orderparamres.D_avg.resize(num_res);
  conf.special().orderparamres.Q_winavg.resize(num_res);
  conf.special().orderparamres.D_winavg.resize(num_res);
  conf.special().orderparamres.S2_avg.resize(num_res);
  conf.special().orderparamres.energy.resize(num_res);

  if (!quiet) {
    os << "ORDER-PARAMETER RESTRAINT INTERACTION" << std::endl;
    switch (sim.param().orderparamrest.orderparamrest) {
      case simulation::oparam_restr_off:
        os << "\trestraining off";
        break;
      case simulation::oparam_restr_av:
        os << "\ttime-averaged restraining using exponential-decay memory function";
        break;
      case simulation::oparam_restr_av_weighted:
        os << "\ttime-averaged restraining, weighted, using exponential-decay memory function";
        break;
      case simulation::oparam_restr_winav:
        os << "\ttime-averaged restraining using window averaging";
        break;
      case simulation::oparam_restr_winav_weighted:
        os << "\ttime-averaged restraining, weighted, using window averaging";
        break;
    }
    os << std::endl;

    os.precision(8);
    os << "  - Number of restraints: " << num_res << std::endl
            << "  - force constant: " << std::setw(15) << sim.param().orderparamrest.K << std::endl;
    os << "  - time-averaging memory relaxation time: " << std::setw(15) << sim.param().orderparamrest.tau << std::endl;
    os << "  - updating average every  " << sim.param().orderparamrest.update_step << " step." << std::endl;
    if (sim.param().orderparamrest.read)
      os << "  - reading initial averages from configuration." << std::endl;
    else
      os << "  - initialising averages to instantaneous value." << std::endl;
    os << std::endl << std::endl;
    
    os << "       N:     I    J    K    L  T1   I    J    K    L  T2 S0      DS0     WOPR" << std::endl;
    os << "  --------------------------------------------------------------------------------" << std::endl;
    // loop over the order parameter restraints
    std::vector<topology::order_parameter_restraint_struct>::iterator
    it = topo.order_parameter_restraints().begin(),
            to = topo.order_parameter_restraints().end();
    for (int l = 1; it != to; ++it, ++l) {
      os << "  " << std::setw(6) << l << ": ";
      for (unsigned int i = 0; i < io::In_Orderparamresspec::MAX_ATOMS; ++i)
        os << std::setw(5) << (int(i) < it->v1.size() ? it->v1.atom(i) + 1 : 0);
      os << std::setw(3) << it->v1.type();
      for (unsigned int i = 0; i < io::In_Orderparamresspec::MAX_ATOMS; ++i)
        os << std::setw(5) << (int(i) < it->v2.size() ? it->v2.atom(i) + 1 : 0);
      os << std::setw(3) << it->v2.type();

      os.precision(4);
      os << std::setw(8) << it->S0 << std::setw(8) << it->dS0 << std::setw(8) << it->w << std::endl;
    }
    os << "END" << std::endl;
  }

  return 0;
};
