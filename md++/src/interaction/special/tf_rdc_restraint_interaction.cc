/**
 * @file tf_rdc_restraint_interaction.cc
 * tensor-free RDC restraint interaction implementation
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

#include <interaction/special/tf_rdc_restraint_interaction.h>

#include <util/template_split.h>
#include <util/debug.h>

#include <io/topology/in_tf_rdc.h>
#include <list>

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE special

/**
 * tensor-free RDC restraint interactions
 */
template<math::boundary_enum B, math::virial_enum V>
int _calculate_tf_rdc_restraint_interactions
(topology::Topology & topo,
        configuration::Configuration & conf,
        simulation::Simulation & sim) {
    // loop over the tensor-free RDC restraints
    std::vector<topology::tf_rdc_restraint_struct>::iterator
    it = topo.tf_rdc_restraints().begin(),
          to = topo.tf_rdc_restraints().end();
    math::Periodicity<B> periodicity(conf.current().box);

    bool weighted =
          sim.param().tfrdc.mode == simulation::tfrdc_restr_av_weighted;

    double exptaur = 0.0;                                       // [-]
    double exptaut = 0.0;                                       // [-]
    if (sim.param().tfrdc.taur > 0.0)
      exptaur = exp(-sim.time_step_size() / sim.param().tfrdc.taur);    // [-]
    if (sim.param().tfrdc.taut > 0.0)
      exptaut = exp(-sim.time_step_size() / sim.param().tfrdc.taut);    // [-]
    const double & K = sim.param().tfrdc.K * 1e24;        // [kJ ps^2 / mol]

    // compute constant force terms
    const double dRavedR = 1.0 - exptaur;                       // [-]
    const double dPavedP = 1.0 - exptaut;                       // [-]
    const double q_scale = dRavedR/dPavedP;                     // [-]
    ++conf.special().tfrdc.num_averaged;
    for(unsigned int l = 0; it != to; ++it, ++l) {
        math::Vec dRdr(0.0, 0.0, 0.0), dPdr(0.0, 0.0, 0.0);
        math::Vec force(0.0, 0.0, 0.0);
        double & energy = conf.special().tfrdc.energy[l];       // [kJ / mol]
        double & R_avg = conf.special().tfrdc.R_avg[l];         // [-]
        double & P_avg = conf.special().tfrdc.P_avg[l];         // [-]

        const math::Vec & r_i = it->v1.pos(conf, topo);         // [nm]
        math::Vec r_ij;
        periodicity.nearest_image(r_i, it->v2.pos(conf, topo), r_ij);

        const math::Vec r_j(r_i - r_ij);  // vector r_ij points from j to i
        DEBUG(9, "r_i  :" << math::v2s(r_i));
        DEBUG(9, "r_j  :" << math::v2s(r_j));
        DEBUG(9, "r_ij :" << math::v2s(r_ij));

        const double d_r_ij = math::abs(r_ij);                  // [nm]
        const double d_r_ij_2 = d_r_ij * d_r_ij;                // [nm^2]
        DEBUG(9, "d_r_ij : " << d_r_ij);
        const double inv_r_ij = 1.0 / d_r_ij;                   // [1 / nm]
        const double inv_r_ij_2 = inv_r_ij * inv_r_ij;          // [1 / nm^2]
        const double inv_r_ij_3 = inv_r_ij_2 * inv_r_ij;        // [1 / nm^3]
        const double inv_r_ij_4 = inv_r_ij_2 * inv_r_ij_2;      // [1 / nm^4]
        const double inv_r_ij_5 = inv_r_ij_2 * inv_r_ij_3;      // [1 / nm^5]

        const double t_z_inv_r_ij_2 = 3.0 * r_ij(2) * r_ij(2) * inv_r_ij_2; // [-]

        const double r_0_3 = pow(it->normalisation_distance, 3.0); // [nm^3]

        // compute R and P
        // t_z_inv_r_ij_2 is 3*(z/r)^2
        // R and P therefore dimensionless
        const double R = r_0_3 * inv_r_ij_3;                    // [-]
        const double P = 0.5 * (t_z_inv_r_ij_2 - 1.0);          // [-]

        DEBUG(15, "R: " << R);
        DEBUG(15, "P: " << P);
        DEBUG(15, "read-in R_avg: " << R_avg);
        DEBUG(15, "read-in P_avg: " << P_avg);
        DEBUG(15, "Theta of k1k2: " << acos(r_ij(2)/d_r_ij));

        // v1.atom(0) is the first of four atoms defining one united atom.  It
        // is the atom sequence number of spin 1 minus 1 (which is the array
        // index of spin 1).

        // time-averaging
        if (sim.param().tfrdc.mode == simulation::tfrdc_restr_av ||
            sim.param().tfrdc.mode == simulation::tfrdc_restr_av_weighted) {
                // initialise average?
                if (sim.steps() == 0 && !sim.param().tfrdc.read) {
                    R_avg = R;                                  // [-]
                    P_avg = P;                                  // [-]
                }
                // apply time averaging
                R_avg = dRavedR * R + exptaur * R_avg;          // [-]
                P_avg = dPavedP * P + exptaut * P_avg;          // [-]
        }

        double term = 0.0;
        switch (sim.param().tfrdc.mode) {
            case simulation::tfrdc_restr_off:
                break;
            case simulation::tfrdc_restr_av:
            case simulation::tfrdc_restr_av_weighted:
            {
                // RDC_avg, D0, dD0 in 1/ps, times pow(10,12) back to Hz

                // compute tensor-free RDC
                double & RDC_avg = conf.special().tfrdc.RDC_avg[l]; // [1 / ps]

                DEBUG(9, "RDC_avg (before damping): " << RDC_avg*pow(10,12));
                RDC_avg = interaction::D_c(it) * R_avg * P_avg;     // [1 / ps]
                DEBUG(9, "RDC_avg (after damping): " << RDC_avg*pow(10,12));

                // the not averaged RDC
                conf.special().tfrdc.RDC[l] = interaction::D_c(it) * R * P; 

                // the cumulative average
                double & RDC_cumavg = conf.special().tfrdc.RDC_cumavg[l]; // [1 / ps]
                RDC_cumavg += (conf.special().tfrdc.RDC[l]-RDC_cumavg)/conf.special().tfrdc.num_averaged; 
                if (isnan(RDC_cumavg)) {
                  std::cerr << conf.special().tfrdc.num_averaged << " " << RDC_cumavg
                            << " P: " << P << " R: " << R << " z " << r_ij(2) 
                            << " t_z_inv_r_ij_2 " << t_z_inv_r_ij_2 << std::endl;
                }

                DEBUG(15, " R_avg: " << R_avg);
                DEBUG(15, " P_avg: " << P_avg);

                DEBUG(15, "RDC atom: " << it->v1.atom(0));
                DEBUG(15, "Input D0 [Hz]: " << it->D0*pow(10,12));
                DEBUG(15, "Input dD0 [Hz]: " << it->dD0*pow(10,12));
                DEBUG(15, "D_k^c [Hz]: " << interaction::D_c(it)*pow(10,12));

                if (RDC_avg > it->D0 + it->dD0) {
                  term = RDC_avg - it->D0 - it->dD0;            // [1 / ps]
                } else if (RDC_avg < it->D0 - it->dD0) {
                  term = RDC_avg - it->D0 + it->dD0;            // [1 / ps]
                } else {
                  DEBUG(9, "TFRDCRES  : restraint fulfilled");
                }

                dRdr = - 3.0 * r_0_3 * inv_r_ij_5 * r_ij;       // [1 / nm]
                for(unsigned int a = 0; a < 3; ++a){
                    if(a<2){
                        dPdr(a) = - t_z_inv_r_ij_2 * inv_r_ij_2 * r_ij(a); // [1 / nm]
                    } else {
                        dPdr(a) = - 3.0 * r_ij(2) * inv_r_ij_4;     // [1 / nm]
                        dPdr(a) *= (r_ij(2) * r_ij(2) - d_r_ij_2);  // [1 / nm]
                    }
                }
                break;

            }
            default:
              io::messages.add("Restraining method not implemented.",
                      "TF_RDC_Restraint_Interaction", io::message::error);
            return 1;
        }

        DEBUG(9, "dR/dr" << math::v2s(dRdr));
        DEBUG(9, "dP/dr" << math::v2s(dPdr));
        DEBUG(9, "dR_avg/dR: " << dRavedR);
        DEBUG(9, "dP_avg/dP: " << dPavedP);

        force = term * interaction::D_c(it) * (P_avg * dRavedR * dRdr + R_avg * dPavedP * dPdr);  // [1 / (ps^2 nm)]

        //TODO: remove the following line:
        math::Vec force_rescaled = term * interaction::D_c(it) * (P_avg * q_scale * dRdr + R_avg * dPdr);  // [1 / (ps^2 nm)]                  // [kJ / mol]
        DEBUG(9, "force           :" << math::v2s(force) );
        DEBUG(9, "force / rescaled:" << math::v2s(force_rescaled) );

        energy = 0.5 * K * term * term;                         // [kJ / mol]
        DEBUG(9, "Term in [1/ps]:" << term);
        DEBUG(9, "Energy before weighting [kJ/mol]:" << energy);
        math::Vec f_i(0.0, 0.0, 0.0), f_j(0.0, 0.0, 0.0);
        f_i = - K * force;                                      // [(kJ ps^2 / mol) * (1 / (ps^2 nm))] = [kJ / (mol nm)]
        f_j = -f_i;                                             // [kJ / (mol nm)]

        // weight if necessary
        if (weighted) {
          energy *= it->w;                                      // [kJ / mol]
          f_i *= it->w;                                         // [kJ / (mol nm)]
          f_j *= it->w;                                         // [kJ / (mol nm)]
        }
        std::cout.precision(15); // useful for debugging
        DEBUG(7, "energy: " << energy);
        DEBUG(7, "f_i: " << math::v2s(f_i));
        DEBUG(7, "f_j: " << math::v2s(f_j));

        // apply energy and force
        conf.current().energies.tfrdc_energy[topo.atom_energy_group()
    					    [it->v1.atom(0)]] += energy;
        it->v1.force(conf, topo, f_i);
        it->v2.force(conf, topo, f_j);


        for (int a = 0; a < 3; ++a) {
          for (int b = 0; b < 3; ++b) { 
            conf.current().virial_tensor(a, b) += r_ij(a) * f_i(b); 
          }
        } 
    }
    return 0;
}

int interaction::TF_RDC_Restraint_Interaction::calculate_interactions
(
        topology::Topology & topo,
        configuration::Configuration & conf,
        simulation::Simulation & sim) {
  m_timer.start(sim);
  SPLIT_VIRIAL_BOUNDARY(_calculate_tf_rdc_restraint_interactions,
          topo, conf, sim);
  m_timer.stop();
  return 0;
}

/**
 * init
*/
int interaction::TF_RDC_Restraint_Interaction::init
(
        topology::Topology &topo,
        configuration::Configuration &conf,
        simulation::Simulation &sim,
        std::ostream &os,
        bool quiet) {

  // resizing the containers in the configuration
  const unsigned int & num_res = topo.tf_rdc_restraints().size();
  conf.special().tfrdc.R_avg.resize(num_res);
  conf.special().tfrdc.P_avg.resize(num_res);
  conf.special().tfrdc.RDC.resize(num_res);
  conf.special().tfrdc.RDC_avg.resize(num_res);
  conf.special().tfrdc.RDC_cumavg.resize(num_res);
  conf.special().tfrdc.energy.resize(num_res);

  if (!sim.param().tfrdc.read) {
    conf.special().tfrdc.num_averaged=0;
  }

  if (!quiet) {
    os << "TENSOR-FREE RDC RESTRAINT INTERACTION" << std::endl;
    switch (sim.param().tfrdc.mode) {
      case simulation::tfrdc_restr_off:
        os << "\trestraining off";
        break;
      case simulation::tfrdc_restr_av:
        os << "\ttime-averaged restraining using exponential-decay memory function";
        break;
      case simulation::tfrdc_restr_av_weighted:
        os << "\ttime-averaged restraining, weighted, using exponential-decay memory function";
        break;
    }
    os << std::endl;

    os.precision(8);
    os << "  - Number of restraints: " << num_res << std::endl
            << "  - force constant: " << std::setw(15) << sim.param().tfrdc.K << std::endl;
    os << "  - time-averaging memory r relaxation time: " << std::setw(15) << sim.param().tfrdc.taur << std::endl;
    os << "  - time-averaging memory theta relaxation time: " << std::setw(15) << sim.param().tfrdc.taut << std::endl;
    if (sim.param().tfrdc.read)
      os << "  - reading initial averages from configuration." << std::endl;
    else
      os << "  - initialising averages to instantaneous value." << std::endl;
    os << std::endl << std::endl;

    os << "       N:     I    J    K    L  T1   I    J    K    L  T2   R0";
    os << "      G1        G2       D0          DD0        WTFRDC" << std::endl;
    os << "  -----------------------------------------------------------------";
    os << "-------------------------------------------------" << std::endl;
    // loop over the tensor-free RDC restraints
    std::vector<topology::tf_rdc_restraint_struct>::iterator
    it = topo.tf_rdc_restraints().begin(),
            to = topo.tf_rdc_restraints().end();
    for (int l = 1; it != to; ++it, ++l) {
      os << "  " << std::setw(6) << l << ": ";
      for (unsigned int i = 0; i < io::In_Tfrdcresspec::MAX_ATOMS; ++i)
        os << std::setw(5) << (int(i) < it->v1.size() ? it->v1.atom(i) + 1 : 0);
      os << std::setw(3) << it->v1.type();
      for (unsigned int i = 0; i < io::In_Tfrdcresspec::MAX_ATOMS; ++i)
        os << std::setw(5) << (int(i) < it->v2.size() ? it->v2.atom(i) + 1 : 0);
      os << std::setw(3) << it->v2.type();

      os.precision(4);
      os << std::setw(10) << it->normalisation_distance << std::setw(10) << it->gyri/0.010364272 << std::setw(9)
      << it->gyrj/0.010364272 << std::setw(12) << it->D0*1000000000000 << std::setw(11)
      << it->dD0*1000000000000 << std::setw(8) << it->w << std::endl;
    }
    os << "END" << std::endl;
  }

  return 0;
};
