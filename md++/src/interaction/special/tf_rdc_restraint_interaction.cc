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
#include <io/topology/in_tf_rdc.h>

#include <util/template_split.h>
#include "../../util/error.h"
#include <util/debug.h>
#include <util/generate_velocities.h>

#include <list>

#include "../../math/random.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE special

template<typename T>
std::vector<double> _linspace(T start_in, T end_in, int num_in)
{

  std::vector<double> linspaced;

  double start = static_cast<double>(start_in);
  double end = static_cast<double>(end_in);
  double num = static_cast<double>(num_in);

  if (num == 0) { return linspaced; }
  if (num == 1) 
    {
      linspaced.push_back(start);
      return linspaced;
    }

  double delta = (end - start) / (num - 1);

  for(int i=0; i < num-1; ++i)
    {
      linspaced.push_back(start + delta * i);
    }
  linspaced.push_back(end); // I want to ensure that start and end
                            // are exactly the same as the input
  return linspaced;
}

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

        double P;
        if (sim.param().tfrdc.nstsd > 0) {
          P = conf.special().tfrdc_mfv.P_avg[l];
        } else {
          P = 0.5 * (t_z_inv_r_ij_2 - 1.0);   // [-]
        }       

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
                    // P_avg = P;                                  // [-]
                    P_avg = it->D0 / interaction::D_c(it); // set initial value so RDC_avg will be approximately the target value to dampen forces in the first steps
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
                if (std::isnan(RDC_cumavg)) {
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
                if (sim.param().tfrdc.nstsd > 0) {
                  dPdr=conf.special().tfrdc_mfv.dPdr_avg[l];
                } else {
                  for(unsigned int a = 0; a < 3; ++a){
                      if(a<2){
                          dPdr(a) = - t_z_inv_r_ij_2 * inv_r_ij_2 * r_ij(a); // [1 / nm]
                      } else {
                          dPdr(a) = - 3.0 * r_ij(2) * inv_r_ij_4;     // [1 / nm]
                          dPdr(a) *= (r_ij(2) * r_ij(2) - d_r_ij_2);  // [1 / nm]
                      }
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


template<math::boundary_enum B>
int _shake_bond(math::Box box, math::VArray &pos, const math::VArray &old_pos, double mass, double d) {
  math::Vec r;
  math::Periodicity<B> periodicity(box);
  periodicity.nearest_image(pos(0), pos(1), r);
  double dist2 = abs2(r);
  double constr_length2 = d  * d;
  double diff = constr_length2-dist2;
  // hardcoded tolerance 1e-4
  int iter=0;
  while (fabs(diff) >= constr_length2 * 1e-4 * 2.0) {
    DEBUG(10, "iteration number " << iter);
    math::Vec ref_r;
    periodicity.nearest_image(old_pos(0), old_pos(1), ref_r);

    double sp = dot(ref_r, r);
    if (sp < constr_length2 * math::epsilon) {
      /*
      io::messages.add("SHAKE error. vectors orthogonal",
          "Shake::???",
          io::message::critical);
      */
      std::cout << "SHAKE ERROR tfrdcres magnetic field vector\n"
              << "\tref i     : " << math::v2s(old_pos(0)) << "\n"
              << "\tref j     : " << math::v2s(old_pos(1)) << "\n"
              << "\tfree i    : " << math::v2s(pos(0)) << "\n"
              << "\tfree j    : " << math::v2s(pos(1)) << "\n"
              << "\tref r     : " << math::v2s(ref_r) << "\n"
              << "\tr         : " << math::v2s(r) << "\n"
              << "\tsp        : " << sp << "\n"
              << "\tconstr    : " << constr_length2 << "\n"
              << "\tdiff      : " << diff << "\n";
      io::messages.add("SHAKE error tfrdcres: vectors orthogonal",
              "Shake::solute",
              io::message::error);

      return E_SHAKE_FAILURE;
    }

    if (iter > 1000) {
      std::cout << "SHAKE error tfrdcres: too many iterations (>1000)\n";
      return E_SHAKE_FAILURE;
    }
    // lagrange multiplier
    double lambda = diff / (sp * 2.0 * 2/mass);

    DEBUG(10, "lagrange multiplier " << lambda);

    const math::Vec cons_force = lambda * ref_r;
    // update positions
    ref_r *= lambda;
    pos(0) += ref_r * 1/mass;
    pos(1) -= ref_r * 1/mass;

    periodicity.nearest_image(pos(0), pos(1), r);
    dist2 = abs2(r);
    diff = constr_length2-dist2;
    iter++;
  } // shake iteration
  return 0;
}

void _add_to_theta_distribution(const math::Vec axis, configuration::Configuration & conf,
        simulation::Simulation & sim, 
        const math::Vec v){

    double d_v = math::abs(v), d_axis = math::abs(axis);

    double theta = acos((axis[0]*v[0]+axis[1]*v[1]+axis[2]*v[2])/(d_axis*d_v));

    // TODO: phi here is not relative to the mfv but to the z-axis
    // either remove or correct
    double dpx = math::dot(v,math::Vec(1, 0, 0));
    double dpy = math::dot(v,math::Vec(0, 1, 0));
    double phi = acos(dpx / sqrt((dpx*dpx)+(dpy*dpy)));
    if (dpy <0) {
      phi=-phi;
    }
   
    while (phi > math::Pi*2) {
      phi-=math::Pi*2;
    }
    while (phi < 0) {
      phi+=math::Pi*2;
    }

    while (theta > math::Pi*2) {
      theta-=math::Pi*2;
    }
    while (theta < 0) {
      theta+=math::Pi*2;
    }

    for (unsigned int i=0; i<sim.param().tfrdc.bins_phi.size()-1; i++) {
      double & bin_upper = sim.param().tfrdc.bins_phi[i+1];
       if (phi < bin_upper) {
         conf.special().tfrdc_mfv.dist_phi.at(i)+=1;
         break;
       }
    }

    for (unsigned int i=0; i<sim.param().tfrdc.bins_theta.size()-1; i++) {
      double & bin_upper = sim.param().tfrdc.bins_theta[i+1];
      if (theta < bin_upper) {
        conf.special().tfrdc_mfv.dist_theta.at(i)+=1;
        break;
      }
    }
}


/**
 * tensor-free RDC restraint interactions
 */
template<math::boundary_enum B, math::virial_enum V>
int _magnetic_field_vector_sd
(topology::Topology & topo,
        configuration::Configuration & conf,
        simulation::Simulation & sim,
        math::RandomGenerator* rng,
        int & err) {

    math::Periodicity<B> periodicity(conf.current().box);

    double exptaut = 0.0;                                       // [-]
    if (sim.param().tfrdc.tauth > 0.0)
      exptaut = exp(-sim.time_step_size() / sim.param().tfrdc.tauth);    // [-]
    const double & K = sim.param().tfrdc.Kmfv * 1e24;        // [kJ ps^2 / mol]
    const double dPavedP = 1.0 - exptaut;                       // [-]


    // initialize
    math::Vec & rh_i = conf.special().tfrdc_mfv.pos(0);         // [nm]
    math::Vec rh_ij;
    periodicity.nearest_image(rh_i, conf.special().tfrdc_mfv.pos(1), rh_ij);
    DEBUG(9, "initial rh_i  :" << math::v2s(rh_i));
    DEBUG(9, "initial rh_ij :" << math::v2s(rh_ij));
    double dh_ij = math::abs(rh_ij);                  // [nm]
    math::VArray r_ij;
    std::vector<double> costheta;
    std::vector<double> d_ij;

    // get the RDC bond vectors once, they won't change
    std::vector<topology::tf_rdc_restraint_struct>::iterator
    it = topo.tf_rdc_restraints().begin(),
          to = topo.tf_rdc_restraints().end();
    for(unsigned int l=0; it != to; ++it, ++l) {
      const math::Vec & r_i1 = it->v1.pos(conf, topo);         // [nm]
      math::Vec r_ij1;
      periodicity.nearest_image(r_i1, it->v2.pos(conf, topo), r_ij1);
      r_ij.push_back(r_ij1);
      d_ij.push_back(math::abs(r_ij1));                  // [nm]
      double ct = (rh_ij(0)*r_ij[l](0)+rh_ij(1)*r_ij[l](1)+rh_ij(2)*r_ij[l](2))/(dh_ij*d_ij[l]);
      costheta.push_back(ct);

      DEBUG(9, "r_i  :" << math::v2s(r_i1));
      DEBUG(9, "r_ij :" << math::v2s(r_ij1));
      // compute  P, P is dimensionless
      const double P = 0.5 * (3 * costheta[l] * costheta[l] - 1.0);          // [-]
      DEBUG(15, "P: " << P);
      DEBUG(15, "Theta of k1k2: " << acos(costheta[l]));

      // time-averaging
      // initialise average?
      double & P_expavg = conf.special().tfrdc_mfv.P_expavg[l];         // [-]
      if (sim.steps() == 0 && !sim.param().tfrdc.continuation) {
          // P_expavg = P;                                 // [-]
          P_expavg = it->D0 / interaction::D_c(it); // set initial value to the target value
      } else {
          DEBUG(15, "read-in P_expavg: " << P_expavg);
      }
      // apply time averaging
      P_expavg = dPavedP * P + exptaut * P_expavg;          // [-]
      DEBUG(15, " P_expavg: " << P_expavg);
    }

    // clear P_avg
    std::fill(conf.special().tfrdc_mfv.P_avg.begin(), conf.special().tfrdc_mfv.P_avg.end(), 0.0);
    std::fill(conf.special().tfrdc_mfv.dPdr_avg.begin(), conf.special().tfrdc_mfv.dPdr_avg.end(), math::Vec(0.0,0.0,0.0));
    for (int sdstep=0; sdstep < sim.param().tfrdc.nstsd; sdstep++) {

      // com translation removal
      if ((sdstep % conf.special().tfrdc_mfv.n_com_translation_removal) == 0) {
        math::Vec com_vel = (conf.special().tfrdc_mfv.vel[0]+conf.special().tfrdc_mfv.vel[1])/2;
        conf.special().tfrdc_mfv.vel[0]-= com_vel;
        conf.special().tfrdc_mfv.vel[1]-= com_vel;
      }
    
      // -----  calculate forces on magn. field vector -----
      math::VArray forces(2, math::Vec(0.0, 0.0, 0.0));

      // loop over the tensor-free RDC restraints
      std::vector<topology::tf_rdc_restraint_struct>::iterator
      it = topo.tf_rdc_restraints().begin(),
            to = topo.tf_rdc_restraints().end();
      for(unsigned int l = 0; it != to; ++it, ++l) {
          math::Vec dPdr(0.0, 0.0, 0.0);
          math::Vec force(0.0, 0.0, 0.0);
          double & P_expavg = conf.special().tfrdc_mfv.P_expavg[l];         // [-]

          // compute tensor-free RDC
          double RDC_avg = interaction::D_c(it) * P_expavg;     // [1 / ps]
          DEBUG(9, "RDC_avg (MFV): " << RDC_avg*pow(10,12));

          // the not averaged RDC
          //double RDC = interaction::D_c(it)  * P;

          DEBUG(15, "RDC atom: " << it->v1.atom(0));
          DEBUG(15, "Input D0 [Hz]: " << it->D0*pow(10,12));
          DEBUG(15, "Input dD0 [Hz]: " << it->dD0*pow(10,12));
          DEBUG(15, "D_k^c [Hz]: " << interaction::D_c(it)*pow(10,12));

          double term = 0.0;
          if (RDC_avg > it->D0 + it->dD0) {
            term = RDC_avg - it->D0 - it->dD0;            // [1 / ps]
          } else if (RDC_avg < it->D0 - it->dD0) {
            term = RDC_avg - it->D0 + it->dD0;            // [1 / ps]
          } else {
            DEBUG(9, "TFRDCRES  : restraint fulfilled");
          }

          for(unsigned int a = 0; a < 3; ++a){
            // eq 27, not to be used; TODO: remove
            //double inv_dh_ij = 1/dh_ij;
            //dPdr(a) = 3 * costheta * inv_dh_ij * (v_r_ij[l](a)/(d_ij) - costheta*rh_ij(a)*inv_dh_ij*inv_dh_ij); // [1 / nm]
            // eq. 32
            double prefix =  3 * costheta[l] / (d_ij[l]*dh_ij);
            dPdr(a) = prefix * r_ij[l](a); // [1 / nm]
          }
          

          DEBUG(9, "dP/dr" << math::v2s(dPdr));
          DEBUG(9, "dP_expavg/dP: " << dPavedP);

          force = term * interaction::D_c(it) * dPavedP * dPdr;  // [1 / (ps^2 nm)]

          forces(0) -= K * force;              // [(kJ ps^2 / mol) * (1 / (ps^2 nm))] = [kJ / (mol nm)]
          forces(1) -= forces(0);              // [kJ / (mol nm)]

          std::cout.precision(15); // useful for debugging
          DEBUG(7, "f_hi: " << math::v2s(forces(0)));
          DEBUG(7, "f_hj: " << math::v2s(forces(1)));
      }

      configuration::Configuration::special_struct::tfrdc_mfv_struct & tfrdc_mfv = conf.special().tfrdc_mfv;

      // Save old velocities and positions
      const math::VArray old_pos(tfrdc_mfv.pos);
      const math::VArray old_vel(tfrdc_mfv.vel);

      // -----  SD velocity 1 -----
      for (unsigned int i=0; i < 2; ++i){
        if(sim.param().tfrdc.cfrich != 0.0) {
          
          //we sample the V'i vector from eq. 2.11.2.20 from a Gaussian with 0.0 mean
          //and width sigma2_sq (sd1)
          rng->stddev(tfrdc_mfv.sd.sd1);
          tfrdc_mfv.sd.vrand1(i) = rng->get_gaussian_vec();
          //we sample the Vi vector to be used in eq. 2.11.2.2 from a Gaussian with 0.0 mean
          //and width rho1_sq (sd2)
          rng->stddev(tfrdc_mfv.sd.sd2);
          tfrdc_mfv.sd.vrand2(i) = rng->get_gaussian_vec();
          
          DEBUG(10, "vrand1=" <<  math::v2s(tfrdc_mfv.sd.vrand1(i)));
          DEBUG(10, "vrand2=" << math::v2s(tfrdc_mfv.sd.vrand2(i)));
                
          DEBUG(9, "old vel(" << i << ") " << math::v2s(old_vel(i)));
          DEBUG(9, "old pos(" << i << ") " << math::v2s(old_pos(i)));
          DEBUG(9, " force(" << i << ") " << math::v2s(forces(i)));
          DEBUG(10, "stochastic integral=" << math::v2s(tfrdc_mfv.stochastic_integral(i)));
          
          // svh .. V in eq. 65; vrand2 .. Yv
          math::Vec svh = tfrdc_mfv.stochastic_integral(i) * tfrdc_mfv.sd.c5 + tfrdc_mfv.sd.vrand2(i);
          //assign vrand to the Stochastic Integral array
          tfrdc_mfv.stochastic_integral(i) = tfrdc_mfv.sd.vrand1(i);
          
          //calculate new velocity using eq. 2.11.2.2 from the GROMOS96 book
          //CC1(NATTOT) = delivered with EXP(-GDT)  (GDT=GAM(J)*DT)
          //CC2(NATTOT) = delivered with (1-EXP(-GDT))/GDT
          //slightly rewritten form of first line of 2.11.2.3 minus
          //the first term in the third line of 2.11.2.3
          
          DEBUG(10, "svh=" << math::v2s(svh));
          DEBUG(10, "cf=" << tfrdc_mfv.sd.cf);
          
          tfrdc_mfv.vel(i) = (old_vel(i) - svh) * tfrdc_mfv.sd.c1;
          // force * m-1 * dt * (1-EXP(-GDT))/GDT, i.e.
          tfrdc_mfv.vel(i) += forces(i) * tfrdc_mfv.sd.cf;
          // last term of 2.11.2.2
          tfrdc_mfv.vel(i) += tfrdc_mfv.sd.vrand1(i);
        
          //we sample the R'i vector from eq. 2.11.2.25 from a Gaussian with 0.0 mean
          //and width rho2_sq (sd3)
          rng->stddev(tfrdc_mfv.sd.sd3);
          tfrdc_mfv.sd.vrand3(i) = rng->get_gaussian_vec();
          //we sample the Ri vector to be used in eq. 2.11.2.26 from a Gaussian with 0.0 mean
          //and width rho1_sq (sd4)
          rng->stddev(tfrdc_mfv.sd.sd4);
          tfrdc_mfv.sd.vrand4(i) = rng->get_gaussian_vec();
          
          DEBUG(10, "vrand3=" << math::v2s(tfrdc_mfv.sd.vrand3(i)));
          DEBUG(10, "vrand4=" << math::v2s(tfrdc_mfv.sd.vrand4(i)));
        } else {
          tfrdc_mfv.vel(i) = old_vel(i) + forces(i) * sim.time_step_size() / topo.mass(i);
          tfrdc_mfv.sd.vrand1(i)=math::Vec(0.0, 0.0, 0.0);
          tfrdc_mfv.sd.vrand2(i)=math::Vec(0.0, 0.0, 0.0);
          tfrdc_mfv.sd.vrand3(i)=math::Vec(0.0, 0.0, 0.0);
          tfrdc_mfv.sd.vrand4(i)=math::Vec(0.0, 0.0, 0.0);
        }
      } // loop over atoms

      // sd_pos1
      for (unsigned int i=0; i < 2; ++i){
        //calc new positions
        //according to step 7 in leap frog for SD (GROMOS96 book)
        tfrdc_mfv.pos(i) = old_pos(i) 
              + tfrdc_mfv.vel(i) * sim.time_step_size() * tfrdc_mfv.sd.c6;    
        DEBUG(9, "old pos(" << i << ") " << math::v2s(old_pos(i)));
        DEBUG(9, "new pos(" << i << ") " << math::v2s(tfrdc_mfv.pos(i)));      
      } // loop over atoms

      // shake (only positions)

      err = _shake_bond<B>(conf.current().box, tfrdc_mfv.pos, old_pos, tfrdc_mfv.mass, tfrdc_mfv.d);
      if (err) break;

      // sd_vel2
      for (unsigned int i=0; i < 2; ++i) {
        double cinv = 1.0 / (tfrdc_mfv.sd.c6 * sim.time_step_size());
        math::Vec r;
        periodicity.nearest_image(tfrdc_mfv.pos(i), old_pos(i), r);
        tfrdc_mfv.vel(i) = r * cinv;
        DEBUG(10, "velocity SHAKEN" << math::v2s(tfrdc_mfv.vel(i)))
      }

      // sd_pos2
      if(sim.param().tfrdc.cfrich != 0.0) {
        for (unsigned int i=0; i < 2; ++i){
          //this is 2.11.2.25
          math::Vec sxh = tfrdc_mfv.stochastic_integral(i) * tfrdc_mfv.sd.c9;
          sxh += tfrdc_mfv.sd.vrand4(i);
          tfrdc_mfv.stochastic_integral(i) = tfrdc_mfv.sd.vrand3(i);  
          //this is 2.11.2.26
          tfrdc_mfv.pos(i) += tfrdc_mfv.sd.vrand3(i) - sxh;
        } // loop over atoms
      }

      // shake (only positions)
      err = _shake_bond<B>(conf.current().box, tfrdc_mfv.pos, old_pos, tfrdc_mfv.mass, tfrdc_mfv.d);
      if (err) break;

      
      periodicity.nearest_image(rh_i, conf.special().tfrdc_mfv.pos(1), rh_ij);
      dh_ij = math::abs(rh_ij);                  // [nm]

    math::Vec axis;
    math::Vec axis1 = topo.tf_rdc_molaxis()[0].pos(conf, topo), 
                      axis2 = topo.tf_rdc_molaxis()[1].pos(conf, topo);
    periodicity.nearest_image(axis1, axis2, axis);
      _add_to_theta_distribution(axis, conf, sim, rh_ij);
      DEBUG(9, "rh_i  :" << math::v2s(rh_i));
      DEBUG(9, "rh_ij :" << math::v2s(rh_ij));
      
      it = topo.tf_rdc_restraints().begin(),
            to = topo.tf_rdc_restraints().end();
      for(unsigned int l=0; it != to; ++it, ++l) {

        costheta[l] = (rh_ij(0)*r_ij[l](0)+rh_ij(1)*r_ij[l](1)+rh_ij(2)*r_ij[l](2))/(dh_ij*d_ij[l]);
 
        const double P = 0.5 * (3 * costheta[l] * costheta[l] - 1.0);          // [-]
        DEBUG(15, "P: " << P);
        DEBUG(15, "Theta of k1k2: " << acos(costheta[l]));

        // apply time averaging
        double & P_expavg = conf.special().tfrdc_mfv.P_expavg[l];         // [-]
        P_expavg = dPavedP * P + exptaut * P_expavg;          // [-]
        DEBUG(15, " P_expavg: " << P_expavg);

        conf.special().tfrdc_mfv.P_avg[l]+=P;
        DEBUG(15, " P_avg: " << conf.special().tfrdc_mfv.P_avg[l]/(sdstep+1));

        double prefix =  3 * costheta[l] / (d_ij[l]*dh_ij);
        // TODO: add case when bond is not constrained
        for(unsigned int a = 0; a < 3; ++a){
            conf.special().tfrdc_mfv.dPdr_avg[l](a)+=prefix*rh_ij(a);
        }
      }

      if (err) return err;
    } // loop sd steps

    // average
    it = topo.tf_rdc_restraints().begin(),
          to = topo.tf_rdc_restraints().end();
    for(unsigned int l=0; it != to; ++it, ++l) {
      conf.special().tfrdc_mfv.P_avg[l]/=sim.param().tfrdc.nstsd;
      for(unsigned int a = 0; a < 3; ++a){
          conf.special().tfrdc_mfv.dPdr_avg[l](a)/=sim.param().tfrdc.nstsd;
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
  int error = 0;
  m_timer.start("mfv");
  if (sim.param().tfrdc.nstsd > 0) {
    SPLIT_VIRIAL_BOUNDARY(_magnetic_field_vector_sd,
          topo, conf, sim, m_rng, error);
    if (error) return error;

    conf.special().tfrdc_mfv.sd.seed=m_rng->seed();
  }
  m_timer.stop("mfv");

  SPLIT_VIRIAL_BOUNDARY(_calculate_tf_rdc_restraint_interactions,
          topo, conf, sim);
  m_timer.stop();
  return error;
}

void _init_mfv_sd(configuration::Configuration & conf, 
                                simulation::Simulation & sim,
                                math::RandomGenerator* rng) {

  configuration::Configuration::special_struct::tfrdc_mfv_struct & tfrdc_mfv = conf.special().tfrdc_mfv;

  // generate SD coefficients
  const double gg = fabs(sim.param().tfrdc.cfrich);
  const double gdt = gg * sim.time_step_size();
  const double gdth = gdt * 0.5;
  
  // TODO: maybe do not use hardcoded values here
  tfrdc_mfv.mass = sim.param().tfrdc.mfv_mass;
  tfrdc_mfv.d = sim.param().tfrdc.mfv_r;
  tfrdc_mfv.n_com_translation_removal = 1000;
  tfrdc_mfv.sd.kToverM = sqrt(math::k_Boltzmann * sim.param().tfrdc.tempsd / tfrdc_mfv.mass);
  DEBUG(12, "kb "<< math::k_Boltzmann << "Temp " << sim.param().tfrdc.tempsd << "mass "<< tfrdc_mfv.mass);

  const double emdth = exp(-gdth);
  const double emdt  = emdth * emdth;
  const double cdth  = gdt - 3.0 + 4.0 * emdth - emdt;

  if (fabs(gdt) > 0.05){
    DEBUG(12, "\tdoing the analytical formulas");
    const double epdth = exp(+gdth);
    const double epdt  = epdth * epdth;
    const double omdt  = 1.0 - emdt;
    const double ddth  = 2.0 - epdth - emdth;
    const double bpdth = gdt * (epdt - 1.0) - 4.0 * (epdth - 1.0) * (epdth - 1.0);
    const double bmdth = gdt * (1.0 - emdt) - 4.0 * (emdth -1.0) * (emdth -1.0);
    tfrdc_mfv.sd.c1 = emdt;
    tfrdc_mfv.sd.c2 = omdt / gdt;
    tfrdc_mfv.sd.c3 = sqrt(fabs(omdt));
    tfrdc_mfv.sd.c4 = sqrt(fabs(bpdth/cdth));
    tfrdc_mfv.sd.c5 = gg * ddth/cdth;
    tfrdc_mfv.sd.c6 = (epdth - emdth) / gdt;
    tfrdc_mfv.sd.c7 = sqrt(fabs(cdth)) / gg;
    tfrdc_mfv.sd.c8 = sqrt(fabs(bmdth/omdt)) / gg;
    tfrdc_mfv.sd.c9 = -ddth/(gg * omdt);
    
  } else {
    DEBUG(12, "\tdoing the power series");
    //this is a power series expansion for the coefficients used
    //in the SD algortihm. it should be used when gamdt < 0.05, otherwise
    //the analytical expression from above is good enough
    const double gdth2 = gdth * gdth;
    const double gdth3 = gdth2 * gdth;
    const double gdth4 = gdth2 * gdth2;
    const double gdth5 = gdth2 * gdth3;
    const double gdth6 = gdth3 * gdth3;
    const double gdth7 = gdth4 * gdth3;

    tfrdc_mfv.sd.c1 = exp(-gdt);

    tfrdc_mfv.sd.c2 = 1.0 - gdth + gdth2 * 2.0/3.0 - gdth3/3.0 
+ gdth4 * 2.0/15.0 - gdth5 * 2.0/45.0 + gdth6 * 4.0/315.0;

    tfrdc_mfv.sd.c3 = sqrt(fabs(tfrdc_mfv.sd.c2) * 2.0 * gdth);

    tfrdc_mfv.sd.c4 = sqrt(fabs(gdth/2.0 + gdth2 * 7.0/8.0 + gdth3 * 367.0/480.0 + gdth4 * 857.0/1920.0 
          + gdth5 * 52813.0/268800.0 + gdth6 * 224881.0/3225600.0 +
          gdth7 * 1341523.0/64512000.0));

    tfrdc_mfv.sd.c5 = -2.0 / sim.time_step_size() *
(1.5 + gdth * 9.0/8.0 + gdth2 * 71.0/160.0 + gdth3 * 81.0/640.0 + gdth4 * 7807.0/268800.0 
  + gdth5 * 1971.0/358400.0 + gdth6 * 56417.0/64512000.0);

    tfrdc_mfv.sd.c6 = 1.0 + gdth2/6.0 + gdth4/10.0 + gdth6/5040.0;

    tfrdc_mfv.sd.c7 = sim.time_step_size() * 0.5 * 
sqrt(fabs(gdth * 2.0/3.0 - gdth2/2.0 + gdth3 * 7.0/30.0 -gdth4/12.0 + gdth5 * 31.0/1260.0 - 
    gdth6/160.0 + gdth7 * 127.0/90720.0));

    tfrdc_mfv.sd.c8 = sim.time_step_size() * 0.5 * 
sqrt(fabs(gdth/6.0 - gdth3/60.0 + gdth5 * 17.0/10080.0 - gdth7 * 31.0/181440.0));

    tfrdc_mfv.sd.c9 = sim.time_step_size() * 0.5 *
(0.5 + gdth/2.0 + gdth2 * 5.0/24.0 + gdth3/24.0 + gdth4/240.0 + gdth5/720.0 + gdth6 * 5.0/8064.0);

  } 



  //CC2(NATTOT) = delivered with (1-EXP(-GDT))/GDT
  tfrdc_mfv.sd.cf = sim.time_step_size() * tfrdc_mfv.sd.c2 / tfrdc_mfv.mass;
  
  //CC3(NATTOT) = delivered with SQRT(1-EXP(-GDT))
  //this is 2.11.2.8
  tfrdc_mfv.sd.sd1 = tfrdc_mfv.sd.c3 * tfrdc_mfv.sd.kToverM;
  
  //CC4(NATTOT) = delivered with SQRT(B(+GDT/2)/C(+GDT/2)) (SEE PAPER)
  // this is 2.11.2.9
  tfrdc_mfv.sd.sd2 = tfrdc_mfv.sd.kToverM * tfrdc_mfv.sd.c4;

  tfrdc_mfv.sd.sd3 = tfrdc_mfv.sd.kToverM * tfrdc_mfv.sd.c7;
  tfrdc_mfv.sd.sd4 = tfrdc_mfv.sd.kToverM * tfrdc_mfv.sd.c8;

  tfrdc_mfv.sd.vrand1.resize(2);
  tfrdc_mfv.sd.vrand2.resize(2);
  tfrdc_mfv.sd.vrand3.resize(2);
  tfrdc_mfv.sd.vrand4.resize(2);
  
  DEBUG(10, "\tkT/m: " << tfrdc_mfv.sd.kToverM <<  " sd1: " << tfrdc_mfv.sd.sd1 << " sd2: " << tfrdc_mfv.sd.sd2);
  DEBUG(10, "sd3: " << tfrdc_mfv.sd.sd3 << " sd4: " << tfrdc_mfv.sd.sd4);
  DEBUG(10, "c5=" << tfrdc_mfv.sd.c5);

  if (!sim.param().tfrdc.continuation){
    tfrdc_mfv.pos.resize(2);
    tfrdc_mfv.vel.resize(2);
    tfrdc_mfv.stochastic_integral.resize(2);
    // set initial positions
    tfrdc_mfv.pos(0)=math::Vec(0,0,0);
    tfrdc_mfv.pos(1)=math::Vec(0,0,tfrdc_mfv.d);

    // initialize velocities, stochastic integrals
    for (unsigned int i = 0; i < 2; ++i){

      const double sd = sqrt(math::k_Boltzmann * sim.param().tfrdc.tempsd / tfrdc_mfv.mass);
      rng->stddev(sd);
      for(int d=0; d<3; ++d){
        tfrdc_mfv.vel(i)(d) = rng->get_gauss();
      }      

      DEBUG(12, "stochastic coefficient for atom i = " << i);
      // just an initial guess. no need(?) for high precision??
      DEBUG(10, "gdt = " << gdt << " in initial stochastic integral");
      if (sim.param().tfrdc.cfrich < math::epsilon){
        tfrdc_mfv.stochastic_integral(i) = 0.0;
      }
      else {
        const double sd = tfrdc_mfv.sd.kToverM / sim.param().tfrdc.cfrich * sqrt(cdth); 
        DEBUG(10, "sigma " << sd << " ktoverm "<< tfrdc_mfv.sd.kToverM);
        rng->stddev(sd);
        tfrdc_mfv.stochastic_integral(i) = rng->get_gaussian_vec();
        DEBUG(10, "initial stochastic integral: " << v2s(tfrdc_mfv.stochastic_integral(i)));
      }
    }
  }
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
  conf.special().tfrdc_mfv.P_expavg.resize(num_res);
  conf.special().tfrdc_mfv.P_avg.resize(num_res);
  conf.special().tfrdc_mfv.dPdr_avg.resize(num_res);


  if (!sim.param().tfrdc.read) {
    conf.special().tfrdc.num_averaged=0;
  }
  if (sim.param().tfrdc.nstsd > 0) {
    // Create random number generator
    m_rng = math::RandomGenerator::create(sim.param(), "0");
    m_rng->mean(0.0);
    if (!sim.param().tfrdc.continuation) {
      // use seed from input file
      std::ostringstream seed;
      seed << sim.param().start.ig;
        m_rng->seed(seed.str());
    } else {
        m_rng->seed(conf.special().tfrdc_mfv.sd.seed);
    }

    _init_mfv_sd(conf, sim, m_rng);
    sim.param().tfrdc.bins_theta = _linspace(0.0,math::Pi,101.0);
    sim.param().tfrdc.bins_phi = _linspace(0.0,2*math::Pi,101.0);

    conf.special().tfrdc_mfv.dist_theta.resize(sim.param().tfrdc.bins_theta.size(), 0.0);
    conf.special().tfrdc_mfv.dist_phi.resize(sim.param().tfrdc.bins_phi.size(), 0.0);
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
    if (sim.param().tfrdc.nstsd > 0) {
      os << "\tsampling magnetic-field vector orientations in " << sim.param().tfrdc.nstsd << " SD steps" << std::endl;
    }

    os.precision(8);
    os << "  - Number of restraints: " << num_res << std::endl
            << "  - force constant: " << std::setw(15) << sim.param().tfrdc.K << std::endl;
    if (sim.param().tfrdc.nstsd > 0) 
      os << "  - force constant mfv " << std::setw(15) << sim.param().tfrdc.Kmfv << std::endl;
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
