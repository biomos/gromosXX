/**
 * @file rdc_restraint_interaction.cc
 */


#include <time.h>

#include <stdheader.h>
#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <simulation/simulation.h>
#include <configuration/configuration.h>
#include <interaction/interaction.h>
#include <interaction/interaction_types.h>
#include <interaction/forcefield/forcefield.h>
#include <interaction/special/rdc_restraint_interaction.h>
#include <math/periodicity.h>
#include <math/random.h>
#include <util/error.h>
#include <util/template_split.h>
#include <create_special.h>
#include <topology/core/body_term.h>

#include <util/debug.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_errno.h>

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE special

//#ifndef FINITE_DIFF
//#define FINITE_DIFF

// FIXME remove
//#ifndef SCRIPT
//#define SCRIPT

using namespace std;

namespace math{
  //E.R. Cohen, T. Cvitas, J.G. Frey, B. Holmström, K. Kuchitsu,
  //R. Marquardt, I. Mills, F. Pavese, M. Quack, J. Stohner,
  //H.L. Strauss, M. Takami, and A.J. Thor, "Quantities, Units and
  //Symbols in Physical Chemistry", IUPAC Green Book, 3rd Edition,
  //2nd Printing, IUPAC & RSC Publishing, Cambridge (2008)
  const double n_avogadro = 6.02214179e23; // [1/mol]
}

namespace interaction {

inline double flat_bottom_pot(const double x, const double x0, const double delta){
#if (__cplusplus > 199711L) // we have c++11 or newer
  return 0.5 * pow(max( {x-x0-delta, -(x-x0+delta), 0.0} ), 2);
#else
  return 0.5 * pow(max(min(x-x0+delta, 0.0), x-x0-delta), 2);
#endif
}
inline double dflat_bottom_pot(const double x, const double x0, const double delta){
#if (__cplusplus > 199711L) // we have c++11 or newer
  return max( {x-x0-delta, -(x-x0+delta), 0.0} );
#else
  return max(min(x-x0+delta, 0.0), x-x0-delta);
#endif
}

/**
 * Helper functions
 */

// This function is for debugging only
void write_mf_forces(const math::VArray &force_vectors_hfield, const char *file_name){
  ofstream mf_vector_file;
  mf_vector_file.open(file_name);//, ofstream::app);

  if(mf_vector_file.is_open()){
    mf_vector_file << "FORCES ON MF VECTORS" << endl;
    for (unsigned int h=0; h<force_vectors_hfield.size(); ++h) {
      for (int i=0; i<3; ++i) {
        mf_vector_file << setw(15) << setprecision(5) << scientific << force_vectors_hfield[h][i];

      }
      mf_vector_file << endl;
    }
    mf_vector_file << "END" << endl;
  }
  else{
    cout << "COULDN'T OPEN FILE" << endl;
  }

  mf_vector_file.close();
}

// This function is for debugging only
void write_tensor_forces(const vector<double> &force_array_tensor, const char *file_name){
  ofstream force_tensor_file;
  force_tensor_file.open(file_name);//, ofstream::app);

  if(force_tensor_file.is_open()){
    force_tensor_file << "FORCES ON TENSOR ELEMENTS" << endl;
    for (int h=0; h<5; ++h) {
      force_tensor_file << setw(15) << setprecision(5) << scientific << force_array_tensor[h];
    }
    force_tensor_file << endl;
    force_tensor_file << "END" << endl;
  }
  else{
    cout << "COULDN'T OPEN FILE" << endl;
  }

  force_tensor_file.close();
}

// This function is for debugging only
void write_clm_forces(const vector<double> &force_array_sh, const char *file_name){
  ofstream force_sh_file;
  force_sh_file.open(file_name);//, ofstream::app);

  if(force_sh_file.is_open()){
    force_sh_file << "FORCES ON SH ELEMENTS" << endl;
    for (int lm=0; lm<5; ++lm) {
      force_sh_file << setw(15) << setprecision(5) << scientific << force_array_sh[lm];
    }
    force_sh_file << endl;
    force_sh_file << "END" << endl;
  }
  else{
    cerr << "COULDN'T OPEN FILE" << endl;
  }

  force_sh_file.close();
}


template<math::boundary_enum B>
struct fit_param {
  vector<vector<topology::rdc_restraint_struct> >::iterator topo_it;
  vector<configuration::Configuration::special_struct::rdc_struct>::iterator conf_it;
  configuration::Configuration * conf;
  simulation::Simulation * sim;
};


//*******************************************************************************************************
//*******************************************************************************************************
//*************************                                                 *****************************
//*************************    M A G N E T I C  F I E L D  V E C T O R S    *****************************
//*************************                                                 *****************************
//*******************************************************************************************************
//*******************************************************************************************************


// calculate forces on magnetic field vectors in the MF representation
template<math::boundary_enum B>
void _calculate_forces_vectors_MF(topology::Topology & topo,
                    configuration::Configuration & conf,
                    const simulation::Simulation & sim,
                  math::VArray & force_vectors_hfield /*empty*/) {

  vector<configuration::Configuration::special_struct::rdc_struct>::iterator
      conf_it = conf.special().rdc.begin(),
      conf_to = conf.special().rdc.end();
  vector<vector<topology::rdc_restraint_struct> >::iterator
      topo_it = topo.rdc_restraints().begin();
  for( ; conf_it!=conf_to; conf_it++, topo_it++){

    // number of magnetic field vectors
    const int number_hvectors = conf_it->MFpoint.size();

    force_vectors_hfield.clear();
    // initialise the vector holding hfield force vectors
    for (int h=0; h<number_hvectors; ++h) {
      force_vectors_hfield.push_back(math::Vec(0.0, 0.0, 0.0));
    }

    // create variables
    math::Periodicity<B> periodicity(conf.current().box);
    math::VArray &pos = conf.current().pos;
    math::Vec ri(0.0, 0.0, 0.0);
    math::Vec rj(0.0, 0.0, 0.0);
    math::Vec rij(0.0, 0.0, 0.0);
    math::Vec rh(0.0, 0.0, 0.0);

    ///////////////
    // MAIN LOOP //
    ///////////////
    int k=0;
    vector<topology::rdc_restraint_struct>::iterator
      it = topo_it->begin(),
      to = topo_it->end();
    for(; it!=to; ++it, ++k) {
      // Get positions of the two atoms involved in this RDC
      ri = pos(it->i);
      rj = pos(it->j);
      periodicity.nearest_image(ri, rj, rij);
      const double dij = math::abs(rij);
      const double dij2 = dij*dij;
      const double dij3 = dij2*dij;
      const double dij5 = dij2*dij3;

      conf_it->curr[k] = 0.0;
      // Loop over all magnetic field vectors
      for(int h=0; h<number_hvectors; ++h) {
        rh = conf_it->MFpoint.at(h); // [nm]
        conf_it->curr[k] +=  pow(math::dot(rij, rh) / dij, 2) - 1.0/3.0; // [1]
      } // Loop over MF
      conf_it->curr[k] *= RDC_max(it)/dij3 * 1.0/number_hvectors * 1.5; // [1/ps]

      // locally calculate average RDC values, use for forces and don't export
      // the average RDC is not updated because that is done when forces on atoms are computed
      double local_average = 0.0;
      if(sim.param().rdc.mode == simulation::rdc_restr_av ||
         sim.param().rdc.mode == simulation::rdc_restr_av_weighted ||
         sim.param().rdc.mode == simulation::rdc_restr_biq ||
         sim.param().rdc.mode == simulation::rdc_restr_biq_weighted){
        const double exp = pow(M_E, -sim.time_step_size()/sim.param().rdc.tau); // [1]
        local_average = exp * conf_it->av[k] + (1.0 - exp) * conf_it->curr[k]; // [1/ps]
        DEBUG(15, "Delta t, tau, e^(- Delta t/tau): " << scientific << sim.time_step_size() << ", " <<  sim.param().rdc.tau << ", " << pow(M_E, -sim.time_step_size()/sim.param().rdc.tau) )
        DEBUG(10, "interaction " << k << ":  R-inst, R-av: " << conf_it->curr[k] << ", " << local_average )
      }

      // K[kJ*ps^2], N_A[1/mol]
      const double force_coefficient = sim.param().rdc.K * it->weight * math::n_avogadro; // [kJ*ps^2/mol]

      // 'forces' on the MF vectors
      // f = - K * (RDC - RDC_0) [ * inner derivative]
      for(int h=0; h<number_hvectors; ++h) {
        rh = conf_it->MFpoint.at(h); // [nm]
        const math::Vec dD_dh = RDC_max(it) * 1.0/number_hvectors * 3.0/dij5 * math::dot(rij, rh) * rij; // [kJ/(nm*mol)] (some |rh| are omitted); //  [kJ/(nm*mol)]
        // forces on MF vectors depending on AV mode
        if(sim.param().rdc.mode == simulation::rdc_restr_inst ||
           sim.param().rdc.mode == simulation::rdc_restr_inst_weighted){
          const double dV_dD = force_coefficient * dflat_bottom_pot(conf_it->curr[k], it->R0, sim.param().rdc.delta); // [kJ/(nm^2*mol)]
          force_vectors_hfield[h] -= dV_dD * dD_dh; // [kJ/(nm*mol)]
          DEBUG(15, "dV_dD; dD_dh: " << dV_dD << ", " << scientific << "(" << dD_dh[0] << ", " << dD_dh[1] << ", " << dD_dh[2] << ")")
        }
        else if(sim.param().rdc.mode == simulation::rdc_restr_av ||
                sim.param().rdc.mode == simulation::rdc_restr_av_weighted){
          const double dVav_dDav = force_coefficient * dflat_bottom_pot(local_average, it->R0, sim.param().rdc.delta); // [kJ/(nm^2*mol)]
          // in the following we set dDav_dD to either 1.0 or [1-e^(-dt/tau)].
          // The latter is correct, the former allows for using the same K as without time AV
          const double dDav_dD = sim.param().rdc.tAVfactor ? (1.0-pow(M_E, -sim.time_step_size()/sim.param().rdc.tau)) : 1.0;
          force_vectors_hfield[h] -= dVav_dDav * dDav_dD * dD_dh; // [kJ/(nm*mol)]
          DEBUG(15, "dVav_dDav; dDav_dD; dD_dh: " << dVav_dDav << ", " << dDav_dD << ", " << scientific << "(" << dD_dh[0] << ", " << dD_dh[1] << ", " << dD_dh[2] << ")")
        }
        else if(sim.param().rdc.mode == simulation::rdc_restr_biq ||
                sim.param().rdc.mode == simulation::rdc_restr_biq_weighted){
          const double dVbq_dD = 2.0 * force_coefficient * dflat_bottom_pot(conf_it->curr[k], it->R0, sim.param().rdc.delta) * flat_bottom_pot(local_average, it->R0, sim.param().rdc.delta); // [kJ/(nm^2*mol)]
          const double dVbq_dDav = 2.0 * force_coefficient * flat_bottom_pot(conf_it->curr[k], it->R0, sim.param().rdc.delta) * dflat_bottom_pot(local_average, it->R0, sim.param().rdc.delta); // [kJ/(nm^2*mol)]
          // in the following we set dDav_dD to either 1.0 or [1-e^(-dt/tau)] or 0.
          // The first is simple, the second correct.
          double dDav_dD = 0;  // case 2
          if (sim.param().rdc.tAVfactor == 0) dDav_dD = 1.0;
          else if (sim.param().rdc.tAVfactor == 1) dDav_dD = 1.0-pow(M_E, -sim.time_step_size()/sim.param().rdc.tau);
          force_vectors_hfield[h] -= (dVbq_dD + dVbq_dDav * dDav_dD) * dD_dh; // [kJ/(nm*mol)]
          DEBUG(15, "dVbq_dD; dVbq_dDav; dDav_dD; dD_dh: " << dVbq_dD << ", " << dVbq_dDav << ", "<< dDav_dD << ", " << scientific << "(" << dD_dh[0] << ", " << dD_dh[1] << ", " << dD_dh[2] << ")")
        }
        DEBUG(10, "mf vector (" << h << "): " << math::v2s(conf_it->MFpoint.at(h)) )
        DEBUG(10, "forces on the mf vectors(" << h << "): " << math::v2s(force_vectors_hfield[h]) )
        DEBUG(10, "velocity of vector component(" << h << "): " << math::v2s(conf_it->MFpointVel[h]))
      } // Loop over MF (Vphys)
    } // Loop over RDC

    for (int i=0; i<number_hvectors; ++i) {
      DEBUG(10, "position of h-vectors(" << i << "): " << math::v2s(conf_it->MFpoint.at(i)))
      DEBUG(10, "forces on h-vectors(" << i << "): " << math::v2s(force_vectors_hfield[i]))
    }
  }// loop over rdc groups

  //write_mf_forces(force_vectors_hfield, "trfmf");
}


// calculate forces on atoms in the MF representation
template<math::boundary_enum B>
void _calculate_forces_atoms_MF(topology::Topology & topo,
                    configuration::Configuration & conf,
                    const simulation::Simulation & sim) {

  vector<configuration::Configuration::special_struct::rdc_struct>::iterator
      conf_it = conf.special().rdc.begin(),
      conf_to = conf.special().rdc.end();
  vector<vector<topology::rdc_restraint_struct> >::iterator
      topo_it = topo.rdc_restraints().begin();
  for( ; conf_it!=conf_to; conf_it++, topo_it++){

    // number of magnetic field vectors
    const int number_hvectors = conf_it->MFpoint.size();

    // create variables
    math::Periodicity<B> periodicity(conf.current().box);
    math::VArray &pos = conf.current().pos;
    math::Vec ri(0.0, 0.0, 0.0);
    math::Vec rj(0.0, 0.0, 0.0);
    math::Vec rij(0.0, 0.0, 0.0);
    math::Vec rh(0.0, 0.0, 0.0);

    ///////////////
    // MAIN LOOP //
    ///////////////
    int k=0;
    vector<topology::rdc_restraint_struct>::iterator
      it = topo_it->begin(),
      to = topo_it->end();
    for(; it!=to; ++it, ++k) {
      // Get positions of the two atoms involved in this RDC
      ri = pos(it->i);
      rj = pos(it->j);
      periodicity.nearest_image(ri, rj, rij);
      const double dij = math::abs(rij);
      const double dij2 = dij*dij;
      const double dij3 = dij2*dij;
      const double dij5 = dij2*dij3;

      conf_it->curr[k] = 0.0;
      // Loop over all magnetic field vectors
      for(int h=0; h<number_hvectors; ++h) {
        rh = conf_it->MFpoint.at(h); // [nm]
        conf_it->curr[k] +=  pow(math::dot(rij, rh) / dij, 2) - 1.0/3.0; // [1]
      } // Loop over MF
      conf_it->curr[k] *= RDC_max(it)/dij3 * 1.0/number_hvectors * 1.5; // [1/ps]

      // locally calculate average RDC values, use for forces and don't export
      // the average RDC is not updated because that is done when forces on atoms are computed
      if(sim.param().rdc.mode == simulation::rdc_restr_av ||
         sim.param().rdc.mode == simulation::rdc_restr_av_weighted ||
         sim.param().rdc.mode == simulation::rdc_restr_biq ||
         sim.param().rdc.mode == simulation::rdc_restr_biq_weighted){
        const double exp = pow(M_E, -sim.time_step_size()/sim.param().rdc.tau); // [1]
        conf_it->av[k] *= exp; // [1/ps]
        conf_it->av[k] += (1.0 - exp) * conf_it->curr[k]; // [1/ps]
        DEBUG(15, "Delta t, tau, e^(- Delta t/tau): " << scientific << sim.time_step_size() << ", " <<  sim.param().rdc.tau << ", " << pow(M_E, -sim.time_step_size()/sim.param().rdc.tau) )
        DEBUG(10, "interaction " << k << ":  R-inst, R-av: " << conf_it->curr[k] << ", " << conf_it->av[k] )
      }

      // K[kJ*ps^2], N_A[1/mol]
      const double force_coefficient = sim.param().rdc.K * it->weight * math::n_avogadro; // [kJ*ps^2/mol]

      if(sim.param().rdc.mode == simulation::rdc_restr_inst ||
         sim.param().rdc.mode == simulation::rdc_restr_inst_weighted){
        // rdc_energy = .5 * K * (RDC - RDC_0)^2
        // split the interaction energy between the two atoms and possibly between two energy groups
        conf.current().energies.rdc_energy[topo.atom_energy_group()[it->i]] += 0.5 * force_coefficient * flat_bottom_pot(conf_it->curr[k], it->R0, sim.param().rdc.delta); // [kJ/mol]
        conf.current().energies.rdc_energy[topo.atom_energy_group()[it->j]] += 0.5 * force_coefficient * flat_bottom_pot(conf_it->curr[k], it->R0, sim.param().rdc.delta); // [kJ/mol]
        DEBUG(15, "energy added to groups " << topo.atom_energy_group()[it->i] << " and " << topo.atom_energy_group()[it->j] )
        DEBUG(15, "group energies are " << conf.current().energies.rdc_energy[topo.atom_energy_group()[it->i]] << " and " << conf.current().energies.rdc_energy[topo.atom_energy_group()[it->j]])
        DEBUG(10, "interaction_energy (inst): " << scientific << force_coefficient * flat_bottom_pot(conf_it->curr[k], it->R0, sim.param().rdc.delta))
      }
      else if(sim.param().rdc.mode == simulation::rdc_restr_av ||
              sim.param().rdc.mode == simulation::rdc_restr_av_weighted){
        // rdc_energy = .5 * K * (<RDC> - RDC_0)^2
        conf.current().energies.rdc_energy[topo.atom_energy_group()[it->i]] += 0.5 * force_coefficient * flat_bottom_pot(conf_it->av[k], it->R0, sim.param().rdc.delta); // [kJ/mol]
        conf.current().energies.rdc_energy[topo.atom_energy_group()[it->j]] += 0.5 * force_coefficient * flat_bottom_pot(conf_it->av[k], it->R0, sim.param().rdc.delta); // [kJ/mol]
        DEBUG(15, "energy added to groups " << topo.atom_energy_group()[it->i] << " and " << topo.atom_energy_group()[it->j] )
        DEBUG(15, "group energies are " << conf.current().energies.rdc_energy[topo.atom_energy_group()[it->i]] << " and " << conf.current().energies.rdc_energy[topo.atom_energy_group()[it->j]])
        DEBUG(10, "interaction_energy (av): " << scientific << force_coefficient * flat_bottom_pot(conf_it->av[k], it->R0, sim.param().rdc.delta))
      }
      else if(sim.param().rdc.mode == simulation::rdc_restr_biq ||
              sim.param().rdc.mode == simulation::rdc_restr_biq_weighted){
        // rdc_energy = .5 * K * (RDC - RDC_0)^2 * (<RDC> - RDC_0)^2
        conf.current().energies.rdc_energy[topo.atom_energy_group()[it->i]] += 0.5 * 2.0 * force_coefficient * flat_bottom_pot(conf_it->curr[k], it->R0, sim.param().rdc.delta) * flat_bottom_pot(conf_it->av[k], it->R0, sim.param().rdc.delta); // [kJ/mol]
        conf.current().energies.rdc_energy[topo.atom_energy_group()[it->j]] += 0.5 * 2.0 * force_coefficient * flat_bottom_pot(conf_it->curr[k], it->R0, sim.param().rdc.delta) * flat_bottom_pot(conf_it->av[k], it->R0, sim.param().rdc.delta); // [kJ/mol]
        DEBUG(15, "energy added to groups " << topo.atom_energy_group()[it->i] << " and " << topo.atom_energy_group()[it->j] )
        DEBUG(15, "group energies are " << conf.current().energies.rdc_energy[topo.atom_energy_group()[it->i]] << " and " << conf.current().energies.rdc_energy[topo.atom_energy_group()[it->j]])
        DEBUG(10, "interaction_energy (biq): " << scientific << 2.0 * force_coefficient * flat_bottom_pot(conf_it->curr[k], it->R0, sim.param().rdc.delta) * flat_bottom_pot(conf_it->av[k], it->R0, sim.param().rdc.delta))
      }

      // forces on atoms
      math::Vec force_vector_atom(0.0, 0.0, 0.0);
      for (int h=0; h<number_hvectors; ++h) {
        rh = conf_it->MFpoint.at(h);
        const double dot_product = math::dot(rij, rh); // [nm^2]
        // the second half of the following expression is zero if rij is of constant length
#ifdef FINITE_DIFF
        force_vector_atom += dot_product * rh  - (2.5 * pow(dot_product/dij, 2) - 0.5) * rij; // [1]
#else
        force_vector_atom += dot_product * rh; // [1]
#endif // FINITE_DIFF

      } // Loop over MF

      // Add the force contribution to the forces on the atoms i and j
      // and here the factor of r^-2 comes in again
      const math::Vec dD_dr = RDC_max(it) * 1.0/number_hvectors * 3.0/dij5 * force_vector_atom; // [kJ/(nm*mol)] (some |rh| are omitted); //  [kJ/(nm*mol)]
      // return forces
      if(sim.param().rdc.mode == simulation::rdc_restr_inst ||
         sim.param().rdc.mode == simulation::rdc_restr_inst_weighted){
        const double dV_dD = force_coefficient * dflat_bottom_pot(conf_it->curr[k], it->R0, sim.param().rdc.delta); // [kJ/(nm^2*mol)]
        conf.current().force(it->i) -= dV_dD * dD_dr; // [kJ/(nm*mol)]
        conf.current().force(it->j) += dV_dD * dD_dr; // [kJ/(nm*mol)]
        DEBUG(15, "dV_dD; dD_dr: " << dV_dD << ", (" << scientific << dD_dr[0] << ", " << dD_dr[1] << ", " << dD_dr[2] << ")")
      }
      else if(sim.param().rdc.mode == simulation::rdc_restr_av ||
              sim.param().rdc.mode == simulation::rdc_restr_av_weighted){
        const double dVav_dDav = force_coefficient * dflat_bottom_pot(conf_it->av[k], it->R0, sim.param().rdc.delta); // [kJ/(nm^2*mol)]
        // in the following we set dDav_dD to either 1.0 or [1-e^(-dt/tau)].
        // The latter is correct, the former allows for using the same K as without time AV
        const double dDav_dD = sim.param().rdc.tAVfactor ? (1.0-pow(M_E, -sim.time_step_size()/sim.param().rdc.tau)) : 1.0;
        conf.current().force(it->i) -= dVav_dDav * dDav_dD * dD_dr; // [kJ/(nm*mol)]
        conf.current().force(it->j) += dVav_dDav * dDav_dD * dD_dr; // [kJ/(nm*mol)]
        DEBUG(15, "dVav_dDav; dDav_dD; dD_dr: " << dVav_dDav << ", " << dDav_dD << ", (" << scientific << dD_dr[0] << ", " << dD_dr[1] << ", " << dD_dr[2] << ")")
      }
      else if(sim.param().rdc.mode == simulation::rdc_restr_biq ||
              sim.param().rdc.mode == simulation::rdc_restr_biq_weighted){
        const double dVbq_dD = 2.0 * force_coefficient * dflat_bottom_pot(conf_it->curr[k], it->R0, sim.param().rdc.delta) * flat_bottom_pot(conf_it->av[k], it->R0, sim.param().rdc.delta); // [kJ/(nm^2*mol)]
        const double dVbq_dDav = 2.0 * force_coefficient * flat_bottom_pot(conf_it->curr[k], it->R0, sim.param().rdc.delta) * dflat_bottom_pot(conf_it->av[k], it->R0, sim.param().rdc.delta); // [kJ/(nm^2*mol)]
        // in the following we set dDav_dD to either 1.0 or [1-e^(-dt/tau)] or 0.
        // The first is simple, the second correct.
        double dDav_dD = 0;  // case 2
        if (sim.param().rdc.tAVfactor == 0) dDav_dD = 1.0;
        else if (sim.param().rdc.tAVfactor == 1) dDav_dD = 1.0-pow(M_E, -sim.time_step_size()/sim.param().rdc.tau);
        conf.current().force(it->i) -= (dVbq_dD + dVbq_dDav * dDav_dD) * dD_dr; // [kJ/(nm*mol)]
        conf.current().force(it->j) += (dVbq_dD + dVbq_dDav * dDav_dD) * dD_dr; // [kJ/(nm*mol)]
        DEBUG(15, "dVbq_dD; dVbq_dDav; dDav_dD; dD_dr: " << dVbq_dD << ", " << dVbq_dDav << ", "<< dDav_dD << ", (" << scientific << dD_dr[0] << ", " << dD_dr[1] << ", " << dD_dr[2] << ")")
      }
    } // Loop over RDC
  } // loop over rdc groups
}


template<math::boundary_enum B>
double mf_potential(const gsl_vector* MF, void* param){
  DEBUG(10, "RDC/MF/EM pot")

  // MF are the variables and param all the fixed parameters per timestep
  // typecast param from void* to fit_param*
  fit_param<B>& parameters = *static_cast<fit_param<B>*>(param);
  vector<vector<topology::rdc_restraint_struct> >::iterator topo_it = parameters.topo_it;
  vector<configuration::Configuration::special_struct::rdc_struct>::iterator conf_it = parameters.conf_it;
  configuration::Configuration & conf = *parameters.conf;
  simulation::Simulation & sim = *parameters.sim;

  // number of magnetic field vectors
  const int number_hcomponents = MF->size;
  const int number_hvectors = number_hcomponents/2.0;

  // Create Array of vectors out of MF components (which is a concatenation of all magnetic field vectors)
  math::VArray h_vectors(number_hvectors);
  for (int h=0; h<number_hvectors; ++h) {
    const double r = 1.0;
    const double theta = gsl_vector_get(MF, 2*h);
    const double phi = gsl_vector_get(MF, 2*h +1);
//    cout << "theta: " << theta << endl;
//    cout << "phi: " << phi << endl;
    h_vectors[h] = math::Vec(r * sin(theta) * cos(phi), r * sin(theta) * sin(phi), r * cos(theta));
  }

  // create variables
  math::Periodicity<B> periodicity(conf.current().box);
  math::VArray &pos = conf.current().pos;
  math::Vec ri(0.0, 0.0, 0.0);
  math::Vec rj(0.0, 0.0, 0.0);
  math::Vec rij(0.0, 0.0, 0.0);

  double potential = 0.0;

  // Loop over all RDC constraints
  int k=0;
  vector<topology::rdc_restraint_struct>::iterator
    it = topo_it->begin(),
    to = topo_it->end();
  for(; it!=to; ++it, ++k) {
    // Get positions of the two atoms involved in this RDC
    ri = pos(it->i);
    rj = pos(it->j);
    periodicity.nearest_image(ri, rj, rij);
    const double dij = math::abs(rij);
    const double dij2 = dij*dij;
    const double dij3 = dij2*dij;

    conf_it->curr[k] = 0.0;
    // Loop over all magnetic field vectors
    for(int h=0; h<number_hvectors; ++h) {
      conf_it->curr[k] +=  pow(math::dot(rij, h_vectors[h]) / dij, 2) - 1.0/3.0; // [1]
    } // Loop over MF
    conf_it->curr[k] *= RDC_max(it)/dij3 * 1.0/number_hvectors * 1.5; // [1/ps]

    const double force_coefficient = sim.param().rdc.K * it->weight * math::n_avogadro; // [kJ*ps^2/mol]

    // potential V = 0.5 * K * (RDC - RDC_0)^2
    potential += force_coefficient * flat_bottom_pot(conf_it->curr[k], it->R0, sim.param().rdc.delta); // [kJ/mol]
    //cout << "potential: " << potential << endl;
  } // loop over RDCs

  return potential;
}


template<math::boundary_enum B>
void mf_gradient(const gsl_vector* MF, void* param, gsl_vector* gradient) {
  DEBUG(10, "RDC/MF/EM grad")
  // typecast param from void* to fit_param*
  fit_param<B>& parameters = *static_cast<fit_param<B>*>(param);
  vector<vector<topology::rdc_restraint_struct> >::iterator topo_it = parameters.topo_it;
  vector<configuration::Configuration::special_struct::rdc_struct>::iterator conf_it = parameters.conf_it;
  configuration::Configuration & conf = *parameters.conf;
  simulation::Simulation & sim = *parameters.sim;

  // number of magnetic field vectors
  const int number_hcomponents = MF->size;
  const int number_hvectors = number_hcomponents/2.0;

  // Create Array of vectors with MF components
  math::VArray h_vectors(number_hvectors);
  for (int h=0; h<number_hvectors; ++h) {
    const double r = 1.0;
    const double theta = gsl_vector_get(MF, 2*h);
    const double phi = gsl_vector_get(MF, 2*h +1);
    h_vectors[h] = math::Vec(r * sin(theta) * cos(phi), r * sin(theta) * sin(phi), r * cos(theta));
  }

  // create variables
  math::Periodicity<B> periodicity(conf.current().box);

  math::VArray &pos = conf.current().pos;
  math::Vec ri(0.0, 0.0, 0.0);
  math::Vec rj(0.0, 0.0, 0.0);
  math::Vec rij(0.0, 0.0, 0.0);
  math::VArray force_vectors_hfield;

  // initialise the vector holding hfield force vectors
  for (int h=0; h<number_hvectors; ++h) {
    force_vectors_hfield.push_back(math::Vec(0.0, 0.0, 0.0));
  }

  ///////////////
  // MAIN LOOP //
  ///////////////
  int k=0;
  vector<topology::rdc_restraint_struct>::iterator
    it = topo_it->begin(),
    to = topo_it->end();
  for(; it!=to; ++it, ++k) {
    // Get positions of the two atoms involved in this RDC
    ri = pos(it->i);
    rj = pos(it->j);
    periodicity.nearest_image(ri, rj, rij);
    const double dij = math::abs(rij);
    const double dij2 = dij*dij; // [nm^2]
    const double dij3 = dij2*dij; // [nm^3]
    const double dij5 = dij2*dij3; // [nm^5]

    conf_it->curr[k] = 0.0;
    // Loop over all magnetic field vectors
    for(int h=0; h<number_hvectors; ++h) {
      conf_it->curr[k] += pow(math::dot(rij, h_vectors[h]) / dij, 2) - 1.0/3.0; // [1]
    } // Loop over MF
    conf_it->curr[k] *= RDC_max(it)/dij3 * 1.0/number_hvectors * 1.5; // [1/ps]

    const double force_coefficient = sim.param().rdc.K * it->weight * math::n_avogadro; // [kJ*ps^2/mol]

    // f = - K * (RDC - RDC_0) [ * inner derivative]
    for(int h=0; h<number_hvectors; ++h) {
      force_vectors_hfield[h] -= 3.0/dij5 * force_coefficient * dflat_bottom_pot(conf_it->curr[k], it->R0, sim.param().rdc.delta) * RDC_max(it) * 1.0/number_hvectors * math::dot(rij, h_vectors[h]) * rij; // [kJ/(nm*mol)]
    } // Loop over MF (Vphys)
//    DEBUG(12, "force on MF foreach rdc:" << math::v2s(force_vectors_hfield[0]) << ", " << math::v2s(force_vectors_hfield[1]) << ", "<< math::v2s(force_vectors_hfield[2]) )
  } // Loop over RDC

  // copy force_vector_hfield (which is [h][3]) to the gsl_vector (which is [2h]) and do a coordinate transformation at the same time
  // the returned gradient function should return the gradient, not the resulting force ... hence the '-'
  for (int h=0; h<number_hvectors; ++h) {
     const double dE_over_dx = -force_vectors_hfield[h][0];
     const double dE_over_dy = -force_vectors_hfield[h][1];
     const double dE_over_dz = -force_vectors_hfield[h][2];

     const double theta = gsl_vector_get(MF, 2*h);
     const double phi = gsl_vector_get(MF, 2*h +1);

     // five of the nine entries of the jacobian matrix [(x,y,z) over (r,theta,phi)]
     // we dont need the three entries that involve r (because r is constant and the derivatives are 0)
     // and dz_dphi is zero
     const double dx_over_dtheta =  cos(theta) * cos(phi);
     const double dy_over_dtheta =  cos(theta) * sin(phi);
     const double dz_over_dtheta = -sin(theta);
     const double dx_over_dphi   = -sin(theta) * sin(phi);
     const double dy_over_dphi   =  sin(theta) * cos(phi);

     gsl_vector_set(gradient, 2*h,    dE_over_dx*dx_over_dtheta + dE_over_dy*dy_over_dtheta + dE_over_dz*dz_over_dtheta );
     gsl_vector_set(gradient, 2*h +1, dE_over_dx*dx_over_dphi   + dE_over_dy*dy_over_dphi );
     //cout << "grad theta: " << dE_over_dx*dx_over_dtheta + dE_over_dy*dy_over_dtheta + dE_over_dz*dz_over_dtheta << endl;
     //cout << "grad phi: " << dE_over_dx*dx_over_dphi   + dE_over_dy*dy_over_dphi << endl;
  } //loop over MF vectors
}


// this is inefficient but is only called once per cycle and hence not relevant
template<math::boundary_enum B>
void mf_potential_gradient(const gsl_vector* MF, void* param, double* f, gsl_vector* gradient) {
  DEBUG(10, "RDC/MF/EM pot+grad")
  *f = mf_potential<B>(MF, param);
  mf_gradient<B>(MF, param, gradient);
}


template<math::boundary_enum B, math::virial_enum V>
int _calculate_interactions_mfield(topology::Topology & topo,
                                   configuration::Configuration & conf,
                                   simulation::Simulation & sim) {

#ifdef DEBUG
  cout.precision(14);
  setw(20);
  cout.setf(ios::fixed, ios::floatfield);
#endif

  DEBUG(5, "Mode chosen: " <<  sim.param().rdc.method)

  switch(sim.param().rdc.method) {
    case simulation::rdc_em: {
      ///////////////////////////////////////////
      //                                       //
      //     E     M I N I M I S A T I O N     //
      //                                       //
      ///////////////////////////////////////////

      const double stepsize = sim.param().rdc.emstepsize; // eg 1e-3
      const double terminate_gradient = sim.param().rdc.emgradient; // eg 1e-5
      const double tol = 1e-1; // accuracy of line minimization, the manual suggests 0.1
      const unsigned int max_iterations = sim.param().rdc.emmaxiter;

      vector<configuration::Configuration::special_struct::rdc_struct>::iterator
          conf_it = conf.special().rdc.begin(),
          conf_to = conf.special().rdc.end();
      vector<vector<topology::rdc_restraint_struct> >::iterator
          topo_it = topo.rdc_restraints().begin();
      for( ; conf_it!=conf_to; conf_it++, topo_it++){

        const int number_magnetic_field_vectors = conf_it->MFpoint.size();
        const int number_magnetic_field_vector_components = 2*number_magnetic_field_vectors;
        // Create the parameters-struct, fill in all needed variabled for the potential function and the derivative
        fit_param<B> parameters;
        parameters.topo_it = topo_it;
        parameters.conf_it = conf_it;
        parameters.conf = &conf;
        parameters.sim = &sim;

        // Create the minimisation function block, define the different functions and parameters
        gsl_multimin_function_fdf potential_function;
        potential_function.n = number_magnetic_field_vector_components;
        potential_function.f = &mf_potential<B>;
        potential_function.df = &mf_gradient<B>;
        potential_function.fdf = &mf_potential_gradient<B>;
        potential_function.params = (void *) &parameters;

        // Starting point
        gsl_vector *MF = gsl_vector_alloc(number_magnetic_field_vector_components);
        for(int h=0; h<number_magnetic_field_vectors; ++h) {
          const double x = conf_it->MFpoint[h][0];
          const double y = conf_it->MFpoint[h][1];
          const double z = conf_it->MFpoint[h][2];
          const double r = sqrt(x*x + y*y + z*z);
          assert(abs(1.0 - r) < 1.0e-5);
          const double theta = acos(z);
          const double phi = atan2(y, x);
          // cout << "x: " << x << "y: " << y << "z: " << z << endl;
          // cout << "theta: " << theta << "phi: " << phi << endl;
          gsl_vector_set(MF, 2*h,    theta);
          gsl_vector_set(MF, 2*h +1, phi);
        }

        size_t iter = 0;
        int status = 0;
        const gsl_multimin_fdfminimizer_type *fT = gsl_multimin_fdfminimizer_conjugate_fr;
        gsl_multimin_fdfminimizer *s = gsl_multimin_fdfminimizer_alloc(fT, number_magnetic_field_vector_components);
        gsl_multimin_fdfminimizer_set(s, &potential_function, MF, stepsize, tol); // runs fdf once
        do {
          ++iter;
          status = gsl_multimin_fdfminimizer_iterate(s);
  //        cout << "EM Loop #" << iter << "  status: " << status << endl;
  //        cout << "MF Vectors: " << scientific << setprecision(6) << endl;
  //        for(int i=0; i<number_magnetic_field_vectors; ++i) {
  //          cout << gsl_vector_get(s->x, 2*i) << " "
  //                    << gsl_vector_get(s->x, 2*i + 1) << endl;
  //        }
  //        cout << "Potential: " << s->f << endl
  //                  << "Gradient: " << scientific << setprecision(6) << endl;
  //        for(int i=0; i<number_magnetic_field_vectors; ++i) {
  //          cout << gsl_vector_get(s->gradient, 2*i)     << " "
  //                    << gsl_vector_get(s->gradient, 2*i + 1) << endl
  //                    << "Norm of Gradient: " << sqrt(
  //                      (gsl_vector_get(s->gradient, 2*i)    * gsl_vector_get(s->gradient, 2*i)) +
  //                      (gsl_vector_get(s->gradient, 2*i +1) * gsl_vector_get(s->gradient, 2*i +1)) )
  //                    << endl;
  //        }

          if(status) {
            break;
          }

#ifdef DEBUG
          for (int i =0; i<number_magnetic_field_vector_components; ++i){
            cout << "grad[" << i << "]: " << gsl_vector_get(s->gradient, i);
          }
          cout << endl;
#endif // DEBUG

          status = gsl_multimin_test_gradient(s->gradient, terminate_gradient);
          //cout << "Status 2: " << status << endl;

        } while (status == GSL_CONTINUE && iter < max_iterations);

#ifndef FINITE_DIFF
        // 1. return magnetic field vectors
        for(int h=0; h<number_magnetic_field_vectors; ++h) {
          const double r = 1.0; // [nm]
          const double theta = gsl_vector_get(s->x, 2*h);
          const double phi = gsl_vector_get(s->x, 2*h +1);
          conf_it->MFpoint.at(h)[0] = r * sin(theta) * cos(phi); // [nm]
          conf_it->MFpoint.at(h)[1] = r * sin(theta) * sin(phi); // [nm]
          conf_it->MFpoint.at(h)[2] = r * cos(theta); // [nm]
          DEBUG(10, "mf[" << h << "]: " << r * sin(theta) * cos(phi) << ", " << r * sin(theta) * sin(phi) << ", " << r * cos(theta))
        }
#endif // FINITE_DIFF

        gsl_multimin_fdfminimizer_free(s);
        gsl_vector_free(MF);

#ifdef SCRIPT
        double denominator=0.0, enumerator=0.0;
        int k=0;
        vector<topology::rdc_restraint_struct>::iterator
          it = topo_it->begin(),
          to = topo_it->end();
        for(; it!=to; ++it, ++k) {
          enumerator += pow(conf_it->curr[k] - it->R0, 2);
          denominator += pow(it->R0, 2);
        }
        double Q = sqrt(enumerator / denominator);
        cout << "SCRIPT " << Q << endl;
#endif

      } // loop over rdc groups

      // 2. write energy of rdc interaction
      // 3. forces due to rdc interaction
      // forces on all atoms are not zero (but their sum is)
      // the following writes back interaction energies and forces for each RDC
      _calculate_forces_atoms_MF<B>(topo, conf, sim);

      break;
    }


    case simulation::rdc_sd: {
      ///////////////////////////////////////////
      //                                       //
      // S T O C H A S T I C   D Y N A M I C S //
      //                                       //
      ///////////////////////////////////////////
      DEBUG(7, "doing RDC/MF/SD")

      // create variables
      math::VArray force_vectors_hfield;

      _calculate_forces_vectors_MF<B>(topo, conf, sim, force_vectors_hfield);

      /////////////////////////
      //      LEAP FROG      //
      /////////////////////////

      vector<configuration::Configuration::special_struct::rdc_struct>::iterator
          conf_it = conf.special().rdc.begin(),
          conf_to = conf.special().rdc.end();
      vector<vector<topology::rdc_restraint_struct> >::iterator
          topo_it = topo.rdc_restraints().begin();
      for( ; conf_it!=conf_to; conf_it++, topo_it++){

      // number of magnetic field vectors
      const int number_hvectors = conf_it->MFpoint.size();

      const double gamma = sim.param().rdc.sdfric;
      const double gamma_abs = fabs(gamma);
      DEBUG(12, "friction coeff: " << gamma)
      DEBUG(12, "step length: " << sim.time_step_size())
      DEBUG(12, "Tref: " << sim.param().rdc.temp)

// FIXME put the rng in conf_it->rng ?
#ifdef DEBUG
      // Create random number generator, set seed to 0 (to make sure it's reproducible)
      static math::RandomGenerator* rng = math::RandomGenerator::create(sim.param(), "0");
#else
      // Create random number generator, choose a remotely random random-seed
      ostringstream ss;
      ss << time(NULL);
      static math::RandomGenerator* rng = math::RandomGenerator::create(sim.param(), ss.str());
#endif // DEBUG

      // Save old velocities and positions
      const math::VArray old_pos(conf_it->MFpoint);
      const math::VArray old_vel(conf_it->MFpointVel);

      vector<double> kT_over_m(number_hvectors);
      for(int h=0; h<number_hvectors; ++h){
        kT_over_m[h] = sqrt(math::k_Boltzmann * sim.param().rdc.temp / conf_it->MFpointMass[h]);
      }

      // Calculate SD coefficients
      const double gdt = gamma * sim.time_step_size();
      const double gdth = gdt * 0.5;

      // depending on the value of gth we choose different expressions to ensure a high accuracy
      double c1 = 0.0, c2 = 0.0, c3 = 0.0, c4 = 0.0, c5 = 0.0, c6 = 0.0, c7 = 0.0, c8 = 0.0, c9 = 0.0;
      if(fabs(gdt) > 0.05){
        DEBUG(10, "Doing the analytical formulas, (gamma * step_size = "<< gdt << " > 0.05)")

        const double emdth = exp(-gdth);
        const double epdth = exp(+gdth);
        const double emdt  = emdth * emdth;
        const double epdt  = epdth * epdth;
        const double omdt  = 1.0 - emdt;
        const double cdth  = gdt - 3.0 + 4.0 * emdth - emdt;
        const double ddth  = 2.0 - epdth - emdth;
        const double bpdth = gdt * (epdt - 1.0) - 4.0 * (epdth - 1.0) * (epdth - 1.0);
        const double bmdth = gdt * (1.0 - emdt) - 4.0 * (emdth - 1.0) * (emdth - 1.0);

        c1 = emdt;
        c2 = omdt / gdt;
        c3 = sqrt(fabs(omdt));
        c4 = sqrt(fabs(bpdth/cdth));
        c5 = gamma_abs * ddth/cdth;
        c6 = (epdth - emdth) / gdt;
        c7 = sqrt(fabs(cdth)) / gamma_abs;
        c8 = sqrt(fabs(bmdth/omdt)) / gamma_abs;
        c9 = -ddth/(gamma_abs * omdt);

        if(sim.steps() == 0) { // and we do all this to make sure the stochastic integral is initialized
          //DEBUG(10, "sim.steps(): " << sim.steps())

          for (int h=0; h<number_hvectors; ++h){
            const double sd = kT_over_m[h] / gamma * sqrt(cdth);
            rng->stddev(sd);
            conf_it->stochastic_integral_mf(h) = rng->get_gaussian_vec();
          }
        }
      } // calculate coefficients using the analytoc formula
      else {
        DEBUG(10, "doing the power series, (gamma * step_size = " << gdt << " < 0.05)")

        //Power Series
        const double gdth2 = gdth * gdth;
        const double gdth3 = gdth2 * gdth;
        const double gdth4 = gdth2 * gdth2;
        const double gdth5 = gdth2 * gdth3;
        const double gdth6 = gdth3 * gdth3;
        const double gdth7 = gdth4 * gdth3;

        c1 = exp(-gdt);
        c2 = 1.0 - gdth + gdth2 * 2.0/3.0 - gdth3/3.0 + gdth4 * 2.0/15.0 - gdth5 * 2.0/45.0 + gdth6 * 4.0/315.0;
        c3 = sqrt(fabs(c2) * 2.0 * gdth);
        c4 = sqrt(fabs(gdth/2.0 + gdth2 * 7.0/8.0 + gdth3 * 367.0/480.0 + gdth4 * 857.0/1920.0 + gdth5 * 52813.0/268800.0 + gdth6 * 224881.0/3225600.0 + gdth7 * 1341523.0/64512000.0));
        c5 = -2.0 / sim.time_step_size() * (1.5 + gdth * 9.0/8.0 + gdth2 * 71.0/160.0 + gdth3 * 81.0/640.0 + gdth4 * 7807.0/268800.0 + gdth5 * 1971.0/358400.0 + gdth6 * 56417.0/64512000.0);
        c6 = 1.0 + gdth2/6.0 + gdth4/10.0 + gdth6/5040.0;
        c7 = sim.time_step_size() * 0.5 * sqrt(fabs(gdth * 2.0/3.0 - gdth2/2.0 + gdth3 * 7.0/30.0 -gdth4/12.0 + gdth5 * 31.0/1260.0 - gdth6/160.0 + gdth7 * 127.0/90720.0));
        c8 = sim.time_step_size() * 0.5 * sqrt(fabs(gdth/6.0 - gdth3/60.0 + gdth5 * 17.0/10080.0 - gdth7 * 31.0/181440.0));
        c9 = sim.time_step_size() * 0.5 * (0.5 + gdth/2.0 + gdth2 * 5.0/24.0 + gdth3/24.0 + gdth4/240.0 + gdth5/720.0 + gdth6 * 5.0/8064.0);

        if(sim.steps() == 0){ // and we do all this to make sure the stochastic integral is initialized
          for (int h=0; h<number_hvectors; ++h){
            if (gamma < math::epsilon) {
              conf_it->stochastic_integral_mf(h) = 0.0;
            }
            else {
              const double emdth = exp(-gdth);
              const double emdt  = emdth * emdth;
              const double cdth  = gdt - 3.0 + 4.0 * emdth - emdt;

              const double sd = kT_over_m[h] / gamma * sqrt(cdth);
              rng->stddev(sd);
              conf_it->stochastic_integral_mf(h) = rng->get_gaussian_vec();
            }
          }
        }
      } // Power series

      DEBUG(10, "c1 through c9: " << c1 << ", " << c2 << ", " << c3 << ", " << c4 << ", " << c5 << ", " << c6 << ", " << c7 << ", " << c8 << ", " << c9)

      conf_it->Ekin = 0.0;

#ifndef FINITE_DIFF
      for(int h=0; h<number_hvectors; ++h){
        // Vel1
        DEBUG(10, "doing stochastic dynamics velocities (1)")
        math::Vec vrand1, vrand2, vrand3, vrand4;
        //2.11.2.16 (rho1 squared)
        const double sd1 = c3 * kT_over_m[h];
        //2.11.2.9 (sigma2 squared)
        const double sd2 = c4 * kT_over_m[h];
        //2.11.2.8 (sigma1 squared)
        const double sd3 = c7 * kT_over_m[h];
        //2.11.2.17 (rho2 squared)
        const double sd4 = c8 * kT_over_m[h];
        // Sample the V' vector from 2.11.2.20 from Gaussian with mean=0.0, sd=sd1
        rng->stddev(sd1);
        vrand1 = rng->get_gaussian_vec();
        // Sample the V vector from 2.11.2.2 from Gaussian with mean=0.0, sd=sd2
        rng->stddev(sd2);
        vrand2 = rng->get_gaussian_vec();
        // Sample the R' vector from 2.11.2.25 from Gaussian with mean=0.0, sd=sd3
        rng->stddev(sd3);
        vrand3 = rng->get_gaussian_vec();
        // Sample the R vector from 2.11.2.26 from Gaussian with mean=0.0, sd=sd4
        rng->stddev(sd4);
        vrand4 = rng->get_gaussian_vec();
        DEBUG(10, "vrand1=" << math::v2s(vrand1))
        DEBUG(10, "vrand2=" << math::v2s(vrand2))
        DEBUG(10, "vrand3=" << math::v2s(vrand3))
        DEBUG(10, "vrand4=" << math::v2s(vrand4))

        //v1, p1, shake, v2, p2, shake

        // Update the velocity
        const math::Vec svh = conf_it->stochastic_integral_mf(h) * c5 + vrand2; // FIXME what is this?
        const double cf = sim.time_step_size() * c2 / conf_it->MFpointMass[h]; // FIXME ditto
        // Assign vrand1 to the stochastic integrals
        conf_it->stochastic_integral_mf(h) = vrand1;
        // 2.11.2.2 (2.11.2.3 minus first term of the third line)
        DEBUG(10, "doing stochastic dynamics: velocities (1)")
        conf_it->MFpointVel[h] = ((conf_it->MFpointVel[h] - svh) * c1)
            + (force_vectors_hfield[h] * cf)
            + vrand1;

        // Pos1
        DEBUG(10, "doing stochastic dynamics: positions (1)")
        DEBUG(10, "old pos(" << h << ") " << math::v2s(conf_it->MFpoint[h]))
        conf_it->MFpoint[h] += conf_it->MFpointVel[h] * sim.time_step_size() * c6;
        DEBUG(10, "new pos(" << h << ") " << math::v2s(conf_it->MFpoint[h]))

        // there should be SHAKE here!
        // instead we just rescale, because it's cheap and the system is simple enough as to not require recursion
        DEBUG(10, "SHAKE(rescale) on position (1)")
        double length_after_displacement = sqrt(pow(conf_it->MFpoint[h][0],2) +  pow(conf_it->MFpoint[h][1],2) + pow(conf_it->MFpoint[h][2],2) );
        DEBUG(10, "length after displacement 1: " << length_after_displacement)
        DEBUG(10, "MFpoint before rescaling: " << math::v2s(conf_it->MFpoint[h]))
        conf_it->MFpoint[h] = conf_it->MFpoint[h] / length_after_displacement;
        DEBUG(10,  "MFpoint afterrescaling: " << math::v2s(conf_it->MFpoint[h]));

        DEBUG(10, "doing stochastic dynamics: velocities (2)")
        // Velocity correction 2.11.2.24
        const double cinv = 1.0 / (c6 * sim.time_step_size());
        const math::Vec r = conf_it->MFpoint[h] - old_pos[h];
        conf_it->MFpointVel[h] = r * cinv;

        DEBUG(10, "doing stochastic dynamics: positions (2)")
        //2.11.2.25
        const math::Vec sxh = conf_it->stochastic_integral_mf(h) * c9 + vrand4;
        conf_it->stochastic_integral_mf(h) = vrand3;
        //2.11.2.26
        conf_it->MFpoint[h] += (vrand3 - sxh);

        // there should be SHAKE here!
        // instead we just rescale, because it's cheap and the system is simple enough as to not require recursion
        DEBUG(10, "SHAKE(rescale) on position (2)")
        length_after_displacement = sqrt(pow(conf_it->MFpoint[h][0],2) +  pow(conf_it->MFpoint[h][1],2) + pow(conf_it->MFpoint[h][2],2) );
        DEBUG(10, "length after displacement 1: " << length_after_displacement)
        DEBUG(10, "MFpoint before rescaling: " << math::v2s(conf_it->MFpoint[h]))
        conf_it->MFpoint[h] /=  length_after_displacement;
        DEBUG(10, "MFpoint after rescaling: " << math::v2s(conf_it->MFpoint[h]))

        DEBUG(10, "SHAKE(rescale) on velocities (1)")
        DEBUG(10, "vel before shake: " << math::v2s(conf_it->MFpointVel[h]))
        conf_it->MFpointVel[h] = (conf_it->MFpoint[h] - old_pos[h]) / sim.time_step_size(); // [nm/ps]
        DEBUG(10, "vel after shake: " << math::v2s(conf_it->MFpointVel[h]))

        // 'Ekin' of the MF vector
        // E = 0.5 * m * v^2 = 0.25 * m * (v_{n-.5}^2 + v_{n+.5}^2)
        // m[u = g/mol], v[nm/ps]
        conf_it->Ekin += 0.25 * conf_it->MFpointMass[h] * (math::abs2(old_vel[h]) + math::abs2(conf_it->MFpointVel[h])); // [kJ/mol]
        DEBUG(10, "Ekin of (" << h << "): " << (math::abs2(old_vel[h]) + math::abs2(conf_it->MFpointVel[h]) ) * 0.25 * conf_it->MFpointMass[h])

      } // Loop over MF
#endif // FINITE_DIFF

#ifdef SCRIPT
      double denominator=0.0, enumerator=0.0;
      int k=0;
      vector<topology::rdc_restraint_struct>::iterator
        it = topo_it->begin(),
        to = topo_it->end();
      for(; it!=to; ++it, ++k) {
        enumerator += pow(conf_it->curr[k] - it->R0, 2);
        denominator += pow(it->R0, 2);
      }
      double Q = sqrt(enumerator / denominator);
      cout << "SCRIPT " << Q << endl;
#endif

      } // loop over rdc groups

      // 1. return potential energy of interaction
      // 2. return forces on atoms
      // the following writes back interaction energies and forces for each RDC
      _calculate_forces_atoms_MF<B>(topo, conf, sim);
      // 3. Ekin is written to the rdc-struct
      // 4. return coordinates of vectors (done implicitly)
      // 5. return velocities (done implicitly)
      // 6. return random state //FIXME

      break;
    }

    case simulation::rdc_md: {
      /////////////////////////////////////////
      //                                     //
      // M O L E C U L A R   D Y N A M I C S //
      //                                     //
      /////////////////////////////////////////

      // create variables
      math::VArray force_vectors_hfield;

      _calculate_forces_vectors_MF<B>(topo, conf, sim, force_vectors_hfield);

      math::Vec vold(0.0, 0.0, 0.0);
      math::Vec rold(0.0, 0.0, 0.0);

      ///////////////
      // LEAP FROG //
      ///////////////

      vector<configuration::Configuration::special_struct::rdc_struct>::iterator
          conf_it = conf.special().rdc.begin(),
          conf_to = conf.special().rdc.end();
      vector<vector<topology::rdc_restraint_struct> >::iterator
          topo_it = topo.rdc_restraints().begin();
      for( ; conf_it!=conf_to; conf_it++, topo_it++){

        // number of magnetic field vectors
        const int number_hvectors = conf_it->MFpoint.size();

        conf_it->Ekin = 0.0;

#ifndef FINITE_DIFF
        for(int h=0; h<number_hvectors; ++h) {
          // Save old values
          vold = conf_it->MFpointVel[h]; // [nm/ps]
          rold = conf_it->MFpoint[h]; // [nm]

          const double vector_length = 1.0;
          const double vector_length_squared = 1.0;

          // Velocity step v = v + f/m*\Delta t // f[kJ/(nm*mol)], dt[ps], m[u]
          // 1u = M_u/N_A = g/(6*10^23) && M_u = 1g/mol && N_A = 6*10^23 1/mol
          // but slightly handwavy we can use u = g/mol
          conf_it->MFpointVel[h] += force_vectors_hfield[h] * sim.time_step_size() / conf_it->MFpointMass[h]; // [nm/ps]
          // Position step r = r + v*\Delta t // r[nm], v[nm/ps], dt[ps]
          conf_it->MFpoint[h] += conf_it->MFpointVel[h] * sim.time_step_size(); // [nm]

          // Apply SHAKE to constrain the length of the vectors to length 1
          //constraint correction
          //TODO check if this is appropriate
          // delta = -[r_new * r_old - 1]*rold
          //        conf_it->MFpoint[h] -= (math::dot(conf_it->MFpoint[h], rold) - 1) * rold;
          const double discriminant = 4 * pow(math::dot(conf_it->MFpoint[h], rold),2) - 4 * math::abs2(rold) * (math::abs2(conf_it->MFpoint[h]) - vector_length_squared);
          //        cout << "disc: " << discriminant << endl;
          if (discriminant >= 0.0){// small displacement, we can use SHAKE (something SHAKE-like anyway ...)
            conf_it->MFpoint[h] -= (2 * math::dot(conf_it->MFpoint[h], rold) - sqrt(discriminant)) / (2 * math::abs2(rold)) * rold;
          }else{// larger displacement, we rescale in direction of r(t+\delta t) which is cheap and works for any length of r(t+\delta t)
            conf_it->MFpoint[h] = vector_length / math::abs(conf_it->MFpoint[h]) * conf_it->MFpoint[h];
            DEBUG(10, "something went wrong:  very large displacement of MFpoint (requiring very large corrections)")
          }
  //        cout <<"abs: " << math::abs(conf_it->MFpoint[h]) << endl;

          // calculate the new (constrained) velocities with respect to the constrained position
          conf_it->MFpointVel[h] = (conf_it->MFpoint[h] - rold) / sim.time_step_size(); // [nm/ps]

          // 'Ekin' of the MF
          // E = 0.5 * m * v^2 = 0.25 * m * (v_{n-.5}^2 + v_{n+.5}^2)
          // m[u = g/mol], v[nm/ps]
          conf_it->Ekin += 0.25 * conf_it->MFpointMass[h] * (math::abs2(vold) + math::abs2(conf_it->MFpointVel[h])); // [kJ/mol]
        } // Loop over MF
#endif // FINITE_DIFF

#ifdef SCRIPT
        double denominator=0.0, enumerator=0.0;
        int k=0;
        vector<topology::rdc_restraint_struct>::iterator
          it = topo_it->begin(),
          to = topo_it->end();
        for(; it!=to; ++it, ++k) {
          enumerator += pow(conf_it->curr[k] - it->R0, 2);
          denominator += pow(it->R0, 2);
        }
        double Q = sqrt(enumerator / denominator);
        cout << "SCRIPT " << Q << endl;
#endif

      } // loop over rdc groups

      // 1. return potential energy of interaction
      // 2. return forces on atoms
      // the following writes back interaction energies and forces for each RDC
      _calculate_forces_atoms_MF<B>(topo, conf, sim);
      // 3. Ekin is written to the rdc-struct
      // 4. return coordinates of vectors (done implicitly)
      // 5. return velocities (done implicitly)

      break;
    }
    default: {
      // should never be the case
      assert(false);
    }
  }
  return 0;
}




//************************************************************************
//************************************************************************
//*********************                           ************************
//*********************     T  E  N  S  O  R      ************************
//*********************                           ************************
//************************************************************************
//************************************************************************


// calculate optimal ah
template<math::boundary_enum B>
void _calculate_ah(topology::Topology & topo,
                   configuration::Configuration & conf,
                   const simulation::Simulation & sim) {

#ifdef DEBUG
  cout.precision(14);
  setw(20);
  cout.setf(ios::fixed, ios::floatfield);
#endif

  // create variables
  math::Periodicity<B> periodicity(conf.current().box);
  math::VArray &pos = conf.current().pos;
  math::Vec ri(0.0, 0.0, 0.0);
  math::Vec rj(0.0, 0.0, 0.0);
  math::Vec rij(0.0, 0.0, 0.0);

  const int n_ah = 5;

  vector<configuration::Configuration::special_struct::rdc_struct>::iterator
      conf_it = conf.special().rdc.begin(),
      conf_to = conf.special().rdc.end();
  vector<vector<topology::rdc_restraint_struct> >::iterator
      topo_it = topo.rdc_restraints().begin();
  for( ; conf_it!=conf_to; conf_it++, topo_it++){

    const int n_rdc = topo_it->size();
    //cout << "n_rdc in calculate ah : " << n_rdc << endl;

    double ck[n_rdc][n_ah];
    for(int k=0; k<n_rdc; ++k){for(int h=0; h<n_ah; ++h){ck[k][h]=0.0;}} // init
    double matrix_array[n_ah*n_ah];
    for(int h=0; h<n_ah*n_ah; ++h){matrix_array[h]=0.0;} // init
    double result_array[n_ah];
    for(int h=0; h<n_ah; ++h){result_array[h]=0.0;} // init

    vector<topology::rdc_restraint_struct>::iterator
        it = topo_it->begin(),
        to = topo_it->end();
    int k=0; // counter of rdcs
    // Loop over all RDC constraints
    for(; it!=to; ++it, ++k){
      // Get positions of the two atoms involved in this RDC, calculate d
      ri = pos(it->i);
      rj = pos(it->j);
      periodicity.nearest_image(ri, rj, rij);
      const double dij = math::abs(rij); // [nm]
      const double dij2 = pow(dij,2); // [nm^2]
      const double dij3 = pow(dij,3); // [nm^3]
      const double x = (ri - rj)[0], y = (ri - rj)[1], z = (ri - rj)[2];

      // Calculate the five tensor factors
      ck[k][0] = (x*x - z*z)/dij2; // [1/nm] ... pseudo nm
      ck[k][1] = (y*y - z*z)/dij2; // [1/nm] ... pseudo nm
      ck[k][2] = (2*x*y)/dij2; // [1/nm] ... pseudo nm
      ck[k][3] = (2*x*z)/dij2; // [1/nm] ... pseudo nm
      ck[k][4] = (2*y*z)/dij2; // [1/nm] ... pseudo nm
      for(int h=0; h<n_ah; ++h){
        DEBUG(15, "ck (" << k << "," << h << "): " << scientific << ck[k][h])
      }

      for(int h_prime=0; h_prime<n_ah; ++h_prime){
        for(int h=0; h<n_ah; ++h){
          matrix_array[h_prime*n_ah + h] += it->weight * ck[k][h] * ck[k][h_prime]; // [1/nm^2]
        }
        result_array[h_prime] += it->weight * (it->R0 * dij3 / RDC_max(it)) * ck[k][h_prime]; // [1/nm]
      }
    } // iterate over rdcs

    for(int i=0; i<n_ah; ++i){
      for(int j=0; j<n_ah; ++j){
        DEBUG(15, "matrix[" << i*n_ah + j << "]: " << scientific << matrix_array[i*n_ah + j])
      }
    }
    for(int j=0; j<n_ah; ++j){
      DEBUG(15, "results vector[" << j << "]: " << scientific << result_array[j])
    }

    // solve the SLE using gsl
    gsl_matrix_view m = gsl_matrix_view_array (matrix_array, n_ah, n_ah);
    gsl_vector_view b = gsl_vector_view_array (result_array, n_ah);
    gsl_vector *x = gsl_vector_alloc (n_ah);
    int signum = 0; // not used
    gsl_permutation * p = gsl_permutation_alloc (n_ah);

    gsl_linalg_LU_decomp (&m.matrix, p, &signum);
    //switch off and save error handler
    gsl_error_handler_t* old_handler = gsl_set_error_handler_off();
    int status = gsl_linalg_LU_solve (&m.matrix, p, &b.vector, x);
    gsl_set_error_handler(old_handler);

    if(status == GSL_SUCCESS){ // copy back a values
      for(int h=0; h<n_ah; ++h){
        conf_it->Tensor[h] = gsl_vector_get(x, h); // [nm]
        DEBUG(15, "a[" << h << "]: " << conf_it->Tensor[h])
      }
    }else if(status == GSL_EDOM){
      // ignore singular matrices which could very rarely occur
      // do not update a_h
      io::messages.add("Singular matrix encountered in gsl_linalg_LU_solve.  This may be harmless.", "rdc_restraint_interaction", io::message::warning);
    }else{
      io::messages.add("Error in gsl_linalg_LU_solve", "rdc_restraint_interaction", io::message::error);
    }

    gsl_permutation_free (p);
    gsl_vector_free (x);

  } // loop over rdc groups
}


// calculate forces on the tensor components
template<math::boundary_enum B>
void _calculate_forces_tensor_T(topology::Topology & topo,
                   configuration::Configuration & conf,
                   const simulation::Simulation & sim,
                   vector<double> &force_array_tensor /* contains all zeros at this point*/) {

  // number of tensor components
  const int n_ah = 5;

  // create variables
  math::Periodicity<B> periodicity(conf.current().box);
  math::VArray &pos = conf.current().pos;
  math::Vec ri(0.0, 0.0, 0.0);
  math::Vec rj(0.0, 0.0, 0.0);
  math::Vec rij(0.0, 0.0, 0.0);

  vector<configuration::Configuration::special_struct::rdc_struct>::iterator
      conf_it = conf.special().rdc.begin(),
      conf_to = conf.special().rdc.end();
  vector<vector<topology::rdc_restraint_struct> >::iterator
      topo_it = topo.rdc_restraints().begin();
  for( ; conf_it!=conf_to; conf_it++, topo_it++){

    int k=0;
    vector<topology::rdc_restraint_struct>::iterator
        it = topo_it->begin(),
        to = topo_it->end();
    // Loop over all RDC constraints
    for(; it!=to; ++it, ++k){
      // Get positions of the two atoms involved in this RDC, calculate d, x, y, z
      ri = pos(it->i);
      rj = pos(it->j);
      periodicity.nearest_image(ri, rj, rij);
      const double dij = math::abs(rij);
      const double dij2 = dij*dij;
      const double dij3 = dij2*dij;
      const double x = (ri - rj)[0], y = (ri - rj)[1], z = (ri - rj)[2];

      double ck[n_ah];
      // Calculate the five tensor factors
      ck[0] = (x*x - z*z)/dij2; // [1]
      ck[1] = (y*y - z*z)/dij2; // [1]
      ck[2] = (2*x*y)/dij2; // [1]
      ck[3] = (2*x*z)/dij2; // [1]
      ck[4] = (2*y*z)/dij2; // [1]

      conf_it->curr[k] = 0.0;
      for(int i=0; i<n_ah; ++i) {
        conf_it->curr[k] += ck[i] * conf_it->Tensor[i]; // [1]
      } // Loop over tensor components
      // multiply with RDC_max to get current RDC value
      conf_it->curr[k] *= RDC_max(it) / dij3; // [1/ps]

      // locally calculate average RDC values, use for forces and don't export
      // the average RDC is not updated because that is done when forces on atoms are computed
      double local_average = 0.0;
      if(sim.param().rdc.mode == simulation::rdc_restr_av ||
         sim.param().rdc.mode == simulation::rdc_restr_av_weighted ||
         sim.param().rdc.mode == simulation::rdc_restr_biq ||
         sim.param().rdc.mode == simulation::rdc_restr_biq_weighted){
        const double exp = pow(M_E, -sim.time_step_size()/sim.param().rdc.tau); // [1]
        local_average = exp * conf_it->av[k] + (1.0 - exp) * conf_it->curr[k]; // [1/ps]
        DEBUG(15, "Delta t, tau, e^(- Delta t/tau): " << scientific << sim.time_step_size() << ", " <<  sim.param().rdc.tau << ", " << pow(M_E, -sim.time_step_size()/sim.param().rdc.tau) )
        DEBUG(10, "interaction " << k << ":  R-inst, R-av: " << conf_it->curr[k] << ", " << local_average )
      }

      const double force_coefficient = sim.param().rdc.K * it->weight * math::n_avogadro; // [kJ*ps^2/mol]

      // forces on the "tensor"
      // f = - K * (RDC - RDC_0) [ * inner derivative]
      for(int i=0; i<n_ah; ++i) {
        // force_array_tensor is summed over all RDCs
        const double dD_dck = RDC_max(it)/dij3 * ck[i]; //  [kJ/(nm*mol)]
        // forces on clm depending on AV mode
        if(sim.param().rdc.mode == simulation::rdc_restr_inst ||
           sim.param().rdc.mode == simulation::rdc_restr_inst_weighted){
          const double dV_dD = force_coefficient * dflat_bottom_pot(conf_it->curr[k], it->R0, sim.param().rdc.delta); // [kJ/(nm^2*mol)]
          force_array_tensor[i] -= dV_dD * dD_dck; // [kJ/(nm*mol)]
          DEBUG(15, "dV_dD; dD_dc: " << dV_dD << ", " << scientific << dD_dck)
        }
        else if(sim.param().rdc.mode == simulation::rdc_restr_av ||
                sim.param().rdc.mode == simulation::rdc_restr_av_weighted){
          const double dVav_dDav = force_coefficient * dflat_bottom_pot(local_average, it->R0, sim.param().rdc.delta); // [kJ/(nm^2*mol)]
          // in the following we set dDav_dD to either 1.0 or [1-e^(-dt/tau)].
          // The latter is correct, the former allows for using the same K as without time AV
          const double dDav_dD = sim.param().rdc.tAVfactor ? (1.0-pow(M_E, -sim.time_step_size()/sim.param().rdc.tau)) : 1.0;
          force_array_tensor[i] -= dVav_dDav * dDav_dD * dD_dck; // [kJ/(nm*mol)]
          DEBUG(15, "dVav_dDav; dDav_dD; dD_dc: " << dVav_dDav << ", " << dDav_dD << ", " << scientific << dD_dck)
        }
        else if(sim.param().rdc.mode == simulation::rdc_restr_biq ||
                sim.param().rdc.mode == simulation::rdc_restr_biq_weighted){
          const double dVbq_dD = 2.0 * force_coefficient * dflat_bottom_pot(conf_it->curr[k], it->R0, sim.param().rdc.delta) * flat_bottom_pot(local_average, it->R0, sim.param().rdc.delta); // [kJ/(nm^2*mol)]
          const double dVbq_dDav = 2.0 * force_coefficient * flat_bottom_pot(conf_it->curr[k], it->R0, sim.param().rdc.delta) * dflat_bottom_pot(local_average, it->R0, sim.param().rdc.delta); // [kJ/(nm^2*mol)]
          // in the following we set dDav_dD to either 1.0 or [1-e^(-dt/tau)] or 0.
          // The first is simple, the second correct.
          double dDav_dD = 0;  // case 2
          if (sim.param().rdc.tAVfactor == 0) dDav_dD = 1.0;
          else if (sim.param().rdc.tAVfactor == 1) dDav_dD = 1.0-pow(M_E, -sim.time_step_size()/sim.param().rdc.tau);
          force_array_tensor[i] -= (dVbq_dD + dVbq_dDav * dDav_dD) * dD_dck; // [kJ/(nm*mol)]
          DEBUG(15, "dVbq_dD; dVbq_dDav; dDav_dD; dD_dc: " << dVbq_dD << ", " << dVbq_dDav << ", "<< dDav_dD << ", " << scientific << dD_dck)
        }
        DEBUG(10, "tensor components(" << i << "): " << conf_it->Tensor[i])
        DEBUG(10, "forces on the tensor(" << i << "): " << force_array_tensor[i])
        DEBUG(10, "velocity of tensor component" << i << "): " << conf_it->TensorVel[i])
      } // Loop over tensor components
    } // Loop over RDCs
  } // loop over rdc groups

  //write_tensor_forces(force_array_tensor, "trft");

}


// calculate RDC interaction and forces on atoms in the tensor representation
template<math::boundary_enum B>
void _calculate_forces_atoms_T(topology::Topology & topo,
                   configuration::Configuration & conf,
                   const simulation::Simulation & sim) {

  const int n_ah = 5;

  // create variables
  math::Periodicity<B> periodicity(conf.current().box);
  math::VArray &pos = conf.current().pos;
  math::Vec ri(0.0, 0.0, 0.0);
  math::Vec rj(0.0, 0.0, 0.0);
  math::Vec rij(0.0, 0.0, 0.0);

  vector<configuration::Configuration::special_struct::rdc_struct>::iterator
      conf_it = conf.special().rdc.begin(),
      conf_to = conf.special().rdc.end();
  vector<vector<topology::rdc_restraint_struct> >::iterator
      topo_it = topo.rdc_restraints().begin();
  for( ; conf_it!=conf_to; conf_it++, topo_it++){

    int k=0;
    vector<topology::rdc_restraint_struct>::iterator
        it = topo_it->begin(),
        to = topo_it->end();
    // Loop over all RDC constraints
    for(; it!=to; ++it, ++k){
      // Get positions of the two atoms involved in this RDC, calculate dij, x, y, z
      ri = pos(it->i);
      rj = pos(it->j);
      periodicity.nearest_image(ri, rj, rij);
      const double dij = math::abs(rij);
      const double dij2 = dij*dij;
      const double dij3 = dij2*dij;
      const double dij5 = dij2*dij3;
      const double x = (ri - rj)[0], y = (ri - rj)[1], z = (ri - rj)[2];

      double ck[n_ah];
      // Calculate the five tensor factors
      ck[0] = (x*x - z*z)/dij2; // [1]
      ck[1] = (y*y - z*z)/dij2; // [1]
      ck[2] = (2*x*y)/dij2; // [1]
      ck[3] = (2*x*z)/dij2; // [1]
      ck[4] = (2*y*z)/dij2; // [1]

      conf_it->curr[k] = 0.0;
      for(int i=0; i<n_ah; ++i) {
        conf_it->curr[k] += ck[i] * conf_it->Tensor[i]; // [1]
      } // Loop over tensor components
      // multiply with RDC_max to get current RDC value
      conf_it->curr[k] *= RDC_max(it) / dij3; // [1/ps]

      // update average RDC values
      if(sim.param().rdc.mode == simulation::rdc_restr_av ||
         sim.param().rdc.mode == simulation::rdc_restr_av_weighted ||
         sim.param().rdc.mode == simulation::rdc_restr_biq ||
         sim.param().rdc.mode == simulation::rdc_restr_biq_weighted){
        const double exp = pow(M_E, -sim.time_step_size()/sim.param().rdc.tau); // [1]
        conf_it->av[k] *= exp; // [1/ps]
        conf_it->av[k] += (1.0 - exp) * conf_it->curr[k]; // [1/ps]
        DEBUG(15, "Delta t, tau, e^(- Delta t/tau): " << scientific << sim.time_step_size() << ", " <<  sim.param().rdc.tau << ", " << pow(M_E, -sim.time_step_size()/sim.param().rdc.tau) )
        DEBUG(10, "interaction " << k << ":  R-inst, R-av: " << conf_it->curr[k] << ", " << conf_it->av[k] )
      }

      const double force_coefficient = sim.param().rdc.K * it->weight * math::n_avogadro; // [kJ*ps^2/mol]

      if(sim.param().rdc.mode == simulation::rdc_restr_inst ||
         sim.param().rdc.mode == simulation::rdc_restr_inst_weighted){
        // rdc_energy = .5 * K * (RDC - RDC_0)^2
        // split the interaction energy between the two atoms and possibly between two energy groups
        conf.current().energies.rdc_energy[topo.atom_energy_group()[it->i]] += 0.5 * force_coefficient * flat_bottom_pot(conf_it->curr[k], it->R0, sim.param().rdc.delta); // [kJ/mol]
        conf.current().energies.rdc_energy[topo.atom_energy_group()[it->j]] += 0.5 * force_coefficient * flat_bottom_pot(conf_it->curr[k], it->R0, sim.param().rdc.delta); // [kJ/mol]
        DEBUG(15, "energy added to groups " << topo.atom_energy_group()[it->i] << " and " << topo.atom_energy_group()[it->j] )
        DEBUG(15, "group energies are " << conf.current().energies.rdc_energy[topo.atom_energy_group()[it->i]] << " and " << conf.current().energies.rdc_energy[topo.atom_energy_group()[it->j]])
        DEBUG(10, "interaction_energy (inst): " << scientific << force_coefficient * flat_bottom_pot(conf_it->curr[k], it->R0, sim.param().rdc.delta))
      }
      else if(sim.param().rdc.mode == simulation::rdc_restr_av ||
              sim.param().rdc.mode == simulation::rdc_restr_av_weighted){
        // rdc_energy = .5 * K * (<RDC> - RDC_0)^2
        conf.current().energies.rdc_energy[topo.atom_energy_group()[it->i]] += 0.5 * force_coefficient * flat_bottom_pot(conf_it->av[k], it->R0, sim.param().rdc.delta); // [kJ/mol]
        conf.current().energies.rdc_energy[topo.atom_energy_group()[it->j]] += 0.5 * force_coefficient * flat_bottom_pot(conf_it->av[k], it->R0, sim.param().rdc.delta); // [kJ/mol]
        DEBUG(15, "energy added to groups " << topo.atom_energy_group()[it->i] << " and " << topo.atom_energy_group()[it->j] )
        DEBUG(15, "group energies are " << conf.current().energies.rdc_energy[topo.atom_energy_group()[it->i]] << " and " << conf.current().energies.rdc_energy[topo.atom_energy_group()[it->j]])
        DEBUG(10, "interaction_energy (av): " << scientific << force_coefficient * flat_bottom_pot(conf_it->av[k], it->R0, sim.param().rdc.delta))
      }
      else if(sim.param().rdc.mode == simulation::rdc_restr_biq ||
              sim.param().rdc.mode == simulation::rdc_restr_biq_weighted){
        // rdc_energy = .5 * K * (RDC - RDC_0)^2 * (<RDC> - RDC_0)^2
        conf.current().energies.rdc_energy[topo.atom_energy_group()[it->i]] += 0.5 * 2.0 * force_coefficient * flat_bottom_pot(conf_it->curr[k], it->R0, sim.param().rdc.delta) * flat_bottom_pot(conf_it->av[k], it->R0, sim.param().rdc.delta); // [kJ/mol]
        conf.current().energies.rdc_energy[topo.atom_energy_group()[it->j]] += 0.5 * 2.0 * force_coefficient * flat_bottom_pot(conf_it->curr[k], it->R0, sim.param().rdc.delta) * flat_bottom_pot(conf_it->av[k], it->R0, sim.param().rdc.delta); // [kJ/mol]
        DEBUG(15, "energy added to groups " << topo.atom_energy_group()[it->i] << " and " << topo.atom_energy_group()[it->j] )
        DEBUG(15, "group energies are " << conf.current().energies.rdc_energy[topo.atom_energy_group()[it->i]] << " and " << conf.current().energies.rdc_energy[topo.atom_energy_group()[it->j]])
        DEBUG(10, "interaction_energy (biq): " << scientific << 2.0 *force_coefficient * flat_bottom_pot(conf_it->curr[k], it->R0, sim.param().rdc.delta) * flat_bottom_pot(conf_it->av[k], it->R0, sim.param().rdc.delta))
      }

      // forces on atoms i and j
      math::VArray dQ_over_dr(n_ah);
      //in the following expressions a factor of r^-2 is omitted
      dQ_over_dr[0] = math::Vec( 2*x * (1 - (x*x - z*z)/dij2), 2*y * (  - (x*x - z*z)/dij2), 2*z * (-1 - (x*x - z*z)/dij2) ); // [nm]
      dQ_over_dr[1] = math::Vec( 2*x * (  - (y*y - z*z)/dij2), 2*y * (1 - (y*y - z*z)/dij2), 2*z * (-1 - (y*y - z*z)/dij2) ); // [nm]
      dQ_over_dr[2] = math::Vec( 2*y * (1 - (2*x*x)/dij2),     2*x * (1 - (2*y*y)/dij2),     2*z * (   - (2*x*y)/dij2)     ); // [nm]
      dQ_over_dr[3] = math::Vec( 2*z * (1 - (2*x*x)/dij2),     2*y * (  - (2*x*z)/dij2),     2*x * ( 1 - (2*z*z)/dij2)     ); // [nm]
      dQ_over_dr[4] = math::Vec( 2*x * (  - (2*y*z)/dij2),     2*z * (1 - (2*y*y)/dij2),     2*y * ( 1 - (2*z*z)/dij2)     ); // [nm]

      math::Vec force_tensor_comp(0.0, 0.0, 0.0);
      for(int h=0; h<n_ah; ++h) {
        force_tensor_comp += conf_it->Tensor[h] * (-3.0 * ck[h] * rij + dQ_over_dr[h]); // [nm]
      }

      // Add the force contribution to the forces on the atoms i and j
      // and here the factor of r^-2 comes in again
      const math::Vec dD_dr = force_tensor_comp * RDC_max(it)/dij5; // [kJ/(nm*mol)]
      // return forces
      if(sim.param().rdc.mode == simulation::rdc_restr_inst ||
         sim.param().rdc.mode == simulation::rdc_restr_inst_weighted){
        const double dV_dD = force_coefficient * dflat_bottom_pot(conf_it->curr[k], it->R0, sim.param().rdc.delta); // [kJ/(nm^2*mol)]
        conf.current().force(it->i) -= dV_dD * dD_dr; // [kJ/(nm*mol)]
        conf.current().force(it->j) += dV_dD * dD_dr; // [kJ/(nm*mol)]
        DEBUG(15, "dV_dD; dD_dr: " << dV_dD << ", (" << scientific << dD_dr[0] << ", " << dD_dr[1] << ", " << dD_dr[2] << ")")
      }
      else if(sim.param().rdc.mode == simulation::rdc_restr_av ||
              sim.param().rdc.mode == simulation::rdc_restr_av_weighted){
        const double dVav_dDav = force_coefficient * dflat_bottom_pot(conf_it->av[k], it->R0, sim.param().rdc.delta); // [kJ/(nm^2*mol)]
        // in the following we set dDav_dD to either 1.0 or [1-e^(-dt/tau)].
        // The latter is correct, the former allows for using the same K as without time AV
        const double dDav_dD = sim.param().rdc.tAVfactor ? (1.0-pow(M_E, -sim.time_step_size()/sim.param().rdc.tau)) : 1.0;
        conf.current().force(it->i) -= dVav_dDav * dDav_dD * dD_dr; // [kJ/(nm*mol)]
        conf.current().force(it->j) += dVav_dDav * dDav_dD * dD_dr; // [kJ/(nm*mol)]
        DEBUG(15, "dVav_dDav; dDav_dD; dD_dr: " << dVav_dDav << ", " << dDav_dD << ", (" << scientific << dD_dr[0] << ", " << dD_dr[1] << ", " << dD_dr[2] << ")")
      }
      else if(sim.param().rdc.mode == simulation::rdc_restr_biq ||
              sim.param().rdc.mode == simulation::rdc_restr_biq_weighted){
        const double dVbq_dD = 2.0 * force_coefficient * dflat_bottom_pot(conf_it->curr[k], it->R0, sim.param().rdc.delta) * flat_bottom_pot(conf_it->av[k], it->R0, sim.param().rdc.delta); // [kJ/(nm^2*mol)]
        const double dVbq_dDav = 2.0 * force_coefficient * flat_bottom_pot(conf_it->curr[k], it->R0, sim.param().rdc.delta) * dflat_bottom_pot(conf_it->av[k], it->R0, sim.param().rdc.delta); // [kJ/(nm^2*mol)]
        // in the following we set dDav_dD to either 1.0 or [1-e^(-dt/tau)] or 0.
        // The first is simple, the second correct.
        double dDav_dD = 0;  // case 2
        if (sim.param().rdc.tAVfactor == 0) dDav_dD = 1.0;
        else if (sim.param().rdc.tAVfactor == 1) dDav_dD = 1.0-pow(M_E, -sim.time_step_size()/sim.param().rdc.tau);
        conf.current().force(it->i) -= (dVbq_dD + dVbq_dDav * dDav_dD) * dD_dr; // [kJ/(nm*mol)]
        conf.current().force(it->j) += (dVbq_dD + dVbq_dDav * dDav_dD) * dD_dr; // [kJ/(nm*mol)]
        DEBUG(15, "dVbq_dD; dVbq_dDav; dDav_dD; dD_dr: " << dVbq_dD << ", " << dVbq_dDav << ", "<< dDav_dD << ", (" << scientific << dD_dr[0] << ", " << dD_dr[1] << ", " << dD_dr[2] << ")")
      }
    } // Loop over RDCs
  } // loop over rdc groups
}


template<math::boundary_enum B, math::virial_enum V>
int _calculate_interactions_tensor(topology::Topology & topo,
                                    configuration::Configuration & conf,
                                    simulation::Simulation & sim) {

#ifdef DEBUG
  cout.precision(14);
  setw(20);
  cout.setf(ios::fixed, ios::floatfield);
#endif

  DEBUG(5, "Mode chosen: " <<  sim.param().rdc.method)

  switch(sim.param().rdc.method) {
    case simulation::rdc_em: {
      /////////////////////////////////////////////
      //                                         //
      //      E    M I N I M I S A T I O N       //
      //                                         //
      /////////////////////////////////////////////

#ifndef FINITE_DIFF
      // 1. return a_h
      _calculate_ah<B>(topo, conf, sim);
#endif // FINITE_DIFF
      // 2. write energy of rdc interaction
      // 3. write back forces on atoms due to RDC interaction
      // the following writes back interaction energies and forces for each RDC
      _calculate_forces_atoms_T<B>(topo, conf, sim);

      vector<configuration::Configuration::special_struct::rdc_struct>::iterator
          conf_it = conf.special().rdc.begin(),
          conf_to = conf.special().rdc.end();
      vector<vector<topology::rdc_restraint_struct> >::iterator
          topo_it = topo.rdc_restraints().begin();
      for( ; conf_it!=conf_to; conf_it++, topo_it++){

        DEBUG(8, "a1-5: " << conf_it->Tensor[0] << ", " << conf_it->Tensor[1] << ", " << conf_it->Tensor[2] << ", " << conf_it->Tensor[3] << ", " << conf_it->Tensor[4])

#ifdef SCRIPT
        double denominator=0.0, enumerator=0.0;
        int k=0;
        vector<topology::rdc_restraint_struct>::iterator
          it = topo_it->begin(),
          to = topo_it->end();
        for(; it!=to; ++it, ++k) {
          enumerator += pow(conf_it->curr[k] - it->R0, 2);
          denominator += pow(it->R0, 2);
        }
        double Q = sqrt(enumerator / denominator);
        cout << "SCRIPT " << conf_it->Tensor[0] << " " << conf_it->Tensor[1] << " " << conf_it->Tensor[2] << " " << conf_it->Tensor[3] << " " << conf_it->Tensor[4] << " " << Q << endl;
#endif
      }

      break;
    }

    case simulation::rdc_sd: {
      ///////////////////////////////////////////
      //                                       //
      // S T O C H A S T I C   D Y N A M I C S //
      //                                       //
      ///////////////////////////////////////////
      DEBUG(7, "doing RDC/T/SD")

#ifndef FINITE_DIFF
      vector<configuration::Configuration::special_struct::rdc_struct>::iterator
          conf_it = conf.special().rdc.begin(),
          conf_to = conf.special().rdc.end();
      for( ; conf_it!=conf_to; conf_it++){
        if(sim.steps() == 0 && abs(conf_it->Tensor[0]+2) < 1.e-10) { // init clm to sensible values in the first cycle, if user didn't input values her- or himself.
          _calculate_ah<B>(topo, conf, sim);
        }
      }
#endif // FINITE_DIFF

      const int n_ah = 5;
      //initialise the vector holding forces on tensor components
      vector<double> force_array_tensor(n_ah, 0.0);

      _calculate_forces_tensor_T<B>(topo, conf, sim, force_array_tensor);

      /////////////////////////
      //      LEAP FROG      //
      /////////////////////////

      const double gamma = sim.param().rdc.sdfric;
      const double gamma_abs = fabs(gamma);

      // Calculate SD coefficients
      const double gdt = gamma * sim.time_step_size();
      const double gdth = gdt * 0.5;

      conf_it = conf.special().rdc.begin();
      conf_to = conf.special().rdc.end();
      vector<vector<topology::rdc_restraint_struct> >::iterator
          topo_it = topo.rdc_restraints().begin();
      for( ; conf_it!=conf_to; conf_it++, topo_it++){

#ifdef DEBUG
        // Create random number generator, set seed to 0 (to make sure it's reproducible)
        static math::RandomGenerator* rng = math::RandomGenerator::create(sim.param(), "0");
#else
        // Create random number generator, choose a remotely random random-seed
        ostringstream ss;
        ss << time(NULL);
        static math::RandomGenerator* rng = math::RandomGenerator::create(sim.param(), ss.str());
#endif // DEBUG

        // Save old velocities and positions
        vector<double> old_Tensor(conf_it->Tensor);
        vector<double> old_TensorVel(conf_it->TensorVel);

        vector<double> kT_over_m(n_ah);
        for (int h=0; h<n_ah; ++h){
           kT_over_m[h] = sqrt(math::k_Boltzmann * sim.param().rdc.temp / conf_it->TensorMass[h]);
           DEBUG(15, "kT_over_m[" << h << "] = " << kT_over_m[h])
        }

        double c1 = 0.0, c2 = 0.0, c3 = 0.0, c4 = 0.0, c5 = 0.0, c6 = 0.0, c7 = 0.0, c8 = 0.0, c9 = 0.0;
        if(fabs(gdt) > 0.05){
          DEBUG(10, "doing the analytical formulas")

          const double emdth = exp(-gdth);
          const double epdth = exp(+gdth);
          const double emdt  = emdth * emdth;
          const double epdt  = epdth * epdth;
          const double omdt  = 1.0 - emdt;
          const double cdth  = gdt - 3.0 + 4.0 * emdth - emdt;
          const double ddth  = 2.0 - epdth - emdth;
          const double bpdth = gdt * (epdt - 1.0) - 4.0 * (epdth - 1.0) * (epdth - 1.0);
          const double bmdth = gdt * (1.0 - emdt) - 4.0 * (emdth - 1.0) * (emdth - 1.0);

          c1 = emdt;
          c2 = omdt / gdt;
          c3 = sqrt(fabs(omdt));
          c4 = sqrt(fabs(bpdth/cdth));
          c5 = gamma_abs * ddth/cdth;
          c6 = (epdth - emdth) / gdt;
          c7 = sqrt(fabs(cdth)) / gamma_abs;
          c8 = sqrt(fabs(bmdth/omdt)) / gamma_abs;
          c9 = -ddth/(gamma_abs * omdt);

          if(sim.steps() == 0){ // and we do all this to make sure the stochastic integral is initialized
            for (int h=0; h<n_ah; ++h){
              const double sd = kT_over_m[h] / gamma * sqrt(cdth);
              rng->stddev(sd);
              conf_it->stochastic_integral_t[h] = rng->get_gauss();
            }
          }
        } // calculate coefficients
        else {
          DEBUG(10, "doing the power series")
          //Power Series
          const double gdth2 = gdth * gdth;
          const double gdth3 = gdth2 * gdth;
          const double gdth4 = gdth2 * gdth2;
          const double gdth5 = gdth2 * gdth3;
          const double gdth6 = gdth3 * gdth3;
          const double gdth7 = gdth4 * gdth3;

          c1 = exp(-gdt);
          c2 = 1.0 - gdth + gdth2 * 2.0/3.0 - gdth3/3.0 + gdth4 * 2.0/15.0 - gdth5 * 2.0/45.0 + gdth6 * 4.0/315.0;
          c3 = sqrt(fabs(c2) * 2.0 * gdth);
          c4 = sqrt(fabs(gdth/2.0 + gdth2 * 7.0/8.0 + gdth3 * 367.0/480.0 + gdth4 * 857.0/1920.0 + gdth5 * 52813.0/268800.0 + gdth6 * 224881.0/3225600.0 + gdth7 * 1341523.0/64512000.0));
          c5 = -2.0 / sim.time_step_size() * (1.5 + gdth * 9.0/8.0 + gdth2 * 71.0/160.0 + gdth3 * 81.0/640.0 + gdth4 * 7807.0/268800.0 + gdth5 * 1971.0/358400.0 + gdth6 * 56417.0/64512000.0);
          c6 = 1.0 + gdth2/6.0 + gdth4/10.0 + gdth6/5040.0;
          c7 = sim.time_step_size() * 0.5 * sqrt(fabs(gdth * 2.0/3.0 - gdth2/2.0 + gdth3 * 7.0/30.0 -gdth4/12.0 + gdth5 * 31.0/1260.0 - gdth6/160.0 + gdth7 * 127.0/90720.0));
          c8 = sim.time_step_size() * 0.5 * sqrt(fabs(gdth/6.0 - gdth3/60.0 + gdth5 * 17.0/10080.0 - gdth7 * 31.0/181440.0));
          c9 = sim.time_step_size() * 0.5 * (0.5 + gdth/2.0 + gdth2 * 5.0/24.0 + gdth3/24.0 + gdth4/240.0 + gdth5/720.0 + gdth6 * 5.0/8064.0);

          if(sim.steps() == 0) { // and we do all this to make sure the stochastic integral is initialized
            for (int h=0; h<n_ah; ++h){
              if (gamma < math::epsilon) {
                conf_it->stochastic_integral_t[h] = 0.0;
              }
              else {
                const double emdth = exp(-gdth);
                const double emdt  = emdth * emdth;
                const double cdth  = gdt - 3.0 + 4.0 * emdth - emdt;

                const double sd = kT_over_m[h] / gamma * sqrt(cdth);
                rng->stddev(sd);
                conf_it->stochastic_integral_t[h] = rng->get_gauss();
              }
            }
          }
        } // Power series

        DEBUG(10, "c1 through c9: " << c1 << ", " << c2 << ", " << c3 << ", " << c4 << ", " << c5 << ", " << c6 << ", " << c7 << ", " << c8 << ", " << c9)

        conf_it->Ekin = 0.0;

#ifndef FINITE_DIFF
        for(int h=0; h<n_ah; ++h){
          // Vel1
          DEBUG(10, "doing stochastic dynamics velocities (1)")
          double vrand1 = 0.0, vrand2 = 0.0, vrand3 = 0.0, vrand4 = 0.0;
          //2.11.2.8 (sigma2 squared)
          const double sd1 = c3 * kT_over_m[h];
          //2.11.2.9 (rho1 squared)
          const double sd2 = c4 * kT_over_m[h];
          // (rho2 squared)
          const double sd3 = c7 * kT_over_m[h];
          // (?? squared)
          const double sd4 = c8 * kT_over_m[h];
          DEBUG(15, "sd1=" << sd1)
          DEBUG(15, "sd2=" << sd2)
          DEBUG(15, "sd3=" << sd3)
          DEBUG(15, "sd4=" << sd4)
          // Sample the V' vector from 2.11.2.20 from Gaussian with mean=0.0, sd=sd1
          // but only 1D in this case
          rng->stddev(sd1);
          vrand1 = rng->get_gauss();
          // Sample the V vector from 2.11.2.2 from Gaussian with mean=0.0, sd=sd2
          // but only 1D in this case
          rng->stddev(sd2);
          vrand2 = rng->get_gauss();
          // Sample the R' vector from 2.11.2.25 from Gaussian with mean=0.0, sd=sd3
          // but only 1D in this case
          rng->stddev(sd3);
          vrand3 = rng->get_gauss();
          // Sample the R vector from 2.11.2.26 from Gaussian with mean=0.0, sd=sd4
          // but only 1D in this case
          rng->stddev(sd4);
          vrand4 = rng->get_gauss();
          DEBUG(10, "vrand1=" << vrand1)
          DEBUG(10, "vrand2=" << vrand2)
          DEBUG(10, "vrand3=" << vrand3)
          DEBUG(10, "vrand4=" << vrand4)

          // Update the velocity
          // and if we wrote to conf.old it would be much simpler
          const double svh = conf_it->stochastic_integral_t[h] * c5 + vrand2; //save the old integral // FIXME really?
          const double cf = sim.time_step_size() * c2 / conf_it->TensorMass[h]; // FIXME what is this
          // Assign vrand1 to the stochastic integrals
          conf_it->stochastic_integral_t[h] = vrand1;
          // 2.11.2.2 (2.11.2.3 minus first term of the third line)
          conf_it->TensorVel[h] = ((conf_it->TensorVel[h] - svh) * c1)
              + (force_array_tensor[h] * cf)
              + vrand1;

          // Pos1
          DEBUG(10, "doing stochastic dynamics: positions (1)")
          DEBUG(10, "old pos(" << h << ") " << conf_it->Tensor[h])
          conf_it->Tensor[h] += conf_it->TensorVel[h] * sim.time_step_size() * c6;
          DEBUG(10, "new pos(" << h << ") " << conf_it->Tensor[h])

          // Vel2
          DEBUG(10, "doing stochastic dynamics: velocities (2)")
          // Velocity correction 2.11.2.24
          const double cinv = 1.0 / (c6 * sim.time_step_size());
          const double r = conf_it->Tensor[h] -  old_Tensor[h];
          conf_it->TensorVel[h] = r * cinv;

          // Pos2
          DEBUG(10, "doing stochastic dynamics: positions (2)")
          //2.11.2.25
          const double sxh = conf_it->stochastic_integral_t[h] * c9 + vrand4;
          conf_it->stochastic_integral_t[h] = vrand3;
          //2.11.2.26
          conf_it->Tensor[h] += (vrand3 - sxh);

          // there is no SHAKE here because the length of tensor components is not constrained (except for extreme values, and then we have other problems)

          // 'Ekin' of the tensor
          // E = 0.5 * m * v^2 = 0.25 * m * (v_{n-.5}^2 + v_{n+.5}^2)
          // m[u = g/mol], v[nm/ps]
          conf_it->Ekin += 0.25 * conf_it->TensorMass[h] * ( pow(old_TensorVel[h],2) + pow(conf_it->TensorVel[h],2) ); // [kJ/mol]

        } // Loop over tensor components (5)
#endif // FINITE_DIFF
        DEBUG(8, "a1-5: " << conf_it->Tensor[0] << ", " << conf_it->Tensor[1] << ", " << conf_it->Tensor[2] << ", " << conf_it->Tensor[3] << ", " << conf_it->Tensor[4])
        DEBUG(8, "Ekin of tensor: " << conf_it->Ekin)

        for (int h=0; h<n_ah; ++h){
          DEBUG(12, "stochastic integral (= vrand3) " << h << ": " << conf_it->stochastic_integral_t[h])
        }

#ifdef SCRIPT
        double denominator=0.0, enumerator=0.0;
        int k=0;
        vector<topology::rdc_restraint_struct>::iterator
          it = topo_it->begin(),
          to = topo_it->end();
        for(; it!=to; ++it, ++k) {
          enumerator += pow(conf_it->curr[k] - it->R0, 2);
          denominator += pow(it->R0, 2);
        }
        double Q = sqrt(enumerator / denominator);
        cout << "SCRIPT " << conf_it->Tensor[0] << " " << conf_it->Tensor[1] << " " << conf_it->Tensor[2] << " " << conf_it->Tensor[3] << " " << conf_it->Tensor[4] << " " <<  conf_it->Ekin << " " << Q << endl;
#endif
      } // loop over rdc groups

      // 1. return potential energy of interaction
      // 2. return forces on atoms
      // the following writes back interaction energies and forces for each RDC
      _calculate_forces_atoms_T<B>(topo, conf, sim);
      // 3. coordinates of tensor (done implicitly)
      // 4. velocities of tensor (done implicitly)
      // 5. return random state // FIXME check
      // 6. "kinetic energy" of the tensor (written to the struct)

      break;
    }

    case simulation::rdc_md: {
      /////////////////////////////////////////
      //                                     //
      // M O L E C U L A R   D Y N A M I C S //
      //                                     //
      /////////////////////////////////////////
      DEBUG(7, "doing RDC/T/MD")

#ifndef FINITE_DIFF
      vector<configuration::Configuration::special_struct::rdc_struct>::iterator
          conf_it = conf.special().rdc.begin(),
          conf_to = conf.special().rdc.end();
      for( ; conf_it!=conf_to; conf_it++){
        if(sim.steps() == 0 && abs(conf_it->Tensor[0]+2) < 1.e-10) { // init clm to sensible values in the first cycle, if user didn't input values her- or himself.
          _calculate_ah<B>(topo, conf, sim);
        }
      }
#endif // FINITE_DIFF

      const int n_ah = 5;
      //initialise the vector holding tensor force vectors
      vector<double> force_array_tensor(n_ah, 0.0);

      _calculate_forces_tensor_T<B>(topo, conf, sim, force_array_tensor);

      ///////////////
      // LEAP FROG //
      ///////////////

      conf_it = conf.special().rdc.begin();
      conf_to = conf.special().rdc.end();
      vector<vector<topology::rdc_restraint_struct> >::iterator
          topo_it = topo.rdc_restraints().begin();
      for( ; conf_it!=conf_to; conf_it++, topo_it++){

        conf_it->Ekin = 0.0;

#ifndef FINITE_DIFF
        for(int h=0; h<n_ah; ++h){

          // tensor components move in 1D
          double vold = conf_it->TensorVel[h]; // [nm/ps]

          // Velocity step v = v + f*dt/m // f[kJ/(nm*mol)], dt[ps], m[u]
          // 1u = M_u/N_A = g/(6*10^23) && M_u = 1g/mol && N_A = 6*10^23 1/mol
          // but slightly handwavy we can use u = g/mol
          conf_it->TensorVel[h] += force_array_tensor[h] * sim.time_step_size() / conf_it->TensorMass[h]; // [nm/ps]
          // Position step r = r + v*dt // r[nm], v[nm/ps], dt[ps]
          conf_it->Tensor[h] += conf_it->TensorVel[h] * sim.time_step_size(); // [nm]

          // there is no SHAKE here because the length of tensor components is not constrained (except for extreme values, and then we have other problems)

          // 'Ekin' of the tensor components
          // E = 0.5 * m * v^2 = 0.25 * m * (v_{n-.5}^2 + v_{n+.5}^2)
          // m[u = g/mol], v[nm/ps]
          conf_it->Ekin += 0.25 * conf_it->TensorMass[h] * (pow(vold,2) + pow(conf_it->TensorVel[h],2)); // [kJ/mol]
        } // Loop over tensor components
#endif // FINITE_DIFF
        DEBUG(8, "a1-5: " << conf_it->Tensor[0] << ", " << conf_it->Tensor[1] << ", " << conf_it->Tensor[2] << ", " << conf_it->Tensor[3] << ", " << conf_it->Tensor[4])
        DEBUG(8, "Ekin of tensor: " << conf_it->Ekin)

#ifdef SCRIPT
        double denominator=0.0, enumerator=0.0;
        int k=0;
        vector<topology::rdc_restraint_struct>::iterator
          it = topo_it->begin(),
          to = topo_it->end();
        for(; it!=to; ++it, ++k) {
          enumerator += pow(conf_it->curr[k] - it->R0, 2);
          denominator += pow(it->R0, 2);
        }
        double Q = sqrt(enumerator / denominator);
        cout << "SCRIPT " << conf_it->Tensor[0] << " " << conf_it->Tensor[1] << " " << conf_it->Tensor[2] << " " << conf_it->Tensor[3] << " " << conf_it->Tensor[4] << " " <<  conf_it->Ekin << " " << Q << endl;
#endif
      } // loop over rdc groups

      // 1. return potential energy of interaction
      // 2. return forces on atoms
      // the following writes back interaction energies and forces for each RDC
      _calculate_forces_atoms_T<B>(topo, conf, sim);
      // 3. coordinates of tensor (implicitly done)
      // 4. velocities of tensor (implicitly done)
      // 5. "kinetic energy" of the tensor (directly copied to the struct)

      break;
    }
    default: {
      // should never be the case
      assert(false);
    }
  }
  return 0;
}




//***************************************************************************************
//***************************************************************************************
//**********************                                           **********************
//**********************   S P H E R I C A L   H A R M O N I C S   **********************
//**********************                                           **********************
//***************************************************************************************
//***************************************************************************************


// Y_{2,-2} & = d_{xy} = i \sqrt{\frac{1}{2}} \left( Y_2^{- 2} - Y_2^2\right) = \frac{1}{2} \sqrt{\frac{15}{\pi}} \cdot \frac{x y}{r^2}
// Y_{2,-1} & = d_{yz} = i \sqrt{\frac{1}{2}} \left( Y_2^{- 1} + Y_2^1 \right) = \frac{1}{2} \sqrt{\frac{15}{\pi}} \cdot \frac{y z}{r^2}
// Y_{20} & = d_{z^2} = Y_2^0 = \frac{1}{4} \sqrt{\frac{5}{\pi}} \cdot \frac{- x^2 - y^2 + 2 z^2}{r^2}
// Y_{21} & = d_{xz} = \sqrt{\frac{1}{2}} \left( Y_2^{- 1} - Y_2^1 \right) = \frac{1}{2} \sqrt{\frac{15}{\pi}} \cdot \frac{z x}{r^2}
// Y_{22} & = d_{x^2-y^2} = \sqrt{\frac{1}{2}} \left( Y_2^{- 2} + Y_2^2 \right) = \frac{1}{4} \sqrt{\frac{15}{\pi}} \cdot \frac{x^2 - y^2 }{r^2}

// hardcode it, as we only ever need the five second order ones
// the following functions are normalized, i.e. \int\int fkt(m,x,y,z)^2 sin(theta) \d\theta\d\phi = 1
inline double real_harmonic_2nd_order(const int m, const double x, const double y, const double z){
  // factors of 1/r^2 are omitted, because r=1
  switch(m){
    case -2: return 0.5  * sqrt(15.0/M_PI) * x * y;
    case -1: return 0.5  * sqrt(15.0/M_PI) * y * z;
    case  0: return 0.25 * sqrt(5.0/M_PI)  * (3.0*z*z - 1.0);
    case  1: return 0.5  * sqrt(15.0/M_PI) * x * z;
    case  2: return 0.25 * sqrt(15.0/M_PI) * (x*x - y*y);
    default: return 0.0;
  }
}


// calculate forces on the clm coefficients
template<math::boundary_enum B>
void _calculate_forces_clm_SH(topology::Topology & topo,
                   configuration::Configuration & conf,
                   const simulation::Simulation & sim,
                   vector<double> &force_array_clm /* contains all zeros at this point*/) {

  DEBUG(10, "getting forces on clm components")

#ifdef DEBUG
  cout.precision(14);
  setw(20);
  cout.setf(ios::fixed, ios::floatfield);
#endif

  // create variables
  math::Periodicity<B> periodicity(conf.current().box);
  math::VArray &pos = conf.current().pos;

  const int l = 2;

  // hardcoded Lebedev Grid of order 5
  const double aux_a = 2.0/30.0 * 4.0*M_PI;
  const double aux_b = 3.0/40.0 * 4.0*M_PI;
  const double aux_c = 1.0/sqrt(3.0);
  const int grid_number_of_points = 14;
  const double grid_weight[grid_number_of_points] = {aux_a, aux_a, aux_a, aux_a, aux_a, aux_a, aux_b,  aux_b,  aux_b,  aux_b,  aux_b,  aux_b,  aux_b,  aux_b};
  const double grid_x[grid_number_of_points] =      {1.0,   -1.0,  0.0,   0.0,   0.0,   0.0,   aux_c, -aux_c,  aux_c, -aux_c,  aux_c, -aux_c,  aux_c, -aux_c};
  const double grid_y[grid_number_of_points] =      {0.0,    0.0,  1.0,  -1.0,   0.0,   0.0,   aux_c,  aux_c, -aux_c, -aux_c,  aux_c,  aux_c, -aux_c, -aux_c};
  const double grid_z[grid_number_of_points] =      {0.0,    0.0,  0.0,   0.0,   1.0,  -1.0,   aux_c,  aux_c,  aux_c,  aux_c, -aux_c, -aux_c, -aux_c, -aux_c};

  vector<configuration::Configuration::special_struct::rdc_struct>::iterator
      conf_it = conf.special().rdc.begin(),
      conf_to = conf.special().rdc.end();
  vector<vector<topology::rdc_restraint_struct> >::iterator
      topo_it = topo.rdc_restraints().begin();
  for( ; conf_it!=conf_to; conf_it++, topo_it++){

    // precalculate all required values of the weight function g
    double g[grid_number_of_points];
    for(int i=0; i<grid_number_of_points; i++){g[i]=0.0;} // init
    for(int m=-l; m<=l; m++){
      //const double c = conf.special().rdc.clm[l+m];
      const double c = conf_it->clm[l+m];
      for(int i=0; i<grid_number_of_points; i++){
        g[i] += c * real_harmonic_2nd_order(m, grid_x[i], grid_y[i], grid_z[i]); // [1]
      }
    }

    vector<topology::rdc_restraint_struct>::iterator
        it = topo_it->begin(),
        to = topo_it->end();
    int k = 0;
    // Loop over all RDC constraints
    for( ; it!=to; ++it, ++k){
      // Get positions of the two atoms involved in this RDC, calculate d, x, y, z
      math::Vec ri(0.0, 0.0, 0.0);
      math::Vec rj(0.0, 0.0, 0.0);
      math::Vec rij(0.0, 0.0, 0.0);
      ri = pos(it->i);
      rj = pos(it->j);
      periodicity.nearest_image(ri, rj, rij);
      const double dij = math::abs(rij); // [nm]
      const double dij2 = dij*dij; // [nm^2]
      const double dij3 = dij2*dij; // [nm^3]

      // weighting is potentially switched off by setting weight to one during initialisation
      // K[kJ*ps^2], weight[1], na[1/mol]
      const double force_coefficient = sim.param().rdc.K * it->weight * math::n_avogadro; // [kJ*ps^2/mol]

      conf_it->curr[k] = 0.0;
      // get current RDC value
      for(int i=0; i<grid_number_of_points; i++){
        const math::Vec rh(grid_x[i], grid_y[i], grid_z[i]); // [nm]
        const double RDC_sampled = RDC_max(it)/dij3 * (1.5 * pow(math::dot(rij, rh) / dij, 2) - 0.5); // [1/ps]
        // save current RDC values
        conf_it->curr[k] += RDC_sampled * g[i] * grid_weight[i]; // [1/ps]
      } // loop over grid points
      // locally calculate average RDC values, use for forces and don't export
      // the average RDC is not updated because that is done when forces on atoms are computed
      double local_average = 0.0;
      if(sim.param().rdc.mode == simulation::rdc_restr_av ||
         sim.param().rdc.mode == simulation::rdc_restr_av_weighted ||
         sim.param().rdc.mode == simulation::rdc_restr_biq ||
         sim.param().rdc.mode == simulation::rdc_restr_biq_weighted){
        const double exp = pow(M_E, -sim.time_step_size()/sim.param().rdc.tau); // [1]
        local_average = exp * conf_it->av[k] + (1.0 - exp) * conf_it->curr[k]; // [1/ps]
        DEBUG(15, "Delta t, tau, e^(- Delta t/tau): " << scientific << sim.time_step_size() << ", " <<  sim.param().rdc.tau << ", " << pow(M_E, -sim.time_step_size()/sim.param().rdc.tau) )
        DEBUG(10, "interaction " << k << ":  R-inst, R-av: " << conf_it->curr[k] << ", " << local_average )
      }

      for(int m=-l; m<=l; m++){
        double dD_dc = 0.0;
        for(int i=0; i<grid_number_of_points; i++){
          const math::Vec rh(grid_x[i], grid_y[i], grid_z[i]); // [nm]
          const double RDC_sampled = RDC_max(it)/dij3 * (1.5 * pow(math::dot(rij, rh) / dij, 2) - 0.5); // [1/ps]
          dD_dc += RDC_sampled * real_harmonic_2nd_order(m, grid_x[i], grid_y[i], grid_z[i]) * grid_weight[i]; // [1/nm*ps] // non-physical nm for consistency
        }
        // forces on clm depending on AV mode
        if(sim.param().rdc.mode == simulation::rdc_restr_inst ||
           sim.param().rdc.mode == simulation::rdc_restr_inst_weighted){
          const double dV_dD = force_coefficient * dflat_bottom_pot(conf_it->curr[k], it->R0, sim.param().rdc.delta); // [kJ/(nm^2*mol)]
          force_array_clm[l+m] -= dV_dD * dD_dc; // [kJ/(nm*mol)]
          DEBUG(15, "dV_dD; dD_dc: " << dV_dD << ", " << scientific << dD_dc)
        }
        else if(sim.param().rdc.mode == simulation::rdc_restr_av ||
                sim.param().rdc.mode == simulation::rdc_restr_av_weighted){
          const double dVav_dDav = force_coefficient * dflat_bottom_pot(local_average, it->R0, sim.param().rdc.delta); // [kJ/(nm^2*mol)]
          // in the following we set dDav_dD to either 1.0 or [1-e^(-dt/tau)].
          // The latter is correct, the former allows for using the same K as without time AV
          const double dDav_dD = sim.param().rdc.tAVfactor ? (1.0-pow(M_E, -sim.time_step_size()/sim.param().rdc.tau)) : 1.0;
          force_array_clm[l+m] -= dVav_dDav * dDav_dD * dD_dc; // [kJ/(nm*mol)]
          DEBUG(15, "dVav_dDav; dDav_dD; dD_dc: " << dVav_dDav << ", " << dDav_dD << ", " << scientific << dD_dc)
        }
        else if(sim.param().rdc.mode == simulation::rdc_restr_biq ||
                sim.param().rdc.mode == simulation::rdc_restr_biq_weighted){
          const double dVbq_dD = 2.0 * force_coefficient * dflat_bottom_pot(conf_it->curr[k], it->R0, sim.param().rdc.delta) * flat_bottom_pot(local_average, it->R0, sim.param().rdc.delta); // [kJ/(nm^2*mol)]
          const double dVbq_dDav = 2.0 * force_coefficient * flat_bottom_pot(conf_it->curr[k], it->R0, sim.param().rdc.delta) * dflat_bottom_pot(local_average, it->R0, sim.param().rdc.delta); // [kJ/(nm^2*mol)]
          // in the following we set dDav_dD to either 1.0 or [1-e^(-dt/tau)] or 0.
          // The first is simple, the second correct.
          double dDav_dD = 0;  // case 2
          if (sim.param().rdc.tAVfactor == 0) dDav_dD = 1.0;
          else if (sim.param().rdc.tAVfactor == 1) dDav_dD = 1.0-pow(M_E, -sim.time_step_size()/sim.param().rdc.tau);
          force_array_clm[l+m] -= (dVbq_dD + dVbq_dDav * dDav_dD) * dD_dc; // [kJ/(nm*mol)]
          DEBUG(15, "dVbq_dD; dVbq_dDav; dDav_dD; dD_dc: " << dVbq_dD << ", " << dVbq_dDav << ", "<< dDav_dD << ", " << scientific << dD_dc)
        }
      }
    } // loop over RDCs in one group
    DEBUG(12, "dV_dclm[1-5]: " << scientific << force_array_clm[0] << ", " << force_array_clm[1] << ", " << force_array_clm[2] << ", " << force_array_clm[3] << ", " << force_array_clm[4])


  } // loop over rdc groups
  //write_clm_forces(force_array_clm, "trfsh");
}



// calculate RDC interaction and forces on atoms in the SH representation
template<math::boundary_enum B>
void _calculate_forces_atoms_SH(topology::Topology & topo,
                   configuration::Configuration & conf,
                   const simulation::Simulation & sim) {

#ifdef DEBUG
  cout.precision(14);
  setw(20);
  cout.setf(ios::fixed, ios::floatfield);
#endif

  // create variables
  math::Periodicity<B> periodicity(conf.current().box);
  math::VArray &pos = conf.current().pos;
  math::Vec ri(0.0, 0.0, 0.0);
  math::Vec rj(0.0, 0.0, 0.0);
  math::Vec rij(0.0, 0.0, 0.0);

  const int l = 2;

  // hardcoded LebedevGrid of order 5
  const double aux_a = 2.0/30.0 * 4.0*M_PI;
  const double aux_b = 3.0/40.0 * 4.0*M_PI;
  const double aux_c = 1.0/sqrt(3.0);
  const int grid_number_of_points = 14;
  const double grid_weight[grid_number_of_points] = {aux_a, aux_a, aux_a, aux_a, aux_a, aux_a, aux_b,  aux_b,  aux_b,  aux_b,  aux_b,  aux_b,  aux_b,  aux_b};
  const double grid_x[grid_number_of_points] =      {1.0,   -1.0,  0.0,   0.0,   0.0,   0.0,   aux_c, -aux_c,  aux_c, -aux_c,  aux_c, -aux_c,  aux_c, -aux_c};
  const double grid_y[grid_number_of_points] =      {0.0,    0.0,  1.0,  -1.0,   0.0,   0.0,   aux_c,  aux_c, -aux_c, -aux_c,  aux_c,  aux_c, -aux_c, -aux_c};
  const double grid_z[grid_number_of_points] =      {0.0,    0.0,  0.0,   0.0,   1.0,  -1.0,   aux_c,  aux_c,  aux_c,  aux_c, -aux_c, -aux_c, -aux_c, -aux_c};

  vector<configuration::Configuration::special_struct::rdc_struct>::iterator
      conf_it = conf.special().rdc.begin(),
      conf_to = conf.special().rdc.end();
  vector<vector<topology::rdc_restraint_struct> >::iterator
      topo_it = topo.rdc_restraints().begin();
  for( ; conf_it!=conf_to; conf_it++, topo_it++){

    // precalculate all required values of the weight function g
    double g[grid_number_of_points];
    for(int i=0; i<grid_number_of_points; i++){g[i]=0.0;} // init
    for(int m=-l; m<=l; m++){
      const double c = conf_it->clm[l+m]; // [nm] // pseudo nm ...
      for(int i=0; i<grid_number_of_points; i++){
        g[i] += c * real_harmonic_2nd_order(m, grid_x[i], grid_y[i], grid_z[i]); // [nm] // pseudo nm ...
      }
    }
    for(int i=0; i<grid_number_of_points; i++){DEBUG(15, "g[" << i << "]: " << scientific << g[i])}

    // Iterate over all RDC constraints
    vector<topology::rdc_restraint_struct>::iterator
        it = topo_it->begin(),
        to = topo_it->end();
    int k=0;
    for(; it!=to; ++it, ++k){
      // Get positions of the two atoms involved in this RDC, calculate dij, x, y, z
      ri = pos(it->i);
      rj = pos(it->j);
      periodicity.nearest_image(ri, rj, rij);
      const double dij = math::abs(rij);
      const double dij3 = pow(dij,3);
      const double dij5 = pow(dij,5);

      conf_it->curr[k] = 0.0;
      // get current RDC values
      for(int i=0; i<grid_number_of_points; i++){
        const math::Vec rh(grid_x[i], grid_y[i], grid_z[i]); // [nm]
        const double RDC_sampled = RDC_max(it)/dij3 * (1.5 * pow(math::dot(rij, rh) / dij, 2) - 0.5); // [1/ps]
        // save current RDC value
        conf_it->curr[k] += RDC_sampled * g[i] * grid_weight[i]; // [1/ps]
      } // loop over grid points
      // update average RDC values
      if(sim.param().rdc.mode == simulation::rdc_restr_av ||
         sim.param().rdc.mode == simulation::rdc_restr_av_weighted ||
         sim.param().rdc.mode == simulation::rdc_restr_biq ||
         sim.param().rdc.mode == simulation::rdc_restr_biq_weighted){
        const double exp = pow(M_E, -sim.time_step_size()/sim.param().rdc.tau); // [1]
        conf_it->av[k] *= exp; // [1/ps]
        conf_it->av[k] += (1.0 - exp) * conf_it->curr[k]; // [1/ps]
        DEBUG(15, "Delta t, tau, e^(- Delta t/tau): " << scientific << sim.time_step_size() << ", " <<  sim.param().rdc.tau << ", " << pow(M_E, -sim.time_step_size()/sim.param().rdc.tau) )
        DEBUG(10, "interaction " << k << ":  R-inst, R-av: " << conf_it->curr[k] << ", " << conf_it->av[k] )
      }

      // weighting is potentially switched off by setting weight to one during initialisation
      const double force_coefficient = sim.param().rdc.K * it->weight * math::n_avogadro; // [kJ*ps^2/mol]

      if(sim.param().rdc.mode == simulation::rdc_restr_inst ||
         sim.param().rdc.mode == simulation::rdc_restr_inst_weighted){
        // rdc_energy = .5 * K * (RDC - RDC_0)^2
        // split the interaction energy between the two atoms and possibly between two energy groups
        conf.current().energies.rdc_energy[topo.atom_energy_group()[it->i]] += 0.5 * force_coefficient * flat_bottom_pot(conf_it->curr[k], it->R0, sim.param().rdc.delta); // [kJ/mol]
        conf.current().energies.rdc_energy[topo.atom_energy_group()[it->j]] += 0.5 * force_coefficient * flat_bottom_pot(conf_it->curr[k], it->R0, sim.param().rdc.delta); // [kJ/mol]
        DEBUG(15, "energy added to groups " << topo.atom_energy_group()[it->i] << " and " << topo.atom_energy_group()[it->j] )
        DEBUG(15, "group energies are " << conf.current().energies.rdc_energy[topo.atom_energy_group()[it->i]] << " and " << conf.current().energies.rdc_energy[topo.atom_energy_group()[it->j]])
        DEBUG(10, "interaction_energy (inst): " << scientific << force_coefficient * flat_bottom_pot(conf_it->curr[k], it->R0, sim.param().rdc.delta))
      }
      else if(sim.param().rdc.mode == simulation::rdc_restr_av ||
              sim.param().rdc.mode == simulation::rdc_restr_av_weighted){
        // rdc_energy = .5 * K * (<RDC> - RDC_0)^2
        conf.current().energies.rdc_energy[topo.atom_energy_group()[it->i]] += 0.5 * force_coefficient * flat_bottom_pot(conf_it->av[k], it->R0, sim.param().rdc.delta); // [kJ/mol]
        conf.current().energies.rdc_energy[topo.atom_energy_group()[it->j]] += 0.5 * force_coefficient * flat_bottom_pot(conf_it->av[k], it->R0, sim.param().rdc.delta); // [kJ/mol]
        DEBUG(15, "energy added to groups " << topo.atom_energy_group()[it->i] << " and " << topo.atom_energy_group()[it->j] )
        DEBUG(15, "group energies are " << conf.current().energies.rdc_energy[topo.atom_energy_group()[it->i]] << " and " << conf.current().energies.rdc_energy[topo.atom_energy_group()[it->j]])
        DEBUG(10, "interaction_energy (av): " << scientific << force_coefficient * flat_bottom_pot(conf_it->av[k], it->R0, sim.param().rdc.delta))
      }
      else if(sim.param().rdc.mode == simulation::rdc_restr_biq ||
              sim.param().rdc.mode == simulation::rdc_restr_biq_weighted){
        // rdc_energy = .5 * K * (RDC - RDC_0)^2 * (<RDC> - RDC_0)^2
        conf.current().energies.rdc_energy[topo.atom_energy_group()[it->i]] += 0.5 * force_coefficient * flat_bottom_pot(conf_it->curr[k], it->R0, sim.param().rdc.delta) * flat_bottom_pot(conf_it->av[k], it->R0, sim.param().rdc.delta); // [kJ/mol]
        conf.current().energies.rdc_energy[topo.atom_energy_group()[it->j]] += 0.5 * force_coefficient * flat_bottom_pot(conf_it->curr[k], it->R0, sim.param().rdc.delta) * flat_bottom_pot(conf_it->av[k], it->R0, sim.param().rdc.delta); // [kJ/mol]
        DEBUG(15, "energy added to groups " << topo.atom_energy_group()[it->i] << " and " << topo.atom_energy_group()[it->j] )
        DEBUG(15, "group energies are " << conf.current().energies.rdc_energy[topo.atom_energy_group()[it->i]] << " and " << conf.current().energies.rdc_energy[topo.atom_energy_group()[it->j]])
        DEBUG(10, "interaction_energy (biq): " << scientific << 2.0 * force_coefficient * flat_bottom_pot(conf_it->curr[k], it->R0, sim.param().rdc.delta) * flat_bottom_pot(conf_it->av[k], it->R0, sim.param().rdc.delta))
      }

      // forces on atoms
      // dD_dr is used for all three cases
      math::Vec dD_dr(0.0, 0.0, 0.0);
      for(int i=0; i<grid_number_of_points; i++){
        const math::Vec rh(grid_x[i], grid_y[i], grid_z[i]); // [nm]
        assert(abs(math::abs(rh) - 1.0) < 1.e-5);
        const double dot_product = math::dot(rij, rh); // [nm^2]
        // the second half of the following expression is zero if rij is of constant length
#ifdef FINITE_DIFF
        dD_dr += g[i] * grid_weight[i] * (dot_product * rh - (2.5 * pow(dot_product/dij, 2) - 0.5) * rij); // [nm]
#else
        dD_dr += g[i] * grid_weight[i] * (dot_product * rh); // [nm]
#endif
      }
      dD_dr *= RDC_max(it) * 3.0/dij5;

      // return forces
      if(sim.param().rdc.mode == simulation::rdc_restr_inst ||
         sim.param().rdc.mode == simulation::rdc_restr_inst_weighted){
        const double dV_dD = force_coefficient * dflat_bottom_pot(conf_it->curr[k], it->R0, sim.param().rdc.delta); // [kJ/(nm^2*mol)]
        conf.current().force(it->i) -= dV_dD * dD_dr; // [kJ/(nm*mol)]
        conf.current().force(it->j) += dV_dD * dD_dr; // [kJ/(nm*mol)]
        DEBUG(15, "dV_dD; dD_dr: " << dV_dD << ", (" << scientific << dD_dr[0] << ", " << dD_dr[1] << ", " << dD_dr[2] << ")")
      }
      else if(sim.param().rdc.mode == simulation::rdc_restr_av ||
              sim.param().rdc.mode == simulation::rdc_restr_av_weighted){
        const double dVav_dDav = force_coefficient * dflat_bottom_pot(conf_it->av[k], it->R0, sim.param().rdc.delta); // [kJ/(nm^2*mol)]
        // in the following we set dDav_dD to either 1.0 or [1-e^(-dt/tau)].
        // The latter is correct, the former allows for using the same K as without time AV
        const double dDav_dD = sim.param().rdc.tAVfactor ? (1.0-pow(M_E, -sim.time_step_size()/sim.param().rdc.tau)) : 1.0;
        conf.current().force(it->i) -= dVav_dDav * dDav_dD * dD_dr; // [kJ/(nm*mol)]
        conf.current().force(it->j) += dVav_dDav * dDav_dD * dD_dr; // [kJ/(nm*mol)]
        DEBUG(15, "dVav_dDav; dDav_dD; dD_dr: " << dVav_dDav << ", " << dDav_dD << ", (" << scientific << dD_dr[0] << ", " << dD_dr[1] << ", " << dD_dr[2] << ")")
      }
      else if(sim.param().rdc.mode == simulation::rdc_restr_biq ||
              sim.param().rdc.mode == simulation::rdc_restr_biq_weighted){
        const double dVbq_dD = 2.0 * force_coefficient * dflat_bottom_pot(conf_it->curr[k], it->R0, sim.param().rdc.delta) * flat_bottom_pot(conf_it->av[k], it->R0, sim.param().rdc.delta); // [kJ/(nm^2*mol)]
        const double dVbq_dDav = 2.0 * force_coefficient * flat_bottom_pot(conf_it->curr[k], it->R0, sim.param().rdc.delta) * dflat_bottom_pot(conf_it->av[k], it->R0, sim.param().rdc.delta); // [kJ/(nm^2*mol)]
        // in the following we set dDav_dD to either 1.0 or [1-e^(-dt/tau)] or 0.
        // The first is simple, the second correct.
        double dDav_dD = 0;  // case 2
        if (sim.param().rdc.tAVfactor == 0) dDav_dD = 1.0;
        else if (sim.param().rdc.tAVfactor == 1) dDav_dD = 1.0-pow(M_E, -sim.time_step_size()/sim.param().rdc.tau);
        conf.current().force(it->i) -= (dVbq_dD + dVbq_dDav * dDav_dD) * dD_dr; // [kJ/(nm*mol)]
        conf.current().force(it->j) += (dVbq_dD + dVbq_dDav * dDav_dD) * dD_dr; // [kJ/(nm*mol)]
        DEBUG(15, "dVbq_dD; dVbq_dDav; dDav_dD; dD_dr: " << dVbq_dD << ", " << dVbq_dDav << ", "<< dDav_dD << ", (" << scientific << dD_dr[0] << ", " << dD_dr[1] << ", " << dD_dr[2] << ")")
      }
    } // Loop over RDCs in group
  } // loop over rdc groups
}


// calculate optimal clm
template<math::boundary_enum B>
void _calculate_clm(topology::Topology & topo,
                   configuration::Configuration & conf,
                   const simulation::Simulation & sim) {

#ifdef DEBUG
  cout.precision(14);
  setw(20);
  cout.setf(ios::fixed, ios::floatfield);
#endif

  // create variables
  math::Periodicity<B> periodicity(conf.current().box);
  math::VArray &pos = conf.current().pos;
  math::Vec ri(0.0, 0.0, 0.0);
  math::Vec rj(0.0, 0.0, 0.0);
  math::Vec rij(0.0, 0.0, 0.0);

  const int l = 2;
  const int n_clm = 5;

  // hardcoded LebedevGrid of order 5
  const double aux_a = 2.0/30.0 * 4.0*M_PI;
  const double aux_b = 3.0/40.0 * 4.0*M_PI;
  const double aux_c = 1.0/sqrt(3.0);
  const int grid_number_of_points = 14;
  const double grid_weight[grid_number_of_points] = {aux_a, aux_a, aux_a, aux_a, aux_a, aux_a, aux_b,  aux_b,  aux_b,  aux_b,  aux_b,  aux_b,  aux_b,  aux_b};
  const double grid_x[grid_number_of_points] =      {1.0,   -1.0,  0.0,   0.0,   0.0,   0.0,   aux_c, -aux_c,  aux_c, -aux_c,  aux_c, -aux_c,  aux_c, -aux_c};
  const double grid_y[grid_number_of_points] =      {0.0,    0.0,  1.0,  -1.0,   0.0,   0.0,   aux_c,  aux_c, -aux_c, -aux_c,  aux_c,  aux_c, -aux_c, -aux_c};
  const double grid_z[grid_number_of_points] =      {0.0,    0.0,  0.0,   0.0,   1.0,  -1.0,   aux_c,  aux_c,  aux_c,  aux_c, -aux_c, -aux_c, -aux_c, -aux_c};

  vector<configuration::Configuration::special_struct::rdc_struct>::iterator
      conf_it = conf.special().rdc.begin(),
      conf_to = conf.special().rdc.end();
  vector<vector<topology::rdc_restraint_struct> >::iterator
      topo_it = topo.rdc_restraints().begin();
  for( ; conf_it!=conf_to; conf_it++, topo_it++){

    const int n_rdc = topo_it->size();

    // alpha_klm = 1/{4\pi} \int^{2\pi}\int^{\pi} D(\theta_k) Y_m^l(\theta\varphi) sin(\theta) d\theta d\varphi [in units of 1/ps]
    double alpha_klm[n_rdc][n_clm];
    for(int k=0; k<n_rdc; ++k){for(int m=-l; m<=l; ++m){alpha_klm[k][l+m]=0.0;}} // init
    double matrix_array[n_clm*n_clm];
    for(int i=0; i<n_clm*n_clm; ++i){matrix_array[i]=0.0;} // init
    double result_array[n_clm];
    for(int i=0; i<n_clm; ++i){result_array[i]=0.0;} // init

    vector<topology::rdc_restraint_struct>::iterator
        it = topo_it->begin(),
        to = topo_it->end();
    int k=0; // counter of rdcs
    // Loop over all RDC constraints
    for(; it!=to; ++it, ++k){
      // Get positions of the two atoms involved in this RDC, calculate d
      ri = pos(it->i);
      rj = pos(it->j);
      periodicity.nearest_image(ri, rj, rij);
      const double dij = math::abs(rij); // [nm]
      const double dij3 =pow(dij,3); // [nm^3]

      for(int m=-l; m<=l; ++m){
        double num_integral = 0.0;
        for(int i=0; i<grid_number_of_points; i++){
          const math::Vec rh(grid_x[i], grid_y[i], grid_z[i]); // [nm]
  //        DEBUG(15, "rh: " << rh[0] << "," << rh[1] << "," << rh[2] << ":  " << sqrt(rh[0]*rh[0] + rh[1]*rh[1] + rh[2]*rh[2]))
  //        DEBUG(15, "rij: " << rij[0] << "," << rij[1] << "," << rij[2] << ":  " << sqrt(rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2]))
          const double RDC_sampled = 1.5 * pow(math::dot(rij, rh) / dij, 2) - 0.5; // [1]
          num_integral += grid_weight[i] * real_harmonic_2nd_order(m, grid_x[i], grid_y[i], grid_z[i]) * RDC_sampled; // [1/nm] // pseudo-nm ...
  //        DEBUG(15, num_integral << " (updated)")
        }
        alpha_klm[k][l+m] = num_integral; // [1/nm] // pseudo-nm ...
        DEBUG(15, "alpha_klm (" << scientific << k << "," << l << "," << m << "): " << num_integral)
      }

      for(int m_prime=-l; m_prime<=l; ++m_prime){
        for(int m=-l; m<=l; ++m){
          // weighting is potentially switched off by setting weight to one during initialisation
          matrix_array[(l+m_prime)*n_clm + l+m] += it->weight * alpha_klm[k][l+m] * alpha_klm[k][l+m_prime]; // [1/nm^2]
        }
        // weighting is potentially switched off by setting weight to one during initialisation
        result_array[l+m_prime] += it->weight * (it->R0 * dij3 / RDC_max(it)) * alpha_klm[k][l+m_prime]; // [1/nm]
      }
    }

    for(int i=0; i<n_clm; ++i){
      for(int j=0; j<n_clm; ++j){
        DEBUG(15, "matrix[" << i*n_clm + j << "]: " << scientific << matrix_array[i*n_clm + j])
      }
    }
    for(int j=0; j<n_clm; ++j){
      DEBUG(15, "results vector[" << j << "]: " << scientific << result_array[j])
    }

    // solve the SLE using gsl
    gsl_matrix_view m = gsl_matrix_view_array (matrix_array, n_clm, n_clm);
    gsl_vector_view b = gsl_vector_view_array (result_array, n_clm);
    gsl_vector *x = gsl_vector_alloc (n_clm);
    int signum = 0; // not used
    gsl_permutation * p = gsl_permutation_alloc (n_clm);

    gsl_linalg_LU_decomp (&m.matrix, p, &signum);
    //switch off and save error handler
    gsl_error_handler_t* old_handler = gsl_set_error_handler_off();
    int status = gsl_linalg_LU_solve (&m.matrix, p, &b.vector, x);
    gsl_set_error_handler(old_handler);

    if(status == GSL_SUCCESS){ // copy back c values
      for(int lm=0; lm<n_clm; ++lm){
        conf_it->clm[lm] = gsl_vector_get(x, lm); // [nm]
        DEBUG(15, "c[" << lm << "]: " << conf_it->clm[lm])
      }
    }else if(status == GSL_EDOM){
      // ignore singular matrices which could very rarely occur
      // do not update c_lm
      io::messages.add("Singular matrix encountered in gsl_linalg_LU_solve.  This may be harmless.", "rdc_restraint_interaction", io::message::warning);
    }else{
      io::messages.add("Error in gsl_linalg_LU_solve", "rdc_restraint_interaction", io::message::error);
    }

    gsl_permutation_free (p);
    gsl_vector_free (x);

  } //loop over rdc groups
}


template<math::boundary_enum B, math::virial_enum V>
int _calculate_interactions_sh(topology::Topology & topo,
                                configuration::Configuration & conf,
                                simulation::Simulation & sim) {

#ifdef DEBUG
  cout.precision(14);
  setw(20);
  cout.setf(ios::fixed, ios::floatfield);
#endif

  DEBUG(5, "Mode chosen: " <<  sim.param().rdc.method)

  switch(sim.param().rdc.method) {
    case simulation::rdc_em: {
      ///////////////////////////////////////////
      //                                       //
      // E N E R G Y   M I N I M I S A T I O N //
      //                                       //
      ///////////////////////////////////////////

#ifndef FINITE_DIFF
      // 1. return clm
      _calculate_clm<B>(topo, conf, sim);
#endif // FINITE_DIFF
      // 2. write energy of rdc interaction
      // 3. forces due to rdc interaction
      // forces on all atoms are not zero (but their sum is)
      // the following writes back interaction energies and forces for each RDC
      _calculate_forces_atoms_SH<B>(topo, conf, sim);

      vector<configuration::Configuration::special_struct::rdc_struct>::iterator
          conf_it = conf.special().rdc.begin(),
          conf_to = conf.special().rdc.end();
      vector<vector<topology::rdc_restraint_struct> >::iterator
          topo_it = topo.rdc_restraints().begin();
      for( ; conf_it!=conf_to; conf_it++, topo_it++){

        DEBUG(8, "clm: " << conf_it->clm[0] << ", " << conf_it->clm[1] << ", " << conf_it->clm[2] << ", " << conf_it->clm[3] << ", " << conf_it->clm[4])

#ifdef SCRIPT
        double denominator=0.0, enumerator=0.0;
        int k=0;
        vector<topology::rdc_restraint_struct>::iterator
          it = topo_it->begin(),
          to = topo_it->end();
        for(; it!=to; ++it, ++k) {
          enumerator += pow(conf_it->curr[k] - it->R0, 2);
          denominator += pow(it->R0, 2);
        }
        double Q = sqrt(enumerator / denominator);
        cout << "SCRIPT " << conf_it->clm[0] << " " << conf_it->clm[1] << " " << conf_it->clm[2] << " " << conf_it->clm[3] << " " << conf_it->clm[4] << " " << Q << endl;
#endif
      } // loop over rdc groups

      break;
    }


    case simulation::rdc_sd: {
      ///////////////////////////////////////////////
      //                                           //
      //   S T O C H A S T I C   D Y N A M I C S   //
      //                                           //
      ///////////////////////////////////////////////
      DEBUG(7, "doing RDC/SH/SD")

#ifndef FINITE_DIFF
      vector<configuration::Configuration::special_struct::rdc_struct>::iterator
          conf_it = conf.special().rdc.begin(),
          conf_to = conf.special().rdc.end();
      for( ; conf_it!=conf_to; conf_it++){
        if(sim.steps() == 0 && abs(conf_it->clm[0]+2) < 1.e-10) { // init clm to sensible values in the first cycle, if user didn't input values her- or himself.
          _calculate_clm<B>(topo, conf, sim);
        }
      }
#endif // FINITE_DIFF

      const int n_clm = 5;
      vector<double> force_array_clm(n_clm, 0.0);

      _calculate_forces_clm_SH<B>(topo, conf, sim, force_array_clm);

      /////////////////////////
      //      LEAP FROG      //
      /////////////////////////

      const double gamma = sim.param().rdc.sdfric;
      const double gamma_abs = fabs(gamma);

      // Calculate SD coefficients
      const double gdt = gamma * sim.time_step_size();
      const double gdth = gdt * 0.5;

      conf_it = conf.special().rdc.begin();
      conf_to = conf.special().rdc.end();
      vector<vector<topology::rdc_restraint_struct> >::iterator
          topo_it = topo.rdc_restraints().begin();
      for( ; conf_it!=conf_to; conf_it++, topo_it++){

#ifdef DEBUG
        // Create random number generator, set seed to 0 (to make sure it's reproducible)
        static math::RandomGenerator* rng = math::RandomGenerator::create(sim.param(), "0");
#else
        // Create random number generator, choose a remotely random random-seed
        ostringstream ss;
        ss << time(NULL);
        static math::RandomGenerator* rng = math::RandomGenerator::create(sim.param(), ss.str());
#endif // DEBUG

        // Save old velocities and positions
        vector<double> old_clm(conf_it->clm);
        vector<double> old_clmVel(conf_it->clmVel);

        vector<double> kT_over_m(n_clm);
        for (int lm=0; lm<n_clm; ++lm){
           kT_over_m[lm] = sqrt(math::k_Boltzmann * sim.param().rdc.temp / conf_it->clmMass[lm]);
           DEBUG(15, "kT_over_m[" << lm << "] = " << kT_over_m[lm])
           DEBUG(18, "kB, temp, mass: " << math::k_Boltzmann << ", " << sim.param().rdc.temp << ", " << conf_it->clmMass[lm])
        }

        double c1 = 0.0, c2 = 0.0, c3 = 0.0, c4 = 0.0, c5 = 0.0, c6 = 0.0, c7 = 0.0, c8 = 0.0, c9 = 0.0;
        if(fabs(gdt) > 0.05){
          DEBUG(10, "doing the analytical formulas")

          const double emdth = exp(-gdth);
          const double epdth = exp(+gdth);
          const double emdt  = emdth * emdth;
          const double epdt  = epdth * epdth;
          const double omdt  = 1.0 - emdt;
          const double cdth  = gdt - 3.0 + 4.0 * emdth - emdt;
          const double ddth  = 2.0 - epdth - emdth;
          const double bpdth = gdt * (epdt - 1.0) - 4.0 * (epdth - 1.0) * (epdth - 1.0);
          const double bmdth = gdt * (1.0 - emdt) - 4.0 * (emdth - 1.0) * (emdth - 1.0);

          c1 = emdt;
          c2 = omdt / gdt;
          c3 = sqrt(fabs(omdt));
          c4 = sqrt(fabs(bpdth/cdth));
          c5 = gamma_abs * ddth/cdth;
          c6 = (epdth - emdth) / gdt;
          c7 = sqrt(fabs(cdth)) / gamma_abs;
          c8 = sqrt(fabs(bmdth/omdt)) / gamma_abs;
          c9 = -ddth/(gamma_abs * omdt);

          if(sim.steps() == 0){ // and we do all this to make sure the stochastic integral is initialized
            for (int lm=0; lm<n_clm; ++lm){
              const double sd = kT_over_m[lm] / gamma * sqrt(cdth);
              rng->stddev(sd);
              conf_it->stochastic_integral_sh[lm] = rng->get_gauss();
            }
          }
        } // calculate coefficients
        else {
          DEBUG(10, "doing the power series")
          //Power Series
          const double gdth2 = gdth * gdth;
          const double gdth3 = gdth2 * gdth;
          const double gdth4 = gdth2 * gdth2;
          const double gdth5 = gdth2 * gdth3;
          const double gdth6 = gdth3 * gdth3;
          const double gdth7 = gdth4 * gdth3;

          c1 = exp(-gdt);
          c2 = 1.0 - gdth + gdth2 * 2.0/3.0 - gdth3/3.0 + gdth4 * 2.0/15.0 - gdth5 * 2.0/45.0 + gdth6 * 4.0/315.0;
          c3 = sqrt(fabs(c2) * 2.0 * gdth);
          c4 = sqrt(fabs(gdth/2.0 + gdth2 * 7.0/8.0 + gdth3 * 367.0/480.0 + gdth4 * 857.0/1920.0 + gdth5 * 52813.0/268800.0 + gdth6 * 224881.0/3225600.0 + gdth7 * 1341523.0/64512000.0));
          c5 = -2.0 / sim.time_step_size() * (1.5 + gdth * 9.0/8.0 + gdth2 * 71.0/160.0 + gdth3 * 81.0/640.0 + gdth4 * 7807.0/268800.0 + gdth5 * 1971.0/358400.0 + gdth6 * 56417.0/64512000.0);
          c6 = 1.0 + gdth2/6.0 + gdth4/10.0 + gdth6/5040.0;
          c7 = sim.time_step_size() * 0.5 * sqrt(fabs(gdth * 2.0/3.0 - gdth2/2.0 + gdth3 * 7.0/30.0 -gdth4/12.0 + gdth5 * 31.0/1260.0 - gdth6/160.0 + gdth7 * 127.0/90720.0));
          c8 = sim.time_step_size() * 0.5 * sqrt(fabs(gdth/6.0 - gdth3/60.0 + gdth5 * 17.0/10080.0 - gdth7 * 31.0/181440.0));
          c9 = sim.time_step_size() * 0.5 * (0.5 + gdth/2.0 + gdth2 * 5.0/24.0 + gdth3/24.0 + gdth4/240.0 + gdth5/720.0 + gdth6 * 5.0/8064.0);

          if(sim.steps() == 0) { // and we do all this to make sure the stochastic integral is initialized
            for (int lm=0; lm<n_clm; ++lm){
              if (gamma < math::epsilon) {
                conf_it->stochastic_integral_sh[lm] = 0.0;
              }
              else {
                const double emdth = exp(-gdth);
                const double emdt  = emdth * emdth;
                const double cdth  = gdt - 3.0 + 4.0 * emdth - emdt;

                const double sd = kT_over_m[lm] / gamma * sqrt(cdth);
                rng->stddev(sd);
                conf_it->stochastic_integral_sh[lm] = rng->get_gauss();
              }
            }
          }
        } // Power series

        DEBUG(10, "c1 through c9: " << c1 << ", " << c2 << ", " << c3 << ", " << c4 << ", " << c5 << ", " << c6 << ", " << c7 << ", " << c8 << ", " << c9)

        conf_it->Ekin = 0.0;

#ifndef FINITE_DIFF
        for(int lm=0; lm<n_clm; ++lm){
          // Vel1
          DEBUG(10, "doing stochastic dynamics velocities (1)")
          double vrand1 = 0.0, vrand2 = 0.0, vrand3 = 0.0, vrand4 = 0.0;
          //2.11.2.8 (sigma2 squared)
          const double sd1 = c3 * kT_over_m[lm];
          //2.11.2.9 (rho1 squared)
          const double sd2 = c4 * kT_over_m[lm];
          // (rho2 squared)
          const double sd3 = c7 * kT_over_m[lm];
          // (?? squared)
          const double sd4 = c8 * kT_over_m[lm];
          DEBUG(15, "sd1=" << sd1)
          DEBUG(15, "sd2=" << sd2)
          DEBUG(15, "sd3=" << sd3)
          DEBUG(15, "sd4=" << sd4)
          // Sample the V' vector from 2.11.2.20 from Gaussian with mean=0.0, sd=sd1
          // but only 1D in this case
          rng->stddev(sd1);
          vrand1 = rng->get_gauss();
          // Sample the V vector from 2.11.2.2 from Gaussian with mean=0.0, sd=sd2
          // but only 1D in this case
          rng->stddev(sd2);
          vrand2 = rng->get_gauss();
          // Sample the R' vector from 2.11.2.25 from Gaussian with mean=0.0, sd=sd3
          // but only 1D in this case
          rng->stddev(sd3);
          vrand3 = rng->get_gauss();
          // Sample the R vector from 2.11.2.26 from Gaussian with mean=0.0, sd=sd4
          // but only 1D in this case
          rng->stddev(sd4);
          vrand4 = rng->get_gauss();
          DEBUG(15, "vrand1=" << vrand1)
          DEBUG(15, "vrand2=" << vrand2)
          DEBUG(15, "vrand3=" << vrand3)
          DEBUG(15, "vrand4=" << vrand4)

          // Update the velocity
          // and if we wrote to conf.old it would be much simpler
          const double svh = conf_it->stochastic_integral_sh[lm] * c5 + vrand2; // save the old integral // FIXME really?
          const double cf = sim.time_step_size() * c2 / conf_it->clmMass[lm]; // FIXME what is this
          // Assign vrand1 to the stochastic integrals
          conf_it->stochastic_integral_sh[lm] = vrand1;
          // 2.11.2.2 (2.11.2.3 minus first term of the third line)
          conf_it->clmVel[lm] = ((conf_it->clmVel[lm] - svh) * c1)
              + (force_array_clm[lm] * cf)
              + vrand1;

          // Pos1
          DEBUG(10, "doing stochastic dynamics: positions (1)")
          DEBUG(10, "old pos(" << lm << ") " << conf_it->clm[lm])
          conf_it->clm[lm] += conf_it->clmVel[lm] * sim.time_step_size() * c6;
          DEBUG(10, "new pos(" << lm << ") " << conf_it->clm[lm])

          // Vel2
          DEBUG(10, "doing stochastic dynamics: velocities (2)")
          // Velocity correction 2.11.2.24
          const double cinv = 1.0 / (c6 * sim.time_step_size());
          const double r = conf_it->clm[lm] -  old_clm[lm];
          conf_it->clmVel[lm] = r * cinv;

          // Pos2
          DEBUG(10, "doing stochastic dynamics: positions (2)")
          //2.11.2.25
          const double sxh = conf_it->stochastic_integral_sh[lm] * c9 + vrand4;
          conf_it->stochastic_integral_sh[lm] = vrand3;
          //2.11.2.26
          conf_it->clm[lm] += (vrand3 - sxh);

          //FIXME constraints should go here

          // 'Ekin' of the clm
          // E = 0.5 * m * v^2 = 0.25 * m * (v_{n-.5}^2 + v_{n+.5}^2)
          // m[u = g/mol], v[nm/ps]
          conf_it->Ekin += 0.25 * conf_it->clmMass[lm] * ( pow(old_clmVel[lm],2) + pow(conf_it->clmVel[lm],2) ); // [kJ/mol]

        } // Loop over clm
#endif // FINITE_DIFF

        if(pow(conf_it->clm[0],2) + pow(conf_it->clm[1],2) + pow(conf_it->clm[2],2) + pow(conf_it->clm[3],2) + pow(conf_it->clm[4],2) > 1.0){
          DEBUG(10, "The norm of the weight function g is exceeding 1.  Please change T, m or K to avoid this.")
        }

        DEBUG(8, "clm: " << conf_it->clm[0] << ", " << conf_it->clm[1] << ", " << conf_it->clm[2] << ", " << conf_it->clm[3] << ", " << conf_it->clm[4])
        DEBUG(8, "Ekin of clm: " << conf_it->Ekin)

        for (int lm=0; lm<n_clm; ++lm){
          DEBUG(12, "stochastic integral (= vrand3) " << lm << ": " << conf_it->stochastic_integral_sh[lm])
        }

#ifdef SCRIPT
        double denominator=0.0, enumerator=0.0;
        int k=0;
        vector<topology::rdc_restraint_struct>::iterator
          it = topo_it->begin(),
          to = topo_it->end();
        for(; it!=to; ++it, ++k) {
          enumerator += pow(conf_it->curr[k] - it->R0, 2);
          denominator += pow(it->R0, 2);
        }
        double Q = sqrt(enumerator / denominator);
        cout << "SCRIPT " << conf_it->clm[0] << " " << conf_it->clm[1] << " " << conf_it->clm[2] << " " << conf_it->clm[3] << " " << conf_it->clm[4] << " " <<  conf_it->Ekin << " " << Q << endl;
#endif
      } // loop over rdc groups

      // 1. return potential energy of interaction
      // 2. return forces on atoms
      // the following writes back interaction energies and forces for each RDC
      _calculate_forces_atoms_SH<B>(topo, conf, sim);
      // 3. values of clm (done implicitly)
      // 4. velocities of clm (done implicitly)
      // 5. return random state // FIXME check
      // 6. "kinetic energy" of the clm (copied to the struct)

      break;
    }


    case simulation::rdc_md: {
      /////////////////////////////////////////
      //                                     //
      // M O L E C U L A R   D Y N A M I C S //
      //                                     //
      /////////////////////////////////////////
      DEBUG(7, "doing RDC/SH/MD")

#ifndef FINITE_DIFF
      vector<configuration::Configuration::special_struct::rdc_struct>::iterator
          conf_it = conf.special().rdc.begin(),
          conf_to = conf.special().rdc.end();
      for( ; conf_it!=conf_to; conf_it++){
        if(sim.steps() == 0 && abs(conf_it->clm[0]+2) < 1.e-10) { // init clm to sensible values in the first cycle, if user didn't input values her- or himself.
          _calculate_clm<B>(topo, conf, sim);
        }
      }
#endif // FINITE_DIFF

      const int n_clm = 5;
      vector<double> force_array_clm(n_clm, 0.0);

      _calculate_forces_clm_SH<B>(topo, conf, sim, force_array_clm);

      /////////////////////////
      //      LEAP FROG      //
      /////////////////////////

      conf_it = conf.special().rdc.begin();
      conf_to = conf.special().rdc.end();
      vector<vector<topology::rdc_restraint_struct> >::iterator
          topo_it = topo.rdc_restraints().begin();
      for( ; conf_it!=conf_to; conf_it++, topo_it++){

        conf_it->Ekin = 0.0;

#ifndef FINITE_DIFF
        for(int lm=0; lm<n_clm; ++lm){
          // clm move in 1D
          double vold = conf_it->clmVel[lm]; // [nm/ps]

          // Velocity step v = v + f*dt/m // f[kJ/(nm*mol)], dt[ps], m[u]
          // 1u = M_u/N_A = g/(6*10^23) && M_u = 1g/mol && N_A = 6*10^23 1/mol
          // but slightly handwavy we can use u = g/mol
          conf_it->clmVel[lm] += force_array_clm[lm] * sim.time_step_size() / conf_it->clmMass[lm]; // [nm/ps]
          DEBUG(12, "force, dt, m: " << force_array_clm[lm] << ", " << sim.time_step_size() << ", " << conf_it->clmMass[lm] << " :: " << force_array_clm[lm] * sim.time_step_size() / conf_it->clmMass[lm])
          // Position step r = r + v*dt // r[nm], v[nm/ps], dt[ps]
          conf_it->clm[lm] += conf_it->clmVel[lm] * sim.time_step_size(); // [nm]
          DEBUG(12, "clm, vel[" << lm << "]: " << conf_it->clm[lm] << ", " << conf_it->clmVel[lm])

          // there is no SHAKE here because the values of clm are not constrained (except for extreme values, and then we have other problems)

          // 'Ekin' of the clm
          // E = 0.5 * m * v^2 = 0.25 * m * (v_{n-.5}^2 + v_{n+.5}^2)
          // m[u = g/mol], v[nm/ps]
          conf_it->Ekin += 0.25 * conf_it->clmMass[lm] * (pow(vold,2) + pow(conf_it->clmVel[lm],2)); // [kJ/mol]
        } // loop over clm
#endif // FINITE_DIFF

        if(pow(conf_it->clm[0],2) + pow(conf_it->clm[1],2) + pow(conf_it->clm[2],2) + pow(conf_it->clm[3],2) + pow(conf_it->clm[4],2) > 1.0){
          DEBUG(10, "The norm of the weight function g is exceeding 1.  Please change T, m or K to avoid this.")
        }

        DEBUG(8, "clm: " << conf_it->clm[0] << ", " << conf_it->clm[1] << ", " << conf_it->clm[2] << ", " << conf_it->clm[3] << ", " << conf_it->clm[4])
        DEBUG(8, "Ekin of clm: " << conf_it->Ekin)

#ifdef SCRIPT
        double denominator=0.0, enumerator=0.0;
        int k=0;
        vector<topology::rdc_restraint_struct>::iterator
          it = topo_it->begin(),
          to = topo_it->end();
        for(; it!=to; ++it, ++k) {
          cout << conf_it->curr[k] << " " << it->R0 << endl;
          enumerator += pow(conf_it->curr[k] - it->R0, 2);
          denominator += pow(it->R0, 2);
        }
        double Q = sqrt(enumerator / denominator);
        cout << "SCRIPT " << conf_it->clm[0] << " " << conf_it->clm[1] << " " << conf_it->clm[2] << " " << conf_it->clm[3] << " " << conf_it->clm[4] << " " <<  conf_it->Ekin << " " << Q << endl;
#endif
      }

      // 1. return potential energy of interaction
      // 2. return forces on atoms
      // the following writes back interaction energies and forces for each RDC
      _calculate_forces_atoms_SH<B>(topo, conf, sim);
      // 3. return Ekin (directly written to the struct)
      // 4. return coordinates of vectors (done implicitly)
      // 5. return velocities (done implicitly)

      break;
    }
    default:
      // should never be the case
      assert(false);
  }
  return 0;
}


/**
 * calculate_interactions
 */
int RDC_Restraint_Interaction::calculate_interactions
        (topology::Topology & topo,
        configuration::Configuration & conf,
        simulation::Simulation & sim) {
  m_timer.start();

  switch(sim.param().rdc.type) {
    case simulation::rdc_mf:
      DEBUG(5, "RDC calculating interactions, magnetic field vector representation")
      SPLIT_VIRIAL_BOUNDARY(_calculate_interactions_mfield, topo, conf, sim)
      break;
    case simulation::rdc_t:
      DEBUG(5, "RDC calculating interactions, alignment tensor representation")
      SPLIT_VIRIAL_BOUNDARY(_calculate_interactions_tensor, topo, conf, sim)
      break;
    case simulation::rdc_sh:
      DEBUG(5, "RDC calculating interactions, spherical harmonics representation")
      SPLIT_VIRIAL_BOUNDARY(_calculate_interactions_sh, topo, conf, sim)
      break;
    default:
      break;
  }

  m_timer.stop();
  return 0;
}

/**
 * init
 */
int RDC_Restraint_Interaction::init(
                topology::Topology &topo,
                configuration::Configuration &conf,
                simulation::Simulation &sim,
                ostream &os,
                bool quiet) {
    if (!quiet) {
    os << "RDC RESTRAINT INTERACTION\n";
    switch (sim.param().rdc.mode) {
      case simulation::rdc_restr_off:
        os << "\trestraining off\n";
        break;
      case simulation::rdc_restr_inst:
        os << "\tinstantaneous restraining\n";
        break;
      case simulation::rdc_restr_inst_weighted:
        os << "\tinstantaneous restraining, weighted\n";
        break;
      case simulation::rdc_restr_av:
        os << "\ttime averaged restraining\n";
        break;
      case simulation::rdc_restr_av_weighted:
        os << "\ttime averaged restraining, weighted\n";
        break;
      case simulation::rdc_restr_biq:
        os << "\tbiquadratic restraining\n";
        break;
      case simulation::rdc_restr_biq_weighted:
        os << "\tbiquadratic restraining, weighted\n";
        break;
    }

//    if (sim.param().rdc.read_av){
//      os << "\treading RDC averages from file\n";
//    }
    os << "END\n";
  }

  return 0;
}

}//namespace interaction

//#endif // FINITE_DIFF

//#endif // SCRIPT
