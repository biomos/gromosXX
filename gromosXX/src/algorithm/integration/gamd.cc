/**
 * @file gamd.cc
 * contains the implementation
 * for the Gaussian accelerated md class
 */

#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"

#include "gamd.h"

#undef MODULE
#undef SUBMODULE
#define MODULE algorithm
#define SUBMODULE integration

/**
 * GAMD step.
 */
int algorithm::GAMD
::apply(topology::Topology & topo,
        configuration::Configuration & conf,
        simulation::Simulation &sim)
 {
  m_timer.start();
  DEBUG(5, "GAMD: algorithm init");
  configuration::Energy ener = conf.current().energies;
  simulation::Parameter::gamd_struct params = sim.param().gamd;
  std::vector<int> used_chargegroups;
  unsigned int num_groups = unsigned(ener.bond_energy.size());
  unsigned int num_atoms = topo.num_atoms();
  // Calculate the energies for each acceleration group
  // loop over all the acceleration groups
  /*
  DEBUG(5, "GAMD: Calculate energie totals");
  for (unsigned int atom=0; atom < num_atoms; atom++){
    unsigned int chargegroup = topo.atom_energy_group()[atom];
    unsigned int gamdgroup = topo.gamd_accel_group(atom);
    // check if the charge group has already been used
    std::vector<int>::iterator it = std::find(used_chargegroups.begin(), used_chargegroups.end(), chargegroup);
    if (it == used_chargegroups.end()){
        DEBUG(7, "GAMD: Calculating energies of charge group " << chargegroup);
        // add the energies to the gamd energies for that accel group
        // dihedral term
        ener.gamd_dihedral_total[gamdgroup] += ener.dihedral_energy[chargegroup];
        ener.gamd_dihedral_total[gamdgroup] += ener.improper_energy[chargegroup];
        ener.gamd_dihedral_total[gamdgroup] += ener.crossdihedral_energy[chargegroup];
        // potential total bonded
        ener.gamd_potential_total[gamdgroup] += ener.dihedral_energy[chargegroup];
        ener.gamd_potential_total[gamdgroup] += ener.improper_energy[chargegroup];
        ener.gamd_potential_total[gamdgroup] += ener.crossdihedral_energy[chargegroup];
        ener.gamd_potential_total[gamdgroup] += ener.bond_energy[chargegroup];
        ener.gamd_potential_total[gamdgroup] += ener.angle_energy[chargegroup];
        // potential total non-bonded
        for(unsigned int chargegroup2=0; chargegroup2<num_groups; chargegroup2++){
            if (chargegroup2 != chargegroup){
                ener.gamd_potential_total[gamdgroup] +=  ener.lj_energy[chargegroup][chargegroup2];
                ener.gamd_potential_total[gamdgroup] +=  ener.crf_energy[chargegroup][chargegroup2];
                }
        }
        //save the used charge group
        used_chargegroups.push_back(chargegroup);              
    } // endif        
  } // loop over atoms 
  */
  DEBUG(5, "GAMD: Interactions calculated now calculate acceleration");

  // total energies have been computed now calculate acceleration
  if (sim.steps() > params.equilibration){
      for (unsigned int gg = 1; gg < params.agroups; gg++){
          switch (params.search){
            case simulation::cmd_search:
                switch (sim.param().gamd.form)
                {
                    case simulation::dih_boost:
                        calc_gamd_std_mean(ener.gamd_dihedral_total[gg], sim.steps(), &params.VmaxD[gg], &params.VminD[gg],
                                           &params.VmeanD[gg], &params.M2D[gg], &params.sigmaVD[gg]);  
                        break;

                    case simulation::tot_boost:
                        calc_gamd_std_mean(ener.gamd_potential_total[gg], sim.steps(), &params.VmaxT[gg], &params.VminT[gg],
                                           &params.VmeanT[gg], &params.M2T[gg], &params.sigmaVT[gg]);  
                        break;

                    case simulation::dual_boost:
                        calc_gamd_std_mean(ener.gamd_dihedral_total[gg], sim.steps(), &params.VmaxD[gg], &params.VminD[gg],
                                           &params.VmeanD[gg], &params.M2D[gg], &params.sigmaVD[gg]);
                        calc_gamd_std_mean(ener.gamd_potential_total[gg] - ener.gamd_dihedral_total[gg], sim.steps(), &params.VmaxT[gg], &params.VminT[gg],
                                           &params.VmeanT[gg], &params.M2T[gg], &params.sigmaVT[gg]);                                           
                    break;          
                default:
                    io::messages.add("Unknown functional form of gaussian accelerated md boosting potential. Should be 1 (dual acceleration), 2 (total potential energy acceleration), or 3 (dihedral acceleration)",
                                     "Forcefield", io::message::critical);
                } 
            break; // end case simulation::cmd_search

            case simulation::gamd_search:
                switch (sim.param().gamd.form){
                    case simulation::dih_boost:
                        calc_gamd_std_mean(ener.gamd_dihedral_total[gg], sim.steps(), &params.VmaxD[gg], &params.VminD[gg],
                                           &params.VmeanD[gg], &params.M2D[gg], &params.sigmaVD[gg]);

                        if (calc_gamd_E_K(params.thresh, params.dihstd, params.VmaxD[gg], params.VminD[gg], params.VmeanD[gg],
                            params.sigmaVD[gg], &params.k0D[gg], &params.kD[gg], &params.ED[gg])){
                                io::messages.add("gamd: k0 < 0 or k0 > 1, switched to lower bound threshold ","Forcefield", io::message::warning);
                                params.thresh = simulation::lower_bound;
                                }
                        break;

                    case simulation::tot_boost:
                        calc_gamd_std_mean(ener.gamd_potential_total[gg], sim.steps(), &params.VmaxT[gg], &params.VminT[gg],
                                           &params.VmeanT[gg], &params.M2T[gg], &params.sigmaVT[gg]);

                        if (calc_gamd_E_K(params.thresh, params.totstd, params.VmaxT[gg], params.VminT[gg], params.VmeanT[gg],
                            params.sigmaVT[gg], &params.k0T[gg], &params.kT[gg], &params.ET[gg])){
                                io::messages.add("gamd: k0 < 0 or k0 > 1, switched to lower bound threshold ","Forcefield", io::message::warning);
                                params.thresh = simulation::lower_bound;
                                }
                        break;

                    case simulation::dual_boost:
                        calc_gamd_std_mean(ener.gamd_dihedral_total[gg], sim.steps(), &params.VmaxD[gg], &params.VminD[gg],
                                           &params.VmeanD[gg], &params.M2D[gg], &params.sigmaVD[gg]);

                        calc_gamd_std_mean(ener.gamd_potential_total[gg] - ener.gamd_dihedral_total[gg], sim.steps(), &params.VmaxT[gg], &params.VminT[gg],
                                           &params.VmeanT[gg], &params.M2T[gg], &params.sigmaVT[gg]);

                        if (calc_gamd_E_K(params.thresh, params.dihstd, params.VmaxD[gg], params.VminD[gg], params.VmeanD[gg],
                            params.sigmaVD[gg], &params.k0D[gg], &params.kD[gg], &params.ED[gg])){
                                io::messages.add("gamd: k0 < 0 or k0 > 1, switched to lower bound threshold ","Forcefield", io::message::warning);
                                params.thresh = simulation::lower_bound;
                                }

                        if (calc_gamd_E_K(params.thresh, params.totstd, params.VmaxT[gg], params.VminT[gg], params.VmeanT[gg],
                            params.sigmaVT[gg], &params.k0T[gg], &params.kT[gg], &params.ET[gg])){
                                io::messages.add("gamd: k0 < 0 or k0 > 1, switched to lower bound threshold ","Forcefield", io::message::warning);
                                params.thresh = simulation::lower_bound;
                                }
                        break;
                    default:
                        io::messages.add("Unknown functional form of gaussian accelerated md boosting potential. Should be 1 (dual acceleration), 2 (total potential energy acceleration), or 3 (dihedral acceleration)",
                                         "Forcefield", io::message::critical);
                    } 
            break; // end case simulation::cmd_search
            case simulation::no_search:
                /* do nothing */
            break; // end case no_search
            default:
                io::messages.add("Unknown functional form of gaussian accelerated md search. Should be 0 (gamd with fixed parameters), 1 (cmd search), or 2 (gamd search)",
                                 "Forcefield", io::message::critical);;
        } // end switch search form
     } // for loop over gamd groups
  } // endif 

  // statistics updated now accelerate
  DEBUG(5, "GAMD: Accelerating");
  double prefactor;
  if (sim.param().gamd.search != simulation::cmd_search){
      for (unsigned int accelgroup = 1; accelgroup < params.agroups;  accelgroup++){

          switch (sim.param().gamd.form)
          {
              case simulation::dih_boost:
              {
                  double VE = ener.gamd_dihedral_total[accelgroup] - params.ED[accelgroup];
                  //if V < E apply boost
                  if (VE < 0){
                      prefactor = (params.kD[accelgroup] * VE) + 1;
                      ener.gamd_DV[accelgroup] = prefactor * VE/2;
                      // loop over atoms
                       for (unsigned int atom=0; atom < num_atoms; atom++){
                              conf.current().force(atom) += conf.special().gamd.dihe_force[accelgroup](atom) * (prefactor - 1);
                              //int chargegroup = topo.atom_energy_group()[atom];
                              // to virial
                              for (int a = 0; a < 3; ++a) {
                                for (int b = 0; b < 3; ++b) {
                                    conf.current().virial_tensor(b, a) += conf.special().gamd.virial_tensor_dihe[accelgroup](b, a) * (prefactor -1);
                                }
                              } // end virial **/
                      } // end loop over atoms
                  } // endif
                  break;
              }

          case simulation::tot_boost:
          {
                  double VE = ener.gamd_potential_total[accelgroup] - params.ET[accelgroup];
                  //if V < E apply boost
                  if (VE < 0){
                      prefactor = (params.kT[accelgroup] * VE) + 1;
                      ener.gamd_DV[accelgroup] = prefactor * VE/2; // 
                      for (unsigned int accelgroup2 = 0; accelgroup2 < params.agroups;  accelgroup2++){
                        // choose the correct interaction factor to scale forces between acceleration groups
                        double interaction_factor = 1.0;
                        calc_interaction_factor(accelgroup, accelgroup2, &interaction_factor);
                        //loop over the atoms
                        for (unsigned int atom=0; atom < num_atoms; atom++){

                                conf.current().force(atom) += conf.special().gamd.total_force[accelgroup][accelgroup2](atom) * (prefactor -1);
                                conf.current().force(atom) += conf.special().gamd.total_force[accelgroup2][accelgroup](atom) * (prefactor -1);

                              // to virial 
                              for (int a = 0; a < 3; ++a) {
                                for (int b = 0; b < 3; ++b) {
                                    conf.current().virial_tensor(b, a) += conf.special().gamd.virial_tensor[accelgroup][accelgroup2](b, a) * (prefactor -1);
                                    conf.current().virial_tensor(b, a) += conf.special().gamd.virial_tensor[accelgroup2][accelgroup](b, a) * (prefactor -1); 
                                }
                              } // end virial **/
                        } // end loop over atoms
                      }
                  } // endif
                  break;
            }

          case simulation::dual_boost:
          {
                  // first total potential
                  double VET = ener.gamd_potential_total[accelgroup] - params.ET[accelgroup] - ener.gamd_dihedral_total[accelgroup];
                  //if V < E apply boost
                  if (VET < 0){
                      prefactor = (params.kT[accelgroup] * VET) + 1;
                      ener.gamd_DV[accelgroup] = prefactor * VET/2; // 
                      for (unsigned int accelgroup2 = 0; accelgroup2 < params.agroups;  accelgroup2++){
                        // choose the correct interaction factor to scale forces between acceleration groups
                        double interaction_factor = 1.0;
                        calc_interaction_factor(accelgroup, accelgroup2, &interaction_factor);
                        //loop over the atoms
                        for (unsigned int atom=0; atom < num_atoms; atom++){

                                conf.current().force(atom) += conf.special().gamd.total_force[accelgroup][accelgroup2](atom) * (prefactor -1);
                                conf.current().force(atom) += conf.special().gamd.total_force[accelgroup2][accelgroup](atom) * (prefactor -1);

                              // to virial 
                              for (int a = 0; a < 3; ++a) {
                                for (int b = 0; b < 3; ++b) {
                                    conf.current().virial_tensor(b, a) += conf.special().gamd.virial_tensor[accelgroup][accelgroup2](b, a) * (prefactor -1);
                                    conf.current().virial_tensor(b, a) += conf.special().gamd.virial_tensor[accelgroup2][accelgroup](b, a) * (prefactor -1); 
                                }
                              } // end virial **/
                        }
                      } // end loop over atoms
                  } // endif
                  // dihedral term
                  double VED = ener.gamd_dihedral_total[accelgroup] - params.ED[accelgroup];
                  //if V < E apply boost
                  if (VED < 0){
                      prefactor = (params.kD[accelgroup] * VED) + 1;
                      ener.gamd_DV[accelgroup] = prefactor * VED/2;
                      // loop over atoms
                       for (unsigned int atom=0; atom < num_atoms; atom++){
                              conf.current().force(atom) += conf.special().gamd.dihe_force[accelgroup](atom) * (prefactor - 1);
                              //int chargegroup = topo.atom_energy_group()[atom];
                              // to virial
                              for (int a = 0; a < 3; ++a) {
                                for (int b = 0; b < 3; ++b) {
                                    conf.current().virial_tensor(b, a) += conf.special().gamd.virial_tensor_dihe[accelgroup](b, a) * (prefactor -1);
                                }
                              } // end virial **/
                        } // end loop over atoms
                    } // endif
                  break; 
            }         
          default:
            io::messages.add("Unknown functional form of gaussian accelerated md boosting potential. Should be 1 (dual acceleration), 2 (total potential energy acceleration), or 3 (dihedral acceleration)",
                             "Forcefield", io::message::critical);
          } // end switch
      } // end loop over acceleration groups
  } // end if
  return 0;
};
// add init in which the initial acceleration values are calculated if needed

// helper gamd methods
void algorithm::GAMD::calc_gamd_std_mean(double V, int step, double *Vmax, double *Vmin, double *Vmean, double *M2, double *sigmaV){
        double delta;
        if(V > *Vmax){
                *Vmax = V;
        }
        else if(V < *Vmin){
                *Vmin = V;
        }

        delta = V - *Vmean;
        *Vmean += delta /(double)step;
        *M2 += delta * (V - (*Vmean));

        *sigmaV = sqrt(*M2/step);
};

int algorithm::GAMD::calc_gamd_E_K(simulation::gamd_thresh_enum Eform, double sigma0, double Vmax, double Vmin, double Vmean, double sigmaV,
                   double *k0, double *k, double *E)
        {
                if (Eform == simulation::upper_bound){
                        *k0 = (1 - sigma0/sigmaV) * (Vmax-Vmin) / (Vmean-Vmin);
                        // k0 should be between 0 and 1 if not switch to lower bound energy threshold
                        if(*k0 > 1.0 || *k0 <= 0.0){
                                *k0 = (sigma0 / sigmaV) * (Vmax-Vmin) / (Vmax-Vmean);
                                *E = Vmax;
                                *k = *k0 / (Vmax - Vmin);
                                return 1;
                        }
                        else{
                                *E = Vmin + (Vmax-Vmin)/(*k0);
                                *k = *k0 / (Vmax - Vmin);
                                return 0;
                        }

                }
                else if(Eform == simulation::lower_bound){
                                *k0 = (sigma0 / sigmaV) * (Vmax-Vmin) / (Vmax-Vmean);
                                *E = Vmax;
                                *k = *k0 / (Vmax - Vmin);
                                return 0;  

                }
        };


 void algorithm::GAMD::calc_interaction_factor(int accelerationgroup, int accelerationgroup2, double *interaction_factor){
     // different ways of treating overlaping regions should be added here
     if (accelerationgroup2 > 0){
         *interaction_factor = 0.5;
     }
 };
