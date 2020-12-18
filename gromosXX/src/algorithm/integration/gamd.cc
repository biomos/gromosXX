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
 * GAMD init
 **/
int algorithm::GAMD
::init(topology::Topology &topo, 
		     configuration::Configuration &conf,
		     simulation::Simulation &sim,
		     std::ostream &os,
		     bool quiet)
 {
     // Calculate initial acceleration
     if (sim.param().gamd.search == simulation::gamd_search && sim.param().gamd.ntisearch == 0){
      for (unsigned int gg = 1; gg < sim.param().gamd.agroups; gg++){
        if (sim.param().gamd.kD[gg] == 0 && sim.param().gamd.kT[gg] == 0){
            switch (sim.param().gamd.form){
                case simulation::dih_boost:
                    if (calc_gamd_E_K(sim.param().gamd.thresh, sim.param().gamd.dihstd, sim.param().gamd.VmaxD[gg], sim.param().gamd.VminD[gg], sim.param().gamd.VmeanD[gg],
                        sim.param().gamd.sigmaVD[gg], &sim.param().gamd.k0D[gg], &sim.param().gamd.kD[gg], &sim.param().gamd.ED[gg])){
                            io::messages.add("gamd: k0 < 0 or k0 > 1, switched to lower bound threshold ","Forcefield", io::message::warning);
                            sim.param().gamd.thresh = simulation::lower_bound;
                            }
                    break;

                case simulation::tot_boost:
                    if (calc_gamd_E_K(sim.param().gamd.thresh, sim.param().gamd.totstd, sim.param().gamd.VmaxT[gg], sim.param().gamd.VminT[gg], sim.param().gamd.VmeanT[gg],
                        sim.param().gamd.sigmaVT[gg], &sim.param().gamd.k0T[gg], &sim.param().gamd.kT[gg], &sim.param().gamd.ET[gg])){
                            io::messages.add("gamd: k0 < 0 or k0 > 1, switched to lower bound threshold ","Forcefield", io::message::warning);
                            sim.param().gamd.thresh = simulation::lower_bound;
                            }
                    break;

                case simulation::dual_boost:
                    if (calc_gamd_E_K(sim.param().gamd.thresh, sim.param().gamd.dihstd, sim.param().gamd.VmaxD[gg], sim.param().gamd.VminD[gg], sim.param().gamd.VmeanD[gg],
                        sim.param().gamd.sigmaVD[gg], &sim.param().gamd.k0D[gg], &sim.param().gamd.kD[gg], &sim.param().gamd.ED[gg])){
                            io::messages.add("gamd: k0 < 0 or k0 > 1, switched to lower bound threshold ","Forcefield", io::message::warning);
                            sim.param().gamd.thresh = simulation::lower_bound;
                            }

                    if (calc_gamd_E_K(sim.param().gamd.thresh, sim.param().gamd.totstd, sim.param().gamd.VmaxT[gg], sim.param().gamd.VminT[gg], sim.param().gamd.VmeanT[gg],
                        sim.param().gamd.sigmaVT[gg], &sim.param().gamd.k0T[gg], &sim.param().gamd.kT[gg], &sim.param().gamd.ET[gg])){
                            io::messages.add("gamd: k0 < 0 or k0 > 1, switched to lower bound threshold ","Forcefield", io::message::warning);
                            sim.param().gamd.thresh = simulation::lower_bound;
                            }
                    break;
                default:
                    io::messages.add("Unknown functional form of gaussian accelerated md boosting potential. Should be 1 (dual acceleration), 2 (total potential energy acceleration), or 3 (dihedral acceleration)",
                                        "Forcefield", io::message::critical);
        }
        } // end if
      } // loop over acceleration groups
    } // end if
    else if (sim.param().gamd.search == simulation::gamd_search && sim.param().gamd.ntisearch == 1){
        io::messages.add("GAMD search should be started from an cMD search with ntisearch = 0. With ntisearch = 1 the statistics will be estimated only from the GAMD search run",
         "Forcefield", io::message::warning);
    }
    return 0;
 }



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
  std::vector<int> used_chargegroups;
  unsigned int num_groups = unsigned(ener.bond_energy.size());
  unsigned int num_atoms = topo.num_atoms();
  DEBUG(5, "GAMD: Interactions calculated now calculate acceleration");

  // total energies have been computed now calculate acceleration
  if (sim.steps() > sim.param().gamd.equilibration){
      for (unsigned int gg = 1; gg < sim.param().gamd.agroups; gg++){
          switch (sim.param().gamd.search){
            case simulation::cmd_search:
                switch (sim.param().gamd.form)
                {
                    case simulation::dih_boost:
                        calc_gamd_std_mean(ener.gamd_dihedral_total[gg], sim.steps(), &sim.param().gamd.VmaxD[gg], &sim.param().gamd.VminD[gg],
                                           &sim.param().gamd.VmeanD[gg], &sim.param().gamd.M2D[gg], &sim.param().gamd.sigmaVD[gg]);  
                        break;

                    case simulation::tot_boost:
                        calc_gamd_std_mean(ener.gamd_potential_total[gg], sim.steps(), &sim.param().gamd.VmaxT[gg], &sim.param().gamd.VminT[gg],
                                           &sim.param().gamd.VmeanT[gg], &sim.param().gamd.M2T[gg], &sim.param().gamd.sigmaVT[gg]);  
                        break;

                    case simulation::dual_boost:
                        calc_gamd_std_mean(ener.gamd_dihedral_total[gg], sim.steps(), &sim.param().gamd.VmaxD[gg], &sim.param().gamd.VminD[gg],
                                           &sim.param().gamd.VmeanD[gg], &sim.param().gamd.M2D[gg], &sim.param().gamd.sigmaVD[gg]);
                        calc_gamd_std_mean(ener.gamd_potential_total[gg] - ener.gamd_dihedral_total[gg], sim.steps(), &sim.param().gamd.VmaxT[gg], &sim.param().gamd.VminT[gg],
                                           &sim.param().gamd.VmeanT[gg], &sim.param().gamd.M2T[gg], &sim.param().gamd.sigmaVT[gg]); 
                        DEBUG(1, "group " << gg << " VmeanD " << sim.param().gamd.VmeanD[gg] <<  " VmeanT " << sim.param().gamd.VmeanT[gg] << " VmaxD " << sim.param().gamd.VmaxD[gg] <<  " VmaxT " << sim.param().gamd.VmaxT[gg]);                                       
                    break;          
                default:
                    io::messages.add("Unknown functional form of gaussian accelerated md boosting potential. Should be 1 (dual acceleration), 2 (total potential energy acceleration), or 3 (dihedral acceleration)",
                                     "Forcefield", io::message::critical);
                } 
            break; // end case simulation::cmd_search

            case simulation::gamd_search:
                switch (sim.param().gamd.form){
                    case simulation::dih_boost:
                        calc_gamd_std_mean(ener.gamd_dihedral_total[gg], sim.steps(), &sim.param().gamd.VmaxD[gg], &sim.param().gamd.VminD[gg],
                                           &sim.param().gamd.VmeanD[gg], &sim.param().gamd.M2D[gg], &sim.param().gamd.sigmaVD[gg]);

                        if (calc_gamd_E_K(sim.param().gamd.thresh, sim.param().gamd.dihstd, sim.param().gamd.VmaxD[gg], sim.param().gamd.VminD[gg], sim.param().gamd.VmeanD[gg],
                            sim.param().gamd.sigmaVD[gg], &sim.param().gamd.k0D[gg], &sim.param().gamd.kD[gg], &sim.param().gamd.ED[gg])){
                                io::messages.add("gamd: k0 < 0 or k0 > 1, switched to lower bound threshold ","Forcefield", io::message::warning);
                                sim.param().gamd.thresh = simulation::lower_bound;
                                }
                        break;

                    case simulation::tot_boost:
                        calc_gamd_std_mean(ener.gamd_potential_total[gg], sim.steps(), &sim.param().gamd.VmaxT[gg], &sim.param().gamd.VminT[gg],
                                           &sim.param().gamd.VmeanT[gg], &sim.param().gamd.M2T[gg], &sim.param().gamd.sigmaVT[gg]);

                        if (calc_gamd_E_K(sim.param().gamd.thresh, sim.param().gamd.totstd, sim.param().gamd.VmaxT[gg], sim.param().gamd.VminT[gg], sim.param().gamd.VmeanT[gg],
                            sim.param().gamd.sigmaVT[gg], &sim.param().gamd.k0T[gg], &sim.param().gamd.kT[gg], &sim.param().gamd.ET[gg])){
                                io::messages.add("gamd: k0 < 0 or k0 > 1, switched to lower bound threshold ","Forcefield", io::message::warning);
                                sim.param().gamd.thresh = simulation::lower_bound;
                                }
                        break;

                    case simulation::dual_boost:
                        calc_gamd_std_mean(ener.gamd_dihedral_total[gg], sim.steps(), &sim.param().gamd.VmaxD[gg], &sim.param().gamd.VminD[gg],
                                           &sim.param().gamd.VmeanD[gg], &sim.param().gamd.M2D[gg], &sim.param().gamd.sigmaVD[gg]);

                        calc_gamd_std_mean(ener.gamd_potential_total[gg] - ener.gamd_dihedral_total[gg], sim.steps(), &sim.param().gamd.VmaxT[gg], &sim.param().gamd.VminT[gg],
                                           &sim.param().gamd.VmeanT[gg], &sim.param().gamd.M2T[gg], &sim.param().gamd.sigmaVT[gg]);

                        if (calc_gamd_E_K(sim.param().gamd.thresh, sim.param().gamd.dihstd, sim.param().gamd.VmaxD[gg], sim.param().gamd.VminD[gg], sim.param().gamd.VmeanD[gg],
                            sim.param().gamd.sigmaVD[gg], &sim.param().gamd.k0D[gg], &sim.param().gamd.kD[gg], &sim.param().gamd.ED[gg])){
                                io::messages.add("gamd: k0 < 0 or k0 > 1, switched to lower bound threshold ","Forcefield", io::message::warning);
                                sim.param().gamd.thresh = simulation::lower_bound;
                                }

                        if (calc_gamd_E_K(sim.param().gamd.thresh, sim.param().gamd.totstd, sim.param().gamd.VmaxT[gg], sim.param().gamd.VminT[gg], sim.param().gamd.VmeanT[gg],
                            sim.param().gamd.sigmaVT[gg], &sim.param().gamd.k0T[gg], &sim.param().gamd.kT[gg], &sim.param().gamd.ET[gg])){
                                io::messages.add("gamd: k0 < 0 or k0 > 1, switched to lower bound threshold ","Forcefield", io::message::warning);
                                sim.param().gamd.thresh = simulation::lower_bound;
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
      for (unsigned int accelgroup = 1; accelgroup < sim.param().gamd.agroups;  accelgroup++){

          switch (sim.param().gamd.form)
          {
              case simulation::dih_boost:
              {
                  double VE = ener.gamd_dihedral_total[accelgroup] - sim.param().gamd.ED[accelgroup];
                  //if V < E apply boost
                  if (VE < 0){
                      prefactor = (sim.param().gamd.kD[accelgroup] * VE) + 1;
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
                  double VE = ener.gamd_potential_total[accelgroup] - sim.param().gamd.ET[accelgroup];
                  //if V < E apply boost
                  if (VE < 0){
                      prefactor = (sim.param().gamd.kT[accelgroup] * VE) + 1;
                      ener.gamd_DV[accelgroup] = prefactor * VE/2; // 
                      for (unsigned int accelgroup2 = 0; accelgroup2 < sim.param().gamd.agroups;  accelgroup2++){
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
                  double VET = ener.gamd_potential_total[accelgroup] - sim.param().gamd.ET[accelgroup] - ener.gamd_dihedral_total[accelgroup];
                  //if V < E apply boost
                  if (VET < 0){
                      prefactor = (sim.param().gamd.kT[accelgroup] * VET) + 1;
                      ener.gamd_DV[accelgroup] = prefactor * VET/2; // 
                      for (unsigned int accelgroup2 = 0; accelgroup2 < sim.param().gamd.agroups;  accelgroup2++){
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
                  double VED = ener.gamd_dihedral_total[accelgroup] - sim.param().gamd.ED[accelgroup];
                  //if V < E apply boost
                  if (VED < 0){
                      prefactor = (sim.param().gamd.kD[accelgroup] * VED) + 1;
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
  m_timer.stop();
  return 0;
};
// add init in which the initial acceleration values are calculated if needed

// helper gamd methods
void algorithm::GAMD::calc_gamd_std_mean(double V, int step, double *Vmax, double *Vmin, double *Vmean, double *M2, double *sigmaV){
        double delta;
        if(*Vmax == 0 && *Vmin == 0){
          *Vmin = V;
          *Vmax = V;
        }
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
