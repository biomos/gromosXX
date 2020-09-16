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
  configuration::Energy ener = conf.current().energies;
  simulation::Parameter::gamd_struct params = sim.param().gamd;
  std::vector<int> used_chargegroups;
  unsigned int num_groups = unsigned(ener.bond_energy.size());
  // to do: add asserts, move energie calculation to energie calculate totals
  // Calculate the energies for each acceleration group
  // loop over all the acceleration groups
  for (unsigned int gamdgroup = 0; gamdgroup < topo.gamd_atoms().size(); gamdgroup++){
          // loop over the atoms to be accelerated in each group
          for (unsigned int atom : topo.gamd_atoms()[gamdgroup]){
                  // check if the charge group has already been used
                  std::vector<int>::iterator it = std::find(used_chargegroups.begin(), used_chargegroups.end(), atom);
                  if (it == used_chargegroups.end()){
                        unsigned int chargegroup = topo.atom_energy_group()[atom];
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
                        for(unsigned int j=0; j<num_groups; j++){
                                if (j != chargegroup){
                                        ener.gamd_potential_total[gamdgroup] +=  ener.lj_energy[chargegroup][j];
                                        ener.gamd_potential_total[gamdgroup] +=  ener.crf_energy[chargegroup][j];
                                }
                        }
                        //save the used charge group
                        used_chargegroups.push_back(chargegroup);
                        
                } // endif        
          } // loop over atoms 
  } // loop over acceleration groups

  // total energies have been computed now calculate acceleration
  if (sim.steps() > params.equilibration){
      for (unsigned int gg = 0; gg < topo.gamd_atoms().size(); gg++){
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
  double factor;
  double prefactor;
  if (sim.param().gamd.search != simulation::cmd_search){
      for (unsigned int gg = 0; gg < topo.gamd_atoms().size(); gg++){
          switch (sim.param().gamd.form)
          {
              case simulation::dih_boost:
                  double VE = ener.gamd_dihedral_total[gg] - params.ED[gg];
                  //if V < E apply boost
                  if (VE < 0){
                      prefactor = params.kD[gg] * params.ED[gg];
                      factor = prefactor + 1.0;
                      ener.gamd_DV[gg] = prefactor * VE/2;
                      // loop over atoms
                      for (unsigned int atom : topo.gamd_atoms()[gg]){
                              conf.current().force(atom) += conf.special().gamd.dihe_force(atom) * (factor - 1);
                              int chargegroup = topo.atom_energy_group()[atom];
                              /**
                              // to virial to do: save value for charge group??
                              for (int a = 0; a < 3; ++a) {
                                for (int b = 0; b < 3; ++b) {
                                    conf.current().virial_tensor(b, a) += conf.special().gamd.virial_tensor_dihe[chargegroup](b, a) * (prefactor -1);
                                }
                              } // end virial **/
                      } // end loop over atoms
                  } // endif
                  break;

          case simulation::tot_boost:
                  double VE = ener.gamd_potential_total[gg] - params.ET[gg];
                  //if V < E apply boost
                  if (VE < 0){
                      prefactor = params.kT[gg] * params.ET[gg];
                      factor = prefactor + 1.0;
                      ener.gamd_DV[gg] = prefactor * VE/2;
                      // loop over atoms
                      for (unsigned int atom : topo.gamd_atoms()[gg]){
                              conf.current().force(atom) *= factor;
                              int chargegroup = topo.atom_energy_group()[atom];
                              /**
                              // to virial to do: save value for charge group check charge group is not used
                              for (int a = 0; a < 3; ++a) {
                                for (int b = 0; b < 3; ++b) {
                                    conf.current().virial_tensor(b, a) += conf.special().gamd.virial_tensor[chargegroup](b, a) * (prefactor -1); 
                                }
                              } // end virial **/
                      } // end loop over atoms
                  } // endif
                  break;

          case simulation::dual_boost:
                  double VET = ener.gamd_potential_total[gg] - params.ET[gg] - ener.gamd_dihedral_total[gg];
                  //if V < E apply boost
                  if (VET < 0){
                      prefactor = params.kT[gg] * params.ET[gg];
                      factor = prefactor + 1.0;
                      ener.gamd_DV[gg] = prefactor * VET/2;
                      // loop over atoms
                      for (unsigned int atom : topo.gamd_atoms()[gg]){
                              conf.current().force(atom) *= factor;
                              int chargegroup = topo.atom_energy_group()[atom];
                              /**
                              // to virial to do: save value for charge group check charge group is not used
                              for (int a = 0; a < 3; ++a) {
                                for (int b = 0; b < 3; ++b) {
                                    conf.current().virial_tensor(b, a) += conf.special().gamd.virial_tensor[chargegroup](b, a) * (prefactor -1); 
                                }
                              } // end virial **/
                      } // end loop over atoms
                  } // endif
                  else{
                      factor = 1.0;
                  }

                  double VED = ener.gamd_dihedral_total[gg] - params.ED[gg];
                  //if V < E apply boost
                  if (VED < 0){
                      prefactor = params.kD[gg] * params.ED[gg];
                      double factor_dihedral = prefactor + 1.0;
                      ener.gamd_DV[gg] += prefactor * VED/2;
                      // loop over atoms
                      for (unsigned int atom : topo.gamd_atoms()[gg]){
                              conf.current().force(atom) += conf.special().gamd.dihe_force(atom) * (factor_dihedral - factor);
                              int chargegroup = topo.atom_energy_group()[atom];
                              /**
                              // to virial to do: save value for charge group??
                              for (int a = 0; a < 3; ++a) {
                                for (int b = 0; b < 3; ++b) {
                                    conf.current().virial_tensor(b, a) += conf.special().gamd.virial_tensor_dihe[chargegroup](b, a) * (prefactor -1);
                                }
                              } // end virial **/
                      } // end loop over atoms
                  } // endif

                  break;          
          default:
            io::messages.add("Unknown functional form of gaussian accelerated md boosting potential. Should be 1 (dual acceleration), 2 (total potential energy acceleration), or 3 (dihedral acceleration)",
                             "Forcefield", io::message::critical);
          } // end switch
      } // end loop over acceleration groups
  } // end if

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
                                return 1;
                        }
                        else{
                                *E = Vmin + (Vmax-Vmin)/(*k0);
                                return 0;
                        }

                }
                else if(Eform == simulation::lower_bound){
                                *k0 = (sigma0 / sigmaV) * (Vmax-Vmin) / (Vmax-Vmean);
                                *E = Vmax;
                                return 0;  

                }
        };