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
     if (sim.param().gamd.search != simulation::cmd_search && sim.param().gamd.ntisearch == 0){
      for (unsigned int gg = 1; gg < sim.param().gamd.igroups; gg++){
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
        io::messages.add("GAMD search should be started from an cMD search with ntisearch = 0. With ntisearch = 1 the statistics will be estimated only from the actual GAMD search run",
         "Forcefield", io::message::warning);
    }
    return 0;
 }


/**
 * GAMD step old.
 */

int algorithm::GAMD
::apply(topology::Topology & topo,
        configuration::Configuration & conf,
        simulation::Simulation &sim)
 {
  m_timer.start(sim);
  DEBUG(5, "GAMD: algorithm init");
  configuration::Energy ener = conf.current().energies;
  std::vector<int> used_chargegroups;
  unsigned int num_groups = unsigned(ener.bond_energy.size());
  unsigned int num_atoms = topo.num_atoms();
  double gamd_energy;
  DEBUG(5, "GAMD: Interactions calculated now calculate acceleration");

  // total energies have been computed now calculate acceleration
  if (sim.steps() >= sim.param().gamd.equilibration){
      sim.param().gamd.stepsdone += 1;
      for (unsigned int gg = 1; gg < sim.param().gamd.igroups; gg++){
          // This part reset the statistics so that they can be obtained form the new simulation, this part is necessary every time that one changes form search mode
          // The idea is that because after each search mode a short equilibration is performed, sim.param().gamd.equilibration will be higher than 0
          if (sim.param().gamd.equilibration == sim.steps() && sim.param().gamd.equilibration > 0){
             // reset parameters after equilibration
             sim.param().gamd.VmaxT[gg] = sim.param().gamd.VminT[gg] = sim.param().gamd.VmeanT[gg];
             sim.param().gamd.VmaxD[gg] = sim.param().gamd.VminD[gg] = sim.param().gamd.VmeanD[gg];
             sim.param().gamd.M2T[gg] = sim.param().gamd.M2D[gg] = sim.param().gamd.sigmaVT[gg] = sim.param().gamd.sigmaVD[gg] = 0.0;
          }
          DEBUG(15, "GAMD KT " << sim.param().gamd.kT[gg]);
          DEBUG(15, "GAMD KD " << sim.param().gamd.kD[gg]);
          DEBUG(15, "GAMD ET " << sim.param().gamd.ET[gg]);
          DEBUG(15, "GAMD ED " << sim.param().gamd.ED[gg]);
          DEBUG(15, "GAMD VmaxT " << sim.param().gamd.VmaxT[gg]);
          DEBUG(15, "GAMD VmaxD " << sim.param().gamd.VmaxD[gg]);
          DEBUG(15, "GAMD VminT " << sim.param().gamd.VminT[gg]);
          DEBUG(15, "GAMD VminD " << sim.param().gamd.VminD[gg]);
          DEBUG(15, "GAMD VmeanD " << sim.param().gamd.VmeanD[gg]);
          DEBUG(15, "GAMD VmeanT " << sim.param().gamd.VmeanT[gg]);
          DEBUG(15, "GAMD steps " << sim.param().gamd.stepsdone);
          switch (sim.param().gamd.search){
            case simulation::cmd_search:
                switch (sim.param().gamd.form)
                {
                    case simulation::dih_boost:
                        // reset the search values of gamd each, gamd_window steps
                        if (sim.param().gamd.gamd_window > 0 && sim.param().gamd.stepsdone % sim.param().gamd.gamd_window == 0){
                            if (ener.gamd_dihedral_total[gg] > sim.param().gamd.VmaxD[gg]) sim.param().gamd.VmaxD[gg] = ener.gamd_dihedral_total[gg];
                            if (ener.gamd_dihedral_total[gg] < sim.param().gamd.VminD[gg]) sim.param().gamd.VminD[gg] = ener.gamd_dihedral_total[gg];
                            sim.param().gamd.M2D[gg] = sim.param().gamd.sigmaVD[gg] = 0.0;
                            sim.param().gamd.VmeanD[gg] = ener.gamd_dihedral_total[gg];
                        }
                        else{  
                        calc_gamd_std_mean(ener.gamd_dihedral_total[gg], sim.param().gamd.stepsdone, &sim.param().gamd.VmaxD[gg], &sim.param().gamd.VminD[gg],
                                           &sim.param().gamd.VmeanD[gg], &sim.param().gamd.M2D[gg], &sim.param().gamd.sigmaVD[gg]);
                        }
                        break;

                    case simulation::tot_boost:
                        // reset the search values of gamd each, gamd_window steps
                        if (sim.param().gamd.gamd_window > 0 && sim.param().gamd.stepsdone % sim.param().gamd.gamd_window == 0){
                            if (ener.gamd_potential_total[gg] + ener.gamd_dihedral_total[gg] > sim.param().gamd.VmaxT[gg]){
                                 sim.param().gamd.VmaxT[gg] = ener.gamd_potential_total[gg] + ener.gamd_dihedral_total[gg];
                            }
                            if (ener.gamd_potential_total[gg] + ener.gamd_dihedral_total[gg] < sim.param().gamd.VminT[gg]){
                                 sim.param().gamd.VminT[gg] = ener.gamd_potential_total[gg] + ener.gamd_dihedral_total[gg];
                            }
                            sim.param().gamd.M2T[gg] = sim.param().gamd.sigmaVT[gg] = 0.0;
                            sim.param().gamd.VmeanT[gg] = ener.gamd_potential_total[gg] + ener.gamd_dihedral_total[gg];
                        }
                        else{  
                        calc_gamd_std_mean(ener.gamd_potential_total[gg] + ener.gamd_dihedral_total[gg], sim.param().gamd.stepsdone, &sim.param().gamd.VmaxT[gg], &sim.param().gamd.VminT[gg],
                                           &sim.param().gamd.VmeanT[gg], &sim.param().gamd.M2T[gg], &sim.param().gamd.sigmaVT[gg]);
                        }  
                        break;

                    case simulation::dual_boost:
                        // reset the search values of gamd each, gamd_window steps
                        if (sim.param().gamd.gamd_window > 0 && sim.param().gamd.stepsdone % sim.param().gamd.gamd_window == 0){
                            if (ener.gamd_dihedral_total[gg] > sim.param().gamd.VmaxD[gg]) sim.param().gamd.VmaxD[gg] = ener.gamd_dihedral_total[gg];
                            if (ener.gamd_dihedral_total[gg] < sim.param().gamd.VminD[gg]) sim.param().gamd.VminD[gg] = ener.gamd_dihedral_total[gg];
                            sim.param().gamd.M2D[gg] = sim.param().gamd.sigmaVD[gg] = 0.0;
                            sim.param().gamd.VmeanD[gg] = ener.gamd_dihedral_total[gg];
                        }
                        else{  
                        calc_gamd_std_mean(ener.gamd_dihedral_total[gg], sim.param().gamd.stepsdone, &sim.param().gamd.VmaxD[gg], &sim.param().gamd.VminD[gg],
                                           &sim.param().gamd.VmeanD[gg], &sim.param().gamd.M2D[gg], &sim.param().gamd.sigmaVD[gg]);
                        }
                        // reset the search values of gamd each, gamd_window steps
                        if (sim.param().gamd.gamd_window > 0 && sim.param().gamd.stepsdone % sim.param().gamd.gamd_window == 0){
                            if (ener.gamd_potential_total[gg] > sim.param().gamd.VmaxT[gg]) sim.param().gamd.VmaxT[gg] = ener.gamd_potential_total[gg];
                            if (ener.gamd_potential_total[gg] < sim.param().gamd.VminT[gg]) sim.param().gamd.VminT[gg] = ener.gamd_potential_total[gg];
                            sim.param().gamd.M2T[gg] = sim.param().gamd.sigmaVT[gg] = 0.0;
                            sim.param().gamd.VmeanT[gg] = ener.gamd_potential_total[gg];
                        }
                        else{  
                        calc_gamd_std_mean(ener.gamd_potential_total[gg], sim.param().gamd.stepsdone, &sim.param().gamd.VmaxT[gg], &sim.param().gamd.VminT[gg],
                                           &sim.param().gamd.VmeanT[gg], &sim.param().gamd.M2T[gg], &sim.param().gamd.sigmaVT[gg]); 
                        DEBUG(10, "group " << gg << " VmeanD " << sim.param().gamd.VmeanD[gg] <<  " VmeanT " << sim.param().gamd.VmeanT[gg] << " VmaxD " << sim.param().gamd.VmaxD[gg] <<  " VmaxT " << sim.param().gamd.VmaxT[gg]);
                        }                                     
                    break;          
                default:
                    io::messages.add("Unknown functional form of gaussian accelerated md boosting potential. Should be 1 (dual acceleration), 2 (total potential energy acceleration), or 3 (dihedral acceleration)",
                                     "Forcefield", io::message::critical);
                } 
            break; // end case simulation::cmd_search

            case simulation::gamd_search:
                switch (sim.param().gamd.form){
                    case simulation::dih_boost:
                        // reset the search values of gamd each, gamd_window steps
                        if (sim.param().gamd.gamd_window > 0 && sim.param().gamd.stepsdone%sim.param().gamd.gamd_window == 0){
                            if (ener.gamd_dihedral_total[gg] > sim.param().gamd.VmaxD[gg]) sim.param().gamd.VmaxD[gg] = ener.gamd_dihedral_total[gg];
                            if (ener.gamd_dihedral_total[gg] < sim.param().gamd.VminD[gg]) sim.param().gamd.VminD[gg] = ener.gamd_dihedral_total[gg];
                            sim.param().gamd.M2D[gg] = sim.param().gamd.sigmaVD[gg] = 0.0;
                            sim.param().gamd.VmeanD[gg] = ener.gamd_dihedral_total[gg];
                        }
                        else{
                            calc_gamd_std_mean(ener.gamd_dihedral_total[gg], sim.param().gamd.stepsdone, &sim.param().gamd.VmaxD[gg], &sim.param().gamd.VminD[gg],
                                           &sim.param().gamd.VmeanD[gg], &sim.param().gamd.M2D[gg], &sim.param().gamd.sigmaVD[gg]);
                        }

                        if (calc_gamd_E_K(sim.param().gamd.thresh, sim.param().gamd.dihstd, sim.param().gamd.VmaxD[gg], sim.param().gamd.VminD[gg], sim.param().gamd.VmeanD[gg],
                            sim.param().gamd.sigmaVD[gg], &sim.param().gamd.k0D[gg], &sim.param().gamd.kD[gg], &sim.param().gamd.ED[gg])){
                                io::messages.add("gamd: k0 < 0 or k0 > 1, switched to lower bound threshold ","Forcefield", io::message::warning);
                                sim.param().gamd.thresh = simulation::lower_bound;
                                }
                        break;

                    case simulation::tot_boost:
                        // reset the search values of gamd each, gamd_window steps
                        if (sim.param().gamd.gamd_window > 0 && sim.param().gamd.stepsdone%sim.param().gamd.gamd_window == 0){
                            if (ener.gamd_potential_total[gg] + ener.gamd_dihedral_total[gg] > sim.param().gamd.VmaxT[gg]){
                                 sim.param().gamd.VmaxT[gg] = ener.gamd_potential_total[gg] + ener.gamd_dihedral_total[gg];
                            }
                            if (ener.gamd_potential_total[gg] + ener.gamd_dihedral_total[gg] < sim.param().gamd.VminT[gg]){
                                 sim.param().gamd.VminT[gg] = ener.gamd_potential_total[gg] + ener.gamd_dihedral_total[gg];
                            }
                            sim.param().gamd.M2T[gg] = sim.param().gamd.sigmaVT[gg] = 0.0;
                            sim.param().gamd.VmeanT[gg] = ener.gamd_potential_total[gg] + ener.gamd_dihedral_total[gg];
                        }
                        else{
                        calc_gamd_std_mean(ener.gamd_potential_total[gg] + ener.gamd_dihedral_total[gg], sim.param().gamd.stepsdone, &sim.param().gamd.VmaxT[gg], &sim.param().gamd.VminT[gg],
                                           &sim.param().gamd.VmeanT[gg], &sim.param().gamd.M2T[gg], &sim.param().gamd.sigmaVT[gg]);
                        }
                        if (calc_gamd_E_K(sim.param().gamd.thresh, sim.param().gamd.totstd, sim.param().gamd.VmaxT[gg], sim.param().gamd.VminT[gg], sim.param().gamd.VmeanT[gg],
                            sim.param().gamd.sigmaVT[gg], &sim.param().gamd.k0T[gg], &sim.param().gamd.kT[gg], &sim.param().gamd.ET[gg])){
                                io::messages.add("gamd: k0 < 0 or k0 > 1, switched to lower bound threshold ","Forcefield", io::message::warning);
                                sim.param().gamd.thresh = simulation::lower_bound;
                                }
                        break;

                    case simulation::dual_boost:
                        // reset the search values of gamd each, gamd_window steps
                        if (sim.param().gamd.gamd_window > 0 && sim.param().gamd.stepsdone%sim.param().gamd.gamd_window == 0){
                            if (ener.gamd_dihedral_total[gg] > sim.param().gamd.VmaxD[gg]) sim.param().gamd.VmaxD[gg] = ener.gamd_dihedral_total[gg];
                            if (ener.gamd_dihedral_total[gg] < sim.param().gamd.VminD[gg]) sim.param().gamd.VminD[gg] = ener.gamd_dihedral_total[gg];
                            sim.param().gamd.M2D[gg] = sim.param().gamd.sigmaVD[gg] = 0.0;
                            sim.param().gamd.VmeanD[gg] = ener.gamd_dihedral_total[gg];
                        }
                        else{
                        calc_gamd_std_mean(ener.gamd_dihedral_total[gg],sim.param().gamd.stepsdone, &sim.param().gamd.VmaxD[gg], &sim.param().gamd.VminD[gg],
                                           &sim.param().gamd.VmeanD[gg], &sim.param().gamd.M2D[gg], &sim.param().gamd.sigmaVD[gg]);
                        }
                        // reset the search values of gamd each, gamd_window steps
                        if (sim.param().gamd.gamd_window > 0 && sim.param().gamd.stepsdone%sim.param().gamd.gamd_window == 0){
                            if (ener.gamd_potential_total[gg] > sim.param().gamd.VmaxT[gg]) sim.param().gamd.VmaxT[gg] = ener.gamd_potential_total[gg];
                            if (ener.gamd_potential_total[gg] < sim.param().gamd.VminT[gg]) sim.param().gamd.VminT[gg] = ener.gamd_potential_total[gg];
                            sim.param().gamd.M2T[gg] = sim.param().gamd.sigmaVT[gg] = 0.0;
                            sim.param().gamd.VmeanT[gg] = ener.gamd_potential_total[gg];
                        }
                        else{
                        calc_gamd_std_mean(ener.gamd_potential_total[gg], sim.param().gamd.stepsdone, &sim.param().gamd.VmaxT[gg], &sim.param().gamd.VminT[gg],
                                           &sim.param().gamd.VmeanT[gg], &sim.param().gamd.M2T[gg], &sim.param().gamd.sigmaVT[gg]);
                        }

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
                // do nothing 
            break; // end case no_search
            default:
                io::messages.add("Unknown functional form of gaussian accelerated md search. Should be 0 (gamd with fixed parameters), 1 (cmd search), or 2 (gamd search)",
                                 "Forcefield", io::message::critical);;
        } // end switch search form
     } // for loop over gamd groups
  } // endif 

  // statistics updated now accelerate
  DEBUG(5, "GAMD: Accelerating");
  double prefactor_t;
  double prefactor_h;
  if (sim.param().gamd.search != simulation::cmd_search){
      for (unsigned int gg = 1; gg < sim.param().gamd.igroups;  gg++){

          switch (sim.param().gamd.form)
          {
              case simulation::dih_boost:
              {
                  double VE = ener.gamd_dihedral_total[gg] - sim.param().gamd.ED[gg];
                  //if V < E apply boost
                  if (VE < 0){
                      //check prefactor
                      prefactor_h = (sim.param().gamd.kD[gg] * VE);
                      conf.current().energies.gamd_DV[gg] = prefactor_h * VE/2;
                      // loop over atoms
                       for (unsigned int atom=0; atom < num_atoms; atom++){
                              conf.current().force(atom) += conf.special().gamd.dihe_force[gg](atom) * (prefactor_h);
                      } // end loop over atoms
                      // to virial
                    for (int a = 0; a < 3; ++a) {
                        for (int b = 0; b < 3; ++b) {
                            conf.current().virial_tensor(b, a) += conf.special().gamd.virial_tensor_dihe[gg](b, a) * (prefactor_h);
                        }
                    } // end virial 
                  } // endif
                  break;
              }

          case simulation::tot_boost:
          {
                  double VE = ener.gamd_potential_total[gg] + ener.gamd_dihedral_total[gg] - sim.param().gamd.ET[gg];
                  //if V < E apply boost
                  if (VE < 0){
                      prefactor_t = (sim.param().gamd.kT[gg] * VE);
                      conf.current().energies.gamd_DV[gg] = prefactor_t * VE/2; 
                      // First do the dihedral term which is stored in a different array
                      // loop over atoms
                      for (unsigned int atom=0; atom < num_atoms; atom++){
                            conf.current().force(atom) += conf.special().gamd.dihe_force[gg](atom) * (prefactor_t);
                            conf.current().force(atom) += conf.special().gamd.total_force[gg](atom) * (prefactor_t);
                      } // end loop over atoms
                      // to virial
                      for (int a = 0; a < 3; ++a) {
                        for (int b = 0; b < 3; ++b) {
                            conf.current().virial_tensor(b, a) += conf.special().gamd.virial_tensor_dihe[gg](b, a) * (prefactor_t);
                            conf.current().virial_tensor(b, a) += conf.special().gamd.virial_tensor[gg](b, a) * (prefactor_t);
                        }
                      } // end virial 
                  } // endif
                  break;
            }

          case simulation::dual_boost:
          {
                  // first total potential
                  double VET = ener.gamd_potential_total[gg] - sim.param().gamd.ET[gg];
                  //if V < E apply boost
                  if (VET < 0){
                      prefactor_t = (sim.param().gamd.kT[gg] * VET);
                      conf.current().energies.gamd_DV[gg] = prefactor_t * VET/2; 
                      for (unsigned int atom=0; atom < num_atoms; atom++){
                          conf.current().force(atom) += conf.special().gamd.total_force[gg](atom) * (prefactor_t);
                        }// end loop over atoms
                        // to virial 
                      for (int a = 0; a < 3; ++a) {
                        for (int b = 0; b < 3; ++b) {
                            conf.current().virial_tensor(b, a) += conf.special().gamd.virial_tensor[gg](b, a) * (prefactor_t);
                        }
                      } // end virial 
                  } // endif

                  // dihedral term
                  double VED = ener.gamd_dihedral_total[gg] - sim.param().gamd.ED[gg];
                  //if V < E apply boost
                  if (VED < 0){
                      prefactor_h = (sim.param().gamd.kD[gg] * VED);
                      conf.current().energies.gamd_DV[gg] = prefactor_h * VED/2;
                      // loop over atoms
                       for (unsigned int atom=0; atom < num_atoms; atom++){
                              conf.current().force(atom) += conf.special().gamd.dihe_force[gg](atom) * (prefactor_h);
                              //int chargegroup = topo.atom_energy_group()[atom];
                        } // end loop over atoms
                        for (int a = 0; a < 3; ++a) {
                            for (int b = 0; b < 3; ++b) {
                                conf.current().virial_tensor(b, a) += conf.special().gamd.virial_tensor_dihe[gg](b, a) * (prefactor_h);
                            }
                        } // end virial 
                    } // endif
                  break; 
            }         
          default:
            io::messages.add("Unknown functional form of gaussian accelerated md boosting potential. Should be 1 (dual acceleration), 2 (total potential energy acceleration), or 3 (dihedral acceleration)",
                             "Forcefield", io::message::critical);
          } // end switch
        DEBUG(5, "for acceleration group " << gg << " ET: " << sim.param().gamd.ET[gg] << " ED: " << sim.param().gamd.ED[gg]
              << " KT: " << sim.param().gamd.kT[gg]<< " KD: " << sim.param().gamd.kD[gg]);
        //TO DO: migrate E and K to energys to avoid having to copy them.
        conf.current().energies.gamd_ED[gg] = sim.param().gamd.ED[gg];
        conf.current().energies.gamd_ET[gg] = sim.param().gamd.ET[gg];
        conf.current().energies.gamd_KD[gg] = sim.param().gamd.kD[gg];
        conf.current().energies.gamd_KT[gg] = sim.param().gamd.kT[gg];
        
      } // end loop over acceleration groups
  } // end if
  m_timer.stop();
  return 0;
};


// helper gamd methods
void algorithm::GAMD::calc_gamd_std_mean(double V, int step, double *Vmax, double *Vmin, double *Vmean, long double *M2, long double *sigmaV){
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

int algorithm::GAMD::calc_gamd_E_K(simulation::gamd_thresh_enum Eform, double sigma0, double Vmax, double Vmin, double Vmean,long double sigmaV,
                   double *k0, double *k, double *E)
        {
                if (Vmax == 0.0 && Vmin == 0.0 && Vmean == 0.0){
                    return 0;
                }
                
                if (Eform == simulation::upper_bound){
                        *k0 = (1 - sigma0/sigmaV) * (Vmax-Vmin) / (Vmean-Vmin);
                        // k0 should be between 0 and 1 if not switch to lower bound energy threshold
                        if(*k0 > 1.0 || *k0 <= 0.0){
                                *k0 = (sigma0 / sigmaV) * (Vmax-Vmin) / (Vmax-Vmean);
                                if (*k0 > 1)  *k0 = 1.0;
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
                                if (*k0 > 1)  *k0 = 1.0;
                                *E = Vmax;
                                *k = *k0 / (Vmax - Vmin);
                                return 0;  

                }
        };