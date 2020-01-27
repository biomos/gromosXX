/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   replica_MPI_slave.cpp
 * Author: bschroed
 * 
 * Created on November 8, 2019, 12:48 PM
 */
#ifdef XXMPI
    #include <mpi.h>
#endif

#include <stdheader.h>

#include "replica_MPI_slave.h"
#include <io/argument.h>
#include <util/error.h>
#include <math/volume.h>

#undef MODULE
#undef SUBMODULE
#define MODULE util
#define SUBMODULE replica_exchange

util::replica_MPI_Slave::replica_MPI_Slave(io::Argument _args, int cont, int globalThreadID, 
         simulation::mpi_control_struct replica_mpi_control) : replica_Interface(globalThreadID, replica_mpi_control, _args){
#ifdef XXMPI
    /**
     * Build up replica MPI Slave - reads in the input and tries to link the replica to the Master?
     * 
     * @param _args
     * @param cont
     * @param _ID
     * @param _globalThreadID
     * @param simulation_globalThreadID
     * @param simulation_ID
     * @param simulation_num_threads
     */
    
    MPI_DEBUG(4, "replica_MPI_SLAVE "<< globalThreadID <<":Constructor:\t  "<< globalThreadID <<":\t START");

    /**
     * READ INPUT
     */
    // do continuation run?
    // change name of input coordinates
    if(cont == 1){
      std::multimap< std::string, std::string >::iterator it = args.lower_bound(("conf"));
      size_t pos = (*it).second.find_last_of(".");
      std::stringstream tmp;
      tmp << "_" << (simulation_ID+1);
      (*it).second.insert(pos, tmp.str());
    }
    
    MPI_DEBUG(5, "replica_MPI_SLAVE "<< globalThreadID <<":Constructor:\t  "<< globalThreadID <<":\t start read in");
    //build up system:
    sim.mpi = true;
    sim.mpi_control = replica_mpi_control;

    if (io::read_input(args, topo, conf, sim, md, *os, true)) { 
      io::messages.display(*os);
      std::cerr << "\nErrors during initialization!\n" << std::endl;
      #ifdef XXMPI
        MPI_Abort(MPI_COMM_WORLD, E_INPUT_ERROR);
      #endif
    }
    MPI_DEBUG(5, "replica_MPI_SLAVE "<< globalThreadID <<":Constructor:\t  "<< globalThreadID <<":\t REad in input already");

    MPI_DEBUG(5, "replica_MPI_SLAVE "<< globalThreadID << "  set replica_params ");

    // check whether we have nonbonded interactions
    do_nonbonded = (sim.param().force.nonbonded_crf || 
                         sim.param().force.nonbonded_vdw);
    //DEBUG(5, "Slave "<< globalThreadID << " do_non_bonded? "<<do_nonbonded);

    // let's get the forcefield
    ff = dynamic_cast<interaction::Forcefield *>(md.algorithm("Forcefield"));
    
    if (do_nonbonded && ff == NULL){
      std::cerr << "globalThreadID: "<< globalThreadID<<"\tMPI slave: could not access forcefield\n\t(internal error)" << std::endl;
      MPI::Finalize();
      throw "ReplicaSlave: \ttINITIALisation error in ff is NULL of slave of simulationID "+ simulation_ID;
    }
    
    nb = ff->interaction("NonBonded");
    if (do_nonbonded && nb == NULL){
      std::cerr << "globalThreadID: "<< globalThreadID<<"\tMPI slave: could not get NonBonded interactions from forcefield"
		<< "\n\t(internal error)"
		<< std::endl;
      MPI::Finalize();
      throw "ReplicaSlave: \tINITIALisation error in ff nbs are NULL of slave of simulationID "+ simulation_ID;
    }

    // get shake and check whether we do it for solvent
    do_shake = ((sim.param().system.npm && sim.param().constraint.solute.algorithm == simulation::constr_shake) ||
      (sim.param().constraint.solvent.algorithm == simulation::constr_shake && sim.param().system.nsm))
        && !sim.param().analyze.no_constraints;

    // for stochastic dynamics simulation we need to call SHAKE twice
    do_shake_twice = sim.param().stochastic.sd && !sim.param().analyze.no_constraints;

    shake = dynamic_cast<algorithm::Shake *>(md.algorithm("Shake"));
    if (do_shake && shake == NULL) {
        std::cerr << "globalThreadID: "<< globalThreadID<<"\tMPI slave: could not get Shake algorithm from MD sequence."
                << "\n\t(internal error)"
                << std::endl;
      MPI::Finalize();
      throw "ReplicaSlave: \tINITIALisation error in shake algorithm of slave of simulationID "+ simulation_ID;
    }
    //DEBUG(5, "Slave "<< globalThreadID << " do_shake? "<<do_shake);


    // get m_shake and check whether we do it for solvent
    do_m_shake = sim.param().system.nsm &&
      sim.param().constraint.solvent.algorithm == simulation::constr_m_shake;

    m_shake = dynamic_cast<algorithm::M_Shake *>(md.algorithm("M_Shake"));
    if (do_m_shake && m_shake == NULL) {
        std::cerr << "globalThreadID: "<< globalThreadID<<"\tMPI slave: could not get M_Shake algorithm from MD sequence."
                << "\n\t(internal error)"
                << std::endl;
      MPI::Finalize();
      throw "ReplicaSlave: \tINITIALisation error in m_Shake algorithm of slave of simulationID "+ simulation_ID;
    }
    
    // get chemical monte carlo
    do_cmc = sim.param().montecarlo.mc;
    monte_carlo = dynamic_cast<algorithm::Monte_Carlo *>(md.algorithm("MonteCarlo"));
    if(do_cmc && monte_carlo == NULL){
      std::cerr << "globalThreadID: "<< globalThreadID<<"\tMPI slave: could not get Monte Carlo algorithm from MD sequence."
              << "\n\t(internal error)"
              << std:: endl;
      MPI::Finalize();
      throw "ReplicaSlave: \tINITIALisation error in Montecarlo algorithm of slave of simulationID "+ simulation_ID;
    }
    MPI_DEBUG(4, "replica_MPI_SLAVE "<< globalThreadID <<":Constructor:\t  "<< globalThreadID <<":\t DONE");

#else
    throw "Can not construct Replica_MPI as MPI is not enabled!";
#endif
}


util::replica_MPI_Slave::~replica_MPI_Slave() {
  delete rng;
  delete traj;
  delete os;
}

void util::replica_MPI_Slave::run_MD(){
    int error;
    int next_step = 0 ;

    while ((unsigned int)(sim.steps()) < stepsPerRun + curentStepNumber) {
      // run a step
       DEBUG(4, "slave " << globalThreadID << " waiting for master");
      if (do_nonbonded && (error = nb->calculate_interactions(topo, conf, sim)) != 0){
        std::cout << "Rank: "<< globalThreadID<<"\tMPI slave " << globalThreadID << ": error in nonbonded calculation!\n" << std::endl;
      }

      // DEBUG(10, "slave " << globalThreadID << " step done");
      // (*os) << "step done (it really worked?)" << std::endl;

      if (do_cmc && (error = monte_carlo->apply(topo, conf, sim)) != 0){
        std::cout << "Rank: "<< globalThreadID<<"\tMPI slave " << globalThreadID << ": error in Monte Carlo algorithm!\n" 
                << std::endl;
      }

      if (do_shake && (error = shake->apply(topo, conf, sim)) != 0) {
        std::cout << "Rank: "<< globalThreadID<<"\tMPI slave " << globalThreadID << ": error in Shake algorithm!\n" << std::endl;
      }

      // if stochastic dynamics simulation, then expect second call to SHAKE
      if (do_shake && do_shake_twice && (error = shake->apply(topo, conf, sim)) != 0) {
        std::cout << "Rank: "<< globalThreadID<<"\tMPI slave " << globalThreadID << ": error in Shake algorithm on second call!\n" << std::endl;
      }

      if (do_m_shake && (error = m_shake->apply(topo, conf, sim)) != 0) {
        std::cout << "Rank: "<< globalThreadID<<"\tMPI slave " << globalThreadID << ": error in M-Shake algorithm!\n" << std::endl;
      }


      MPI::COMM_WORLD.Bcast(&next_step, 1, MPI::INT, , sim.mpi_control.simulationMasterThreadID, sim.mpi_control.simulationCOMM);

      if (next_step == 2) {
        (*os) << "globalThreadID: "<< globalThreadID<<"\tMessage from master: Steepest descent: minimum reached." << std::endl;
        error = 0;
        break;
      }
      else if (!next_step) {
        (*os) << "Rank: "<< globalThreadID<<"\tThere was an error in the master. Check output file for details."
              << "Exiting from MD main loop." << std::endl;
        error = 1;
        break;
      }

      ++sim.steps();
      sim.time() = sim.param().step.t0 + sim.steps() * sim.time_step_size();

    }

    if (error){
      (*os) << "Rank: "<< globalThreadID<<"\t\nErrors encountered during run in simulation "<< simulation_ID << " in Slave Thread "<< simulation_globalThreadID <<" - check above!\n" << std::endl;
      MPI_Abort(MPI_COMM_WORLD, E_INPUT_ERROR);
    }
    else{
      (*os) << "Rank: "<< globalThreadID<<"\n" GROMOSXX " MPI slave " << globalThreadID << "  simulation "<< simulation_ID << " in Slave Thread "<< simulation_globalThreadID << "finished successfully\n" << std::endl;
    }
}
    