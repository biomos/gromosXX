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
#define XXMPI

util::replica_MPI_Slave::replica_MPI_Slave(io::Argument _args, int cont, int _ID, int _rank, 
        int simulation_rank, int simulation_ID, int simulation_num_threads) : replica_Interface(_ID, _rank, _args),
         simulation_ID(simulation_ID), simulation_rank(simulation_rank), simulation_num_threads(simulation_num_threads){
#ifdef XXMPI
    /**
     * Build up replica MPI Slave - reads in the input and tries to link the replica to the Master?
     * 
     * @param _args
     * @param cont
     * @param _ID
     * @param _rank
     * @param simulation_rank
     * @param simulation_ID
     * @param simulation_num_threads
     */
    
    MPI_DEBUG(4, "replica_MPI_SLAVE "<< rank <<":Constructor:\t  "<< rank <<":\t START");

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
    
    MPI_DEBUG(5, "replica_MPI_SLAVE "<< rank <<":Constructor:\t  "<< rank <<":\t start read in");
    //build up system:
    sim.mpi = true;
    
    if (io::read_input(args, topo, conf, sim, md, *os, true)) { 
      io::messages.display(*os);
      std::cerr << "\nErrors during initialization!\n" << std::endl;
      #ifdef XXMPI
        MPI_Abort(MPI_COMM_WORLD, E_INPUT_ERROR);
      #endif
    }
    MPI_DEBUG(5, "replica_MPI_SLAVE "<< rank <<":Constructor:\t  "<< rank <<":\t REad in input already");

    MPI_DEBUG(5, "replica_MPI_SLAVE "<< rank << "  set replica_params ");
    
    // adapt system to replica parameters
    maxSteps = sim.param().step.number_of_steps;  //neede?
    run = 0;  //neede?
    total_runs = sim.param().replica.trials + sim.param().replica.equilibrate;  //neede?
    steps = 0;  //neede?
    switched = 0;  //neede?

    const int numT = sim.param().replica.num_T;

    T = sim.param().replica.temperature[simulation_ID % numT];
    l = sim.param().replica.lambda[simulation_ID / numT];
    dt = sim.param().replica.dt[simulation_ID / numT];

    set_lambda();
    set_temp();
    
    MPI_DEBUG(5, "replica_MPI_SLAVE "<< rank << "  done replica_params ");

    MPI_DEBUG(5, "Slave "<< rank << " now make me a slave ");

    // check whether we have nonbonded interactions
    do_nonbonded = (sim.param().force.nonbonded_crf || 
                         sim.param().force.nonbonded_vdw);
    //DEBUG(5, "Slave "<< rank << " do_non_bonded? "<<do_nonbonded);

    // let's get the forcefield
    ff = dynamic_cast<interaction::Forcefield *>(md.algorithm("Forcefield"));
    
    if (do_nonbonded && ff == NULL){
      std::cerr << "RANK: "<< rank<<"\tMPI slave: could not access forcefield\n\t(internal error)" << std::endl;
      MPI::Finalize();
      throw "ReplicaSlave: \ttINITIALisation error in ff is NULL of slave of simulationID "+ simulation_ID;
    }
    
    nb = ff->interaction("NonBonded");
    if (do_nonbonded && nb == NULL){
      std::cerr << "RANK: "<< rank<<"\tMPI slave: could not get NonBonded interactions from forcefield"
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
        std::cerr << "RANK: "<< rank<<"\tMPI slave: could not get Shake algorithm from MD sequence."
                << "\n\t(internal error)"
                << std::endl;
      MPI::Finalize();
      throw "ReplicaSlave: \tINITIALisation error in shake algorithm of slave of simulationID "+ simulation_ID;
    }
    //DEBUG(5, "Slave "<< rank << " do_shake? "<<do_shake);


    // get m_shake and check whether we do it for solvent
    do_m_shake = sim.param().system.nsm &&
      sim.param().constraint.solvent.algorithm == simulation::constr_m_shake;

    m_shake = dynamic_cast<algorithm::M_Shake *>(md.algorithm("M_Shake"));
    if (do_m_shake && m_shake == NULL) {
        std::cerr << "RANK: "<< rank<<"\tMPI slave: could not get M_Shake algorithm from MD sequence."
                << "\n\t(internal error)"
                << std::endl;
      MPI::Finalize();
      throw "ReplicaSlave: \tINITIALisation error in m_Shake algorithm of slave of simulationID "+ simulation_ID;
    }
    
    // get chemical monte carlo
    do_cmc = sim.param().montecarlo.mc;
    monte_carlo = dynamic_cast<algorithm::Monte_Carlo *>(md.algorithm("MonteCarlo"));
    if(do_cmc && monte_carlo == NULL){
      std::cerr << "RANK: "<< rank<<"\tMPI slave: could not get Monte Carlo algorithm from MD sequence."
              << "\n\t(internal error)"
              << std:: endl;
      MPI::Finalize();
      throw "ReplicaSlave: \tINITIALisation error in Montecarlo algorithm of slave of simulationID "+ simulation_ID;
    }
    MPI_DEBUG(4, "replica_MPI_SLAVE "<< rank <<":Constructor:\t  "<< rank <<":\t DONE");

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
    sim.steps() = steps;
    sim.time() = time;
    int next_step = 0 ;

    while ((unsigned int)(sim.steps()) < maxSteps + steps) {
      // run a step
       DEBUG(4, "slave " << rank << " waiting for master");
      if (do_nonbonded && (error = nb->calculate_interactions(topo, conf, sim)) != 0){
        std::cout << "RANK: "<< rank<<"\tMPI slave " << rank << ": error in nonbonded calculation!\n" << std::endl;
      }

      // DEBUG(10, "slave " << rank << " step done");
      // (*os) << "step done (it really worked?)" << std::endl;

      if (do_cmc && (error = monte_carlo->apply(topo, conf, sim)) != 0){
        std::cout << "RANK: "<< rank<<"\tMPI slave " << rank << ": error in Monte Carlo algorithm!\n" 
                << std::endl;
      }

      if (do_shake && (error = shake->apply(topo, conf, sim)) != 0) {
        std::cout << "RANK: "<< rank<<"\tMPI slave " << rank << ": error in Shake algorithm!\n" << std::endl;
      }

      // if stochastic dynamics simulation, then expect second call to SHAKE
      if (do_shake && do_shake_twice && (error = shake->apply(topo, conf, sim)) != 0) {
        std::cout << "RANK: "<< rank<<"\tMPI slave " << rank << ": error in Shake algorithm on second call!\n" << std::endl;
      }

      if (do_m_shake && (error = m_shake->apply(topo, conf, sim)) != 0) {
        std::cout << "RANK: "<< rank<<"\tMPI slave " << rank << ": error in M-Shake algorithm!\n" << std::endl;
      }


      MPI::COMM_WORLD.Bcast(&next_step, 1, MPI::INT, 0);

      if (next_step == 2) {
        (*os) << "RANK: "<< rank<<"\tMessage from master: Steepest descent: minimum reached." << std::endl;
        error = 0;
        break;
      }
      else if (!next_step) {
        (*os) << "RANK: "<< rank<<"\tThere was an error in the master. Check output file for details."
              << "Exiting from MD main loop." << std::endl;
        error = 1;
        break;
      }

      ++sim.steps();
      sim.time() = sim.param().step.t0 + sim.steps() * sim.time_step_size();

    }

    if (error){
      (*os) << "RANK: "<< rank<<"\t\nErrors encountered during run in simulation "<< simulation_ID << " in Slave Thread "<< simulation_rank <<" - check above!\n" << std::endl;
      MPI_Abort(MPI_COMM_WORLD, E_INPUT_ERROR);
    }
    else{
      (*os) << "RANK: "<< rank<<"\n" GROMOSXX " MPI slave " << rank << "  simulation "<< simulation_ID << " in Slave Thread "<< simulation_rank << "finished successfully\n" << std::endl;
    }
}
    
