/* 
 * File:   replica_exchange_base.cc
 * Author: wissphil, sriniker
 * 
 * Created on April 29, 2011, 2:06 PM
 */

#include "replicaExchange/replica/replica.h"
#include "replicaExchange/replica/replica_MPI_master.h"
#include "replicaExchange/replica/replica_MPI_slave.h"
#include "replicaExchange/replica_exchangers/2D_T_lambda_REPEX/replica_exchange_base.h"

#include <stdheader.h>

#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <simulation/simulation.h>
#include <configuration/configuration.h>

#include <algorithm/algorithm/algorithm_sequence.h>
#include <interaction/interaction.h>
#include <interaction/forcefield/forcefield.h>

#include <io/argument.h>
#include <util/usage.h>
#include <util/error.h>

#include <io/read_input.h>
#include <io/print_block.h>

#include <time.h>
#include <unistd.h>

#include <io/configuration/out_configuration.h>
#include <math/random.h>
#include <math/volume.h>
#include <string>

#ifdef XXMPI
#include <mpi.h>
#endif

#undef MODULE
#undef SUBMODULE
#define MODULE util
#define SUBMODULE replica_exchange

re::replica_exchange_base::replica_exchange_base(io::Argument _args,
                                                   unsigned int cont, 
                                                   unsigned int globalThreadID,
                                                   replica_graph_control & replicaGraphMPIControl,
                                                   simulation::MpiControl & replica_mpi_control) : 
        replica_exchange_base_interface(_args, cont, globalThreadID, replicaGraphMPIControl, replica_mpi_control){
#ifdef XXMPI
  DEBUG(3,"replica_exchange_base "<< globalThreadID <<":Constructor:\t START ");
  DEBUG(3,"replica_exchange_base "<< globalThreadID <<":Constructor:\t SIMULATIONID:  "<< simulationID);

  //construct replica obj
  DEBUG(4,"replica_exchange_base_interface "<< globalThreadID <<":Constructor:\t  createReplica START");
  createReplicas(cont, globalThreadID, replica_mpi_control);
  DEBUG(4,"replica_exchange_base_interface "<< globalThreadID <<":Constructor:\t createdReplica DONE" << replica);

  //RE-Vars
  DEBUG(3,"replica_exchange_base_eds "<< globalThreadID <<":Constructor:\t setParams" );
  setParams();  
  
  DEBUG(3,"replica_exchange_base "<< globalThreadID <<":Constructor:\t Constructor \t DONE");
  #else
    throw "Cannot use send_to_master from replica_exchange_slave_eds without MPI!"; 
  #endif
}


re::replica_exchange_base::~replica_exchange_base() {
    delete replica;
}

void re::replica_exchange_base::setParams(){
    // set some variables
    stepsPerRun = replica->sim.param().step.number_of_steps;
    DEBUG(4,"replica_exchange_base_interface "<< globalThreadID <<":setParams:\t NUMBER OF STEPS "<<stepsPerRun);

    run = 0;
    total_runs = replica->sim.param().replica.trials + replica->sim.param().replica.equilibrate;
    DEBUG(4,"replica_exchange_base_interface "<< globalThreadID <<":setParams:\t NUMBER OF total_runs "<<total_runs);

    partnerReplicaID = simulationID;
    time = replica->sim.time();
    steps = 0;
    switched = 0;
    replica->curentStepNumber=0;
    replica->totalStepNumber = total_runs*stepsPerRun;
    replica->stepsPerRun= stepsPerRun;

    const int numT = replica->sim.param().replica.num_T;
    DEBUG(4,"replica_exchange_base_interface "<< globalThreadID <<":setParams:\t PARAM START");

    T = replica->sim.param().replica.temperature[simulationID % numT];
    DEBUG(4,"replica_exchange_base_interface "<< globalThreadID <<":setParams:\t T START");

    l = replica->sim.param().replica.lambda[simulationID / numT];
    DEBUG(4,"replica_exchange_base_interface "<< globalThreadID <<":setParams:\t Lamba");

    dt = replica->sim.param().replica.dt[simulationID / numT];
    DEBUG(4,"replica_exchange_base_interface "<< globalThreadID <<":setParams:\t PARAM DONE ");

    set_lambda();
    set_temp();
}

//SWAPPING Functions
void re::replica_exchange_base::determine_switch_probabilities(){
    DEBUG(3,"replica_exchange_base "<< globalThreadID <<":determineSwitchPos:\t START");
    swap_replicas_2D(partnerReplicaID);
    DEBUG(3,"replica_exchange_base "<< globalThreadID <<":determineSwitchPos:\t DONE");
}

void re::replica_exchange_base::swap_replicas_2D(const unsigned int partnerReplicaID) {
  DEBUG(4, "replica "<<  globalThreadID <<":swap2D:\t  START");
  DEBUG(4, "replica "<<  globalThreadID <<":swap2D:\t  sim: "<< simulationID << " \t\t "<< partnerReplicaID);

  DEBUG(1,"replica "<<  globalThreadID <<":replica_exchange_base_interface: SWAP_2D\n\n");

  unsigned int partnerReplicaMasterThreadID = partnerReplicaID;
  unsigned int numT = replica->sim.param().replica.num_T;
  unsigned int numL = replica->sim.param().replica.num_l;

  // does partner exist?
  if (partnerReplicaID < numT * numL && partnerReplicaID != simulationID) {
    // the one with lower ID does probability calculation
    if (simulationID < partnerReplicaID) {
        DEBUG(4, "replica "<<  globalThreadID <<":swap2D:\t  swap because smaller :)");

      // posts a MPI_Recv(...) matching the MPI_Send below
      probability = calc_probability(partnerReplicaID);
      const double randNum = rng.get();

      std::vector<double> prob(2);
      prob[0] = probability;
      prob[1] = randNum;
#ifdef XXMPI
      MPI_Send(&prob[0], 2, MPI_DOUBLE, partnerReplicaID, SENDCOORDS, replicaGraphMPIControl().comm);
#endif

      if (randNum < probability) {
        switched = true;
      } else
        switched = false;
    } else {    //The Partner sends his data to The calculating Thread
      //special case if lambda also needs to be exchanged
      bool sameLambda = (l == replica->sim.param().replica.lambda[partnerReplicaID / replica->sim.param().replica.num_T]);
      if(!sameLambda){      //exchange LAMBDA
        // E21: Energy with configuration 2 and lambda 1(of partner)
        const double E21 = calculate_energy(partnerReplicaID);
        // this we can store as the partner energy of the current replica
        epot_partner = E21;
        // E22: Energy with configuration 2 and lambda 2(own one)
#ifdef XXMPI
        const double E22 = epot;
        // send E21 and E22
        double energies[2] = {E22, E21};
        //this send operation is matched in calc_probability()
        DEBUG(1,"\n\nreplica "<<  globalThreadID <<"replica_exchange_base_interface: SWAP_2D before Send\n");
        MPI_Send(&energies[0], 2, MPI_DOUBLE, partnerReplicaID, SWITCHENERGIES,  replicaGraphMPIControl().comm);
        DEBUG(1,"\n\nreplica "<<  globalThreadID <<"replica_exchange_base_interface: SWAP_2D after Send\n");
#endif
      } else { // sameLambda
#ifdef XXMPI
        double energies[2] = {epot, 0.0};
        MPI_Send(&energies[0],2,MPI_DOUBLE, partnerReplicaID, SWITCHENERGIES,  replicaGraphMPIControl().comm);
#endif
     }
      if (replica->sim.param().pcouple.scale != math::pcouple_off) {
#ifdef XXMPI
        math::Box box_replica = replica->conf.current().box;    //exchange box
        MPI_Send(&box_replica(0)[0], 1, MPI_BOX, partnerReplicaID, BOX,  replicaGraphMPIControl().comm);
#endif
      }

#ifdef XXMPI
      MPI_Status status;
#endif
      std::vector<double> prob;
      prob.resize(2);
#ifdef XXMPI
      MPI_Recv(&prob[0], 2, MPI_DOUBLE, partnerReplicaID, SENDCOORDS,  replicaGraphMPIControl().comm, &status);
#endif
      //Have we been exchanged little partner?
      probability = prob[0];
      double randNum = prob[1];

      if (randNum < probability) {
        switched = true;
      } else {
        switched = false;
      }
    }

  } else {//This should be an error!
    //throw "Partner does not exist!";
    DEBUG(1, "replica "<<  globalThreadID <<":swap2D:\t  No swap because edgy");
    switched = false;
    probability = 0.0;
  }
    DEBUG(1, "replica "<< globalThreadID <<":swap2D:\t  DONE");
}


