/* 
 * File:   replica_exchange_master.cc
 * Author: wissphil, sriniker
 * 
 * Created on April 29, 2011, 2:18 PM
 */

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
#include <util/replicaExchange/replica/replica.h>
#include <util/replicaExchange/replica_exchangers/replica_exchange_base_interface.h>
#include <util/replicaExchange/replica_exchangers/replica_exchange_master_interface.h>

#include <string>

#ifdef XXMPI
    #include <mpi.h>
#endif

#undef MODULE
#undef SUBMODULE
#define MODULE util
#define SUBMODULE replica_exchange

util::replica_exchange_master_interface::replica_exchange_master_interface(io::Argument & args,
        unsigned int cont,
        unsigned int globalThreadID,
        std::vector<std::vector<unsigned int> > replica_owned_threads,
        std::map<ID_t, rank_t> & thread_id_replica_map) :
        replica_exchange_base_interface(args, cont, globalThreadID, replica_owned_threads, thread_id_replica_map ),
        repParams(replica->sim.param().replica),
        repdatName(args["repdat"])
{
#ifdef XXMPI
  DEBUG(2,"replica_exchange_master "<< globalThreadID <<":Constructor:\t START");

  assert(globalThreadID == 0);    //TODO: This can be removed in future! bscrhoed
  assert(numReplicas > 0);
  DEBUG(2,"replica_exchange_master "<< globalThreadID <<":Constructor:\t rep_params THERE?");
  DEBUG(2,"replica_exchange_master "<< globalThreadID <<":Constructor:\t" << replica->sim.param().replica.num_l);
  DEBUG(2,"replica_exchange_master "<< globalThreadID <<":Constructor:\t" << replica->sim.param().replica.lambda[0]);

  assert(repParams.num_l > 0);
  
  DEBUG(4,"replica_exchange_master "<< globalThreadID <<":Constructor:\t Init Replicas \t Next");
  replicaData.resize(numReplicas);
  DEBUG(4,"replica_exchange_master "<< globalThreadID <<":Constructor:\t Replica_data type \t " << typeid(replicaData).name());

  //initialize data of replicas
  int ID = 0;
  for (int i = 0; i < repParams.num_l; ++i) {
    for (int j = 0; j < repParams.num_T; ++j) {
      replicaData[ID].ID = ID;
      replicaData[ID].T = repParams.temperature[j];
      DEBUG(4,"replica_exchange_master "<< globalThreadID <<":Constructor:\t Init Replicas \t "<< repParams.temperature[j]);
      replicaData[ID].l = repParams.lambda[i];
      replicaData[ID].dt = repParams.dt[i];
      ++ID;
    }
  }

  // set output file
 DEBUG(2,"replica_exchange_master "<< globalThreadID <<":Constructor:\t DONE");
#else
   throw "Cannot initialize replica_exchange_master without MPI!"; 
#endif
}


util::replica_exchange_master_interface::~replica_exchange_master_interface() {
   repOut.close();
}

void util::replica_exchange_master_interface::receive_from_all_slaves() {
    #ifdef XXMPI
    DEBUG(2,"replica_exchange_master "<< globalThreadID <<":receive_from_all_slaves:\t START\n");
    double start = MPI_Wtime();

    MPI_Status status;
    util::repInfo info;

    // receive all information from slaves
    for (unsigned int slaveReplicaID = 0; slaveReplicaID < numReplicas; ++slaveReplicaID) {
      if (slaveReplicaID != simulationID) {
        unsigned int replicaMasterThreadID = replica_owned_threads[slaveReplicaID][0];
        DEBUG(2,"replica_exchange_master "<< globalThreadID <<":receive_from_all_slaves:\t get_MPI replicaMasterThread: "<< replicaMasterThreadID << "\n");
        MPI_Recv(&info, 1, MPI_REPINFO, replicaMasterThreadID, REPINFO, MPI_COMM_WORLD, &status);
        replicaData[slaveReplicaID].run = info.run;
        replicaData[slaveReplicaID].epot = info.epot;
        replicaData[slaveReplicaID].epot_partner = info.epot_partner;
        replicaData[slaveReplicaID].probability = info.probability;
        replicaData[slaveReplicaID].switched = info.switched;
        replicaData[slaveReplicaID].partner = info.partner;
      }
      DEBUG(2,"replica_exchange_master "<< globalThreadID <<":receive_from_all_slaves:\t got all MPI reps\n");

      // write all information from master node to data structure
      replicaData[simulationID].run = run;
      replicaData[simulationID].partner = partnerReplicaID;
      replicaData[simulationID].epot = epot;
      replicaData[simulationID].epot_partner = epot_partner;
      replicaData[simulationID].probability = probability;
      replicaData[simulationID].switched = switched;
    }
    DEBUG(2,"replica_exchange_master "<< globalThreadID <<":receive_from_all_slaves:\t " << "time used for receiving all messages: " << MPI_Wtime() - start << " seconds\n");
    DEBUG(2,"replica_exchange_master "<< globalThreadID <<":receive_from_all_slaves:\t DONE: \n");
    #else
     throw "Cannot use replica_exchange_master without MPI!"; 
    #endif
}
