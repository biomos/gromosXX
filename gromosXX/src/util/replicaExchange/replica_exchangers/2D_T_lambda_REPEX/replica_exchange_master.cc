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
#include "util/replicaExchange/replica_mpi_tools.h"
#include <util/replicaExchange/replica/replica.h>
#include <util/replicaExchange/replica_exchangers/replica_exchange_master_interface.h>
#include <util/replicaExchange/replica_exchangers/2D_T_lambda_REPEX/replica_exchange_master.h>

#include <string>

#ifdef XXMPI
    #include <mpi.h>
#endif

#undef MODULE
#undef SUBMODULE
#define MODULE util
#define SUBMODULE replica_exchange

util::replica_exchange_master::replica_exchange_master(io::Argument & args,
        unsigned int cont,
        unsigned int globalThreadID,
        replica_graph_mpi_control replicaGraphMPIControl,
        simulation::mpi_control_struct replica_mpi_control) :
        replica_exchange_base_interface(args, cont, globalThreadID, replicaGraphMPIControl, replica_mpi_control),
        replica_exchange_base(args, cont, globalThreadID, replicaGraphMPIControl, replica_mpi_control),
        replica_exchange_master_interface(args, cont, globalThreadID, replicaGraphMPIControl, replica_mpi_control)
{
#ifdef XXMPI
  DEBUG(2,"replica_exchange_master "<< globalThreadID <<":Constructor:\t START");

  assert(replicaGraphMPIControl.replicaGraphMasterID == replicaGraphMPIControl.replicaGraphThisThreadID);    //TODO: This can be removed in future! bscrhoed
  assert(replicaGraphMPIControl.numberOfReplicas > 0);
  DEBUG(2,"replica_exchange_master "<< globalThreadID <<":Constructor:\t rep_params THERE?");
  DEBUG(2,"replica_exchange_master "<< globalThreadID <<":Constructor:\t" << replica->sim.param().replica.num_l);
  DEBUG(2,"replica_exchange_master "<< globalThreadID <<":Constructor:\t" << replica->sim.param().replica.lambda[0]);

  assert(repParams.num_l > 0);
  
  DEBUG(4,"replica_exchange_master "<< globalThreadID <<":Constructor:\t Init Replicas \t Next");
  replicaData.resize(replicaGraphMPIControl.numberOfReplicas);
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


util::replica_exchange_master::~replica_exchange_master() {
   repOut.close();
}

