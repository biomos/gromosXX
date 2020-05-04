/* 
 * File:   replica_exchange_master.cc
 * Author: wissphil, sriniker
 * 
 * Created on April 29, 2011, 2:06 PM
 */

#define REPEX_MPI
#include <util/replicaExchange/repex_mpi.h>
#include "util/replicaExchange/replica_mpi_tools.h"
#include <util/replicaExchange/replica_exchangers/replica_exchange_base_interface.h>
#include <util/replicaExchange/replica_exchangers/replica_exchange_slave_interface.h>

#include <stdheader.h>

#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <simulation/simulation.h>
#include <configuration/configuration.h>

#include <algorithm/algorithm/algorithm_sequence.h>
#include <interaction/interaction.h>
#include <interaction/forcefield/forcefield.h>

#include <io/argument.h>
#include <time.h>
#include <unistd.h>

#include <io/configuration/out_configuration.h>


#ifdef XXMPI
#include <mpi.h>
#endif

util::replica_exchange_slave_interface::replica_exchange_slave_interface(io::Argument & _args,
                                                    unsigned int cont,
                                                    unsigned int globalThreadID,
                                                    replica_graph_mpi_control replicaGraphMPIControl,
                                                    simulation::MpiControl replica_mpi_control) : 
        replica_exchange_base_interface(_args, cont, globalThreadID, replicaGraphMPIControl, replica_mpi_control) {
        MPI_DEBUG(2, "replica_exchange_slave_interface "<< globalThreadID <<":Constructor:\t START");
        
        if(replica_mpi_control.masterID == replica_mpi_control.threadID){
            DEBUG(5, "replica_exchange_slave_interface "<< globalThreadID <<":Constructor:\t RE-Graph Sender");
            not_sender = true;
        }
        else{
            DEBUG(5, "replica_exchange_slave_interface "<< globalThreadID <<":Constructor:\t I'm not an RE-Graph Sender");
            not_sender = false;
        }
        
        MPI_DEBUG(2, "replica_exchange_slave_interface "<< globalThreadID <<":Constructor:\t DONE");
}

util::replica_exchange_slave_interface::~replica_exchange_slave_interface() {
}

void util::replica_exchange_slave_interface::send_to_master() const {
#ifdef XXMPI
  
  if(not_sender){
    DEBUG(2,"replica_exchange_slave_interface " << globalThreadID << ":send_to_master \t START");
      util::repInfo info;
      info.run = run;
      info.epot = epot;
      info.epot_partner = epot_partner;
      info.partner = partnerReplicaID;
      info.probability = probability;
      info.switched = int(switched);
      MPI_Send(&info, 1, MPI_REPINFO, replicaGraphMPIControl.masterID, REPINFO, replicaGraphMPIControl.comm);

    DEBUG(2,"replica_exchange_slave_interface " << globalThreadID << ":\t send_to_master \t Done");
  }
#endif
}
