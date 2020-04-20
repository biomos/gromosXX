/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   replica_exchange_slave_eds.cc
 * Author: bschroed
 * 
 * Created on August 31, 2018, 10:43 AM
 */
#include "util/replicaExchange/replica_mpi_tools.h"
#include <util/replicaExchange/replica_exchangers/1D_S_RE_EDS/replica_exchange_slave_eds.h>


#undef MODULE
#undef SUBMODULE
#define MODULE util
#define SUBMODULE replica_exchange

util::replica_exchange_slave_eds::replica_exchange_slave_eds(io::Argument & _args,
                                                            unsigned int cont,
                                                            unsigned int globalThreadID,
                                                            replica_graph_mpi_control replicaGraphMPIControl,
                                                            simulation::mpi_control_struct replica_mpi_control):
            replica_exchange_base_interface(_args, cont, globalThreadID, replicaGraphMPIControl, replica_mpi_control),
            replica_exchange_base_eds(_args, cont, globalThreadID, replicaGraphMPIControl, replica_mpi_control),
            replica_exchange_slave_interface(_args, cont, globalThreadID, replicaGraphMPIControl, replica_mpi_control)
{

    //MPI_DEBUG(2,"replica_exchange_slave_eds " << globalThreadID << ":Constructor:\t START");

    //MPI_DEBUG(2,"replica_exchange_slave_eds " << globalThreadID << ":Constructor:\t DONE"); 
}

void util::replica_exchange_slave_eds::send_to_master() const{
  #ifdef XXMPI
  DEBUG(2,"replica_exchange_slave_eds " << globalThreadID << ":send_to_master:\t START");
  if(not_sender){
  DEBUG(2,"replica_exchange_slave_eds " << globalThreadID << ":send_to_master:\t Show vPots");

    util::repInfo info;
    std::vector<double> eds_energies;
    info.run = run;
    info.epot = epot;
    info.epot_partner = epot_partner;
    info.partner = partnerReplicaID;
    info.probability = probability;
    info.switched = int(switched);
    
    DEBUG(4,"replica_exchange_slave_eds " << globalThreadID << "send_to_master:\t epotTot\t "<< info.epot);

    DEBUG(4,"replica_exchange_slave_eds " << globalThreadID << ":send_to_master:\t\t send MPI_REPINFO");
    MPI_Send(&info, 1, MPI_REPINFO, replicaGraphMPIControl.masterID, REPINFO, replicaGraphMPIControl.comm);

    DEBUG(4,"replica_exchange_slave_eds " << globalThreadID << ":send_to_master:\t\t send MPI_EDS");

    eds_energies= replica->conf.current().energies.eds_vi;

    MPI_Send(&eds_energies[0], 1, MPI_EDSINFO, replicaGraphMPIControl.masterID, EDSINFO, replicaGraphMPIControl.comm);
    DEBUG(4,"replica_exchange_slave_eds " << globalThreadID << ":send_to_master:\t\t send MPI_EDS \t DONE" );
  }
  DEBUG(2,"replica_exchange_slave_eds " << globalThreadID << ":send_to_master:\t DONE");

  #else
    throw "Cannot use send_to_master from replica_exchange_slave_eds without MPI!"; 
  #endif
}

void util::replica_exchange_slave_eds::swap() { //todo move into interface!?
    MPI_DEBUG(3,"replica_exchange_slave_eds "<< globalThreadID <<":swap:\t START");
    if(not_sender){
            partnerReplicaID = find_partner();
    
        if (partnerReplicaID != simulationID) // different replica?
        {
          //TODO: RENAME 
          swap_replicas_2D(partnerReplicaID);
          if (switched) {
            if (globalThreadID < partnerReplicaID) {
              send_coord(partnerReplicaID);
              receive_new_coord(partnerReplicaID);
              // the averages of current and old are interchanged after calling exchange_state() and have to be switched back
              exchange_averages();
            } else {
              receive_new_coord(partnerReplicaID);
              send_coord(partnerReplicaID);
              // the averages of current and old are interchanged after calling exchange_state() and have to be switched back
              exchange_averages();
            }
          }
        }
        else {  // no exchange with replica itself
          probability = 0.0;
          switched = 0;
        }
        if(switched && replica->sim.param().replica.scale) {
          velscale(partnerReplicaID);
        }

    }   else{
        DEBUG(3, "replica_exchange_slave_eds "<< globalThreadID <<":swap:\tNOT sending")
    }
    MPI_DEBUG(3,"replica_exchange_slave_eds "<< globalThreadID <<":swap:\t DONE");

}



