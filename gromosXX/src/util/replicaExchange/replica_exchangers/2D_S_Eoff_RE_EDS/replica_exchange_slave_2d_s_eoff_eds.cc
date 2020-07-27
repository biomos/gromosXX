/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * File:   replica_exchange_slave_2d_s_eoff_eds.cc
 * Author: theosm
 *
 * Created on March 29, 2020, 11:05 AM
 */
#include "util/replicaExchange/replica_mpi_tools.h"
#include <util/replicaExchange/replica_exchangers/2D_S_Eoff_RE_EDS/replica_exchange_slave_2d_s_eoff_eds.h>


#undef MODULE
#undef SUBMODULE
#define MODULE util
#define SUBMODULE replica_exchange

util::replica_exchange_slave_2d_s_eoff_eds::replica_exchange_slave_2d_s_eoff_eds(io::Argument & _args,
                                                            unsigned int cont,
                                                            unsigned int globalThreadID,
                                                            replica_graph_control &replicaGraphMPIControl,
                                                            simulation::MpiControl &replica_mpi_control):
            replica_exchange_base_interface(_args, cont, globalThreadID, replicaGraphMPIControl, replica_mpi_control),
            replica_exchange_base_2d_s_eoff_eds(_args, cont, globalThreadID, replicaGraphMPIControl, replica_mpi_control),
            replica_exchange_slave_interface(_args, cont, globalThreadID, replicaGraphMPIControl, replica_mpi_control)
{

    MPI_DEBUG(2,"replica_exchange_slave_2d_s_eoff_eds " << globalThreadID << ":Constructor:\t START");

    MPI_DEBUG(2,"replica_exchange_slave_2d_s_eoff_eds " << globalThreadID << ":Constructor:\t DONE");
}

void util::replica_exchange_slave_2d_s_eoff_eds::send_to_master() const{
  #ifdef XXMPI
  if(not_sender){
  DEBUG(2,"replica_exchange_slave_2d_s_eoff_eds " << globalThreadID << ":send_to_master:\t START");
  DEBUG(2,"replica_exchange_slave_2d_s_eoff_eds " << globalThreadID << ":send_to_master:\t Show vPots");

    util::repInfo info;
    std::vector<double> eds_energies;
    info.run = run;
    info.epot = epot;
    info.epot_partner = epot_partner;
    info.partner = partnerReplicaID;
    info.probability = probability;
    info.switched = switched;

    DEBUG(4,"replica_exchange_slave_2d_s_eoff_eds " << globalThreadID << "send_to_master:\t epotTot\t "<< info.epot);

    DEBUG(4,"replica_exchange_slave_2d_s_eoff_eds " << globalThreadID << ":send_to_master:\t\t send MPI_REPINFO");
    MPI_Send(&info, 1, MPI_REPINFO, replicaGraphMPIControl().masterID, REPINFO, replicaGraphMPIControl().comm);

    DEBUG(4,"replica_exchange_slave_2d_s_eoff_eds " << globalThreadID << ":send_to_master:\t\t send MPI_EDS");

    eds_energies= replica->conf.current().energies.eds_vi;

    MPI_Send(&eds_energies[0], 1, MPI_EDSINFO, replicaGraphMPIControl().masterID, EDSINFO, replicaGraphMPIControl().comm);
    DEBUG(4,"replica_exchange_slave_2d_s_eoff_eds " << globalThreadID << ":send_to_master:\t\t send MPI_EDS \t DONE" );
    DEBUG(2,"replica_exchange_slave_2d_s_eoff_eds " << globalThreadID << ":send_to_master:\t DONE");
  }

  #else
    throw "Cannot use send_to_master from replica_exchange_slave_2d_s_eoff_eds without MPI!";
  #endif
}