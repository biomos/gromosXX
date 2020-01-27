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

#include <util/replicaExchange/replica_exchangers/1D_S_RE_EDS/replica_exchange_slave_eds.h>


#undef MODULE
#undef SUBMODULE
#define MODULE util
#define SUBMODULE replica_exchange

util::replica_exchange_slave_eds::replica_exchange_slave_eds(io::Argument & _args,
                                                            unsigned int cont,
                                                            unsigned int globalThreadID,
                                                            std::vector<std::vector<unsigned int> > replica_owned_threads,
                                                            std::map<ID_t, rank_t> & thread_id_replica_map,
                                                            simulation::mpi_control_struct replica_mpi_control):
            replica_exchange_base(_args, cont, globalThreadID, replica_owned_threads, thread_id_replica_map, replica_mpi_control),
            replica_exchange_base_eds(_args, cont, globalThreadID, replica_owned_threads, thread_id_replica_map, replica_mpi_control),
            replica_exchange_slave(_args, cont, globalThreadID, replica_owned_threads, thread_id_replica_map, replica_mpi_control)
{

    DEBUG(2,"replica_exchange_slave_eds " << globalThreadID << ":Constructor:\t START");

    DEBUG(2,"replica_exchange_slave_eds " << globalThreadID << ":Constructor:\t DONE"); 
}

void util::replica_exchange_slave_eds::send_to_master() const{
  #ifdef XXMPI
  DEBUG(2,"replica_exchange_slave_eds " << globalThreadID << ":send_to_master:\t START");
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
    MPI_Send(&info, 1, MPI_REPINFO, globalMasterThreadID, REPINFO, MPI_COMM_WORLD);

    DEBUG(4,"replica_exchange_slave_eds " << globalThreadID << ":send_to_master:\t\t send MPI_EDS");

    eds_energies= replica->conf.current().energies.eds_vi;

    MPI_Send(&eds_energies[0], 1, MPI_EDSINFO, globalMasterThreadID, EDSINFO, MPI_COMM_WORLD);
    DEBUG(4,"replica_exchange_slave_eds " << globalThreadID << ":send_to_master:\t\t send MPI_EDS \t DONE" );

  DEBUG(2,"replica_exchange_slave_eds " << globalThreadID << ":send_to_master:\t DONE");
  #else
    throw "Cannot use send_to_master from replica_exchange_slave_eds without MPI!"; 
  #endif
}

