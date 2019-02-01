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

#include <util/replicaExchange/replica_exchange_slave_eds.h>

#undef MODULE
#undef SUBMODULE
#define MODULE util
#define SUBMODULE replica_exchange


util::replica_exchange_slave_eds::replica_exchange_slave_eds(io::Argument & _args,
            int cont,
            int rank,
            std::vector<int> repIDs,
            std::map<ID_t, rank_t> & repMap): 
            replica_exchange_base(_args, cont, rank, repIDs, repMap),
            replica_exchange_slave(_args, cont, rank, repIDs, repMap),
            replica_exchange_base_eds(_args, cont, rank, repIDs, repMap),
            reedsParam(replicas[0]->sim.param().reeds){

    DEBUG(2,"replica_exchange_slave_eds " << rank << ":\t Constructor \t START");
   
    DEBUG(2,"replica_exchange_slave_eds " << rank << ":\t Constructor \t DONE"); 
}

void util::replica_exchange_slave_eds::send_to_master() const{
  DEBUG(2,"replica_exchange_slave_eds " << rank << ":\t send_to_master \t START");
    DEBUG(2,"replica_exchange_slave_eds " << rank << ":\t send_to_master \t Show vPots");
    for(int s=0;s<reedsParam.eds_para[0].numstates; s++){
        DEBUG(2,"replica_exchange_slave_eds " << rank << ":\t send_to_master \t State "<< s << "En: "<< std::min(replicas[0]->conf.current().energies.eds_vi[s], 10000000.0));
        std::cerr << "EEE: " << std::min(replicas[0]->conf.current().energies.eds_vi[s], 10000000.0);//Output potential energies for each state 
    }
  for (std::vector<util::replica_reeds *>::const_iterator it = replicas.begin(); it < replicas.end(); ++it) {
    util::repInfo info;
    std::vector<double> eds_energies;
    info.run = (*it)->run;
    info.epot = (*it)->epot;
    info.epot_partner = (*it)->epot_partner;
    info.partner = (*it)->partner;
    info.probability = (*it)->probability;
    info.switched = int((*it)->switched);
    DEBUG(2,"replica_exchange_slave_eds " << rank << ":\t send_to_master \t\t send MPI_REPINFO");
    MPI_Send(&info, 1, MPI_REPINFO, 0, REPINFO, MPI_COMM_WORLD);

    DEBUG(2,"replica_exchange_slave_eds " << rank << ":\t send_to_master \t\t send MPI_EDS");
    
    eds_energies=(*it)->conf.current().energies.eds_vi;
    DEBUG(2,"replica_exchange_slave_eds " << rank << ":\t send_to_master \t\t send MPI_EDS \t current energy: " << eds_energies[0]);
    
    MPI_Send(&eds_energies, 1, MPI_EDSINFO, 0, EDSINFO, MPI_COMM_WORLD);
    DEBUG(2,"replica_exchange_slave_eds " << rank << ":\t send_to_master \t\t send MPI_EDS \t recived energy: " << eds_energies[0]);
    DEBUG(2,"replica_exchange_slave_eds " << rank << ":\t send_to_master \t\t send MPI_EDS \t DONE" );
    
  }
  DEBUG(2,"replica_exchange_slave_eds " << rank << ":\t send_to_master \t DONE");
}

