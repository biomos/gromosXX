/* 
 * File:   replica_exchange_master.cc
 * Author: wissphil, sriniker
 * 
 * Created on April 29, 2011, 2:06 PM
 */

#define REPEX_MPI
#include <util/replicaExchange/replica_exchangers/replica_exchange_slave.h>

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

#include <util/replicaExchange/repex_mpi.h>
#include <util/replicaExchange/replica_exchangers/replica_exchange_base.h>

#ifdef XXMPI
#include <mpi.h>
#endif

util::replica_exchange_slave::replica_exchange_slave(io::Argument & _args,
                                                    unsigned int cont,
                                                    unsigned int globalThreadID,
                                                    std::vector<std::vector<unsigned int> > replica_owned_threads,
                                                    std::map<ID_t, rank_t> & thread_id_replica_map) : 
    replica_exchange_base(_args, cont, globalThreadID, replica_owned_threads, thread_id_replica_map) {
        DEBUG(4, "replica_exchange_slave "<< globalThreadID <<":Constructor:\t START");

        DEBUG(4, "replica_exchange_slave "<< globalThreadID <<":Constructor:\t DONE");
}

util::replica_exchange_slave::~replica_exchange_slave() {
}

void util::replica_exchange_slave::send_to_master() const {
#ifdef XXMPI
  DEBUG(2,"replica_exchange_slave " << globalThreadID << ":send_to_master \t START");
    util::repInfo info;
    info.run = run;
    info.epot = epot;
    info.epot_partner = epot_partner;
    info.partner = partnerReplicaID;
    info.probability = probability;
    info.switched = int(switched);
    MPI_Send(&info, 1, MPI_REPINFO, globalMasterThreadID, REPINFO, MPI_COMM_WORLD);

  DEBUG(2,"replica_exchange_slave " << globalThreadID << ":\t send_to_master \t Done");
#endif
}
