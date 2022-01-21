/* 
 * File:   replica_exchange_slave_eds.h
 * Author: bschroed
 *
 * Created on August 31, 2018, 10:43 AM
 * Modified June 18, 2021 - bschroed, srieder
 */

#include <replicaExchange/replica_exchangers/replica_exchange_slave_interface.h>
#include <replicaExchange/replica_exchangers/1D_S_RE_EDS/replica_exchange_base_eds.h>
//for the constructor
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

#ifndef REPLICA_EXCHANGE_SLAVE_EDS_H
#define REPLICA_EXCHANGE_SLAVE_EDS_H

namespace re{
    class replica_exchange_slave_eds: public replica_exchange_base_eds, public replica_exchange_slave_interface {
    public:
        replica_exchange_slave_eds(io::Argument & _args,
                unsigned int cont,
                unsigned int globalThreadID,
                replica_graph_control & replicaGraphMPIControl,
                simulation::MpiControl & replica_mpi_control);
        /**
        * sends information of all replicas to master
        */
        void send_to_master() const override;

    private:
        //replica_exchange_slave_eds(const replica_exchange_slave_eds& orig);
        virtual ~replica_exchange_slave_eds(){};

        //give all information of this node to Master.
        replica_exchange_slave_eds(const replica_exchange_slave_eds& orig); //Todo: Messy method, bschroed
        
        
            };
}
#endif /* REPLICA_EXCHANGE_SLAVE_EDS_H */

