/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   replica_exchange_slave_eds.h
 * Author: bschroed
 *
 * Created on August 31, 2018, 10:43 AM
 */

#include <util/replicaExchange/replica_exchangers/replica_exchange_slave_interface.h>
#include <util/replicaExchange/replica_exchangers/1D_S_RE_EDS/replica_exchange_base_eds.h>
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

namespace util{
    class replica_exchange_slave_eds: public replica_exchange_base_eds, public replica_exchange_slave_interface {
    public:
        replica_exchange_slave_eds(io::Argument & _args,
                unsigned int cont,
                unsigned int globalThreadID,
                replica_graph_mpi_control replicaGraphMPIControl,
                simulation::mpi_control_struct replica_mpi_control);
        /**
        * sends information of all replicas to master
        */
        void send_to_master() const override;

        void swap() override;
        
    private:
        //replica_exchange_slave_eds(const replica_exchange_slave_eds& orig);
        virtual ~replica_exchange_slave_eds(){};

        //give all information of this node to Master.
        replica_exchange_slave_eds(const replica_exchange_slave_eds& orig); //Todo: Messy method, bschroed
            };
}
#endif /* REPLICA_EXCHANGE_SLAVE_EDS_H */

