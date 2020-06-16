/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * File:   replica_exchange_slave_2d_s_eoff_eds.h
 * Author: theosm
 *
 * Created on March 29, 2020, 11:03 AM
 */

#include <util/replicaExchange/replica_exchangers/replica_exchange_slave_interface.h>
#include <util/replicaExchange/replica_exchangers/2D_S_Eoff_RE_EDS/replica_exchange_base_2d_s_eoff_eds.h>
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

#ifndef REPLICA_EXCHANGE_SLAVE_2D_S_EOFF_EDS_H
#define REPLICA_EXCHANGE_SLAVE_2D_S_EOFF_EDS_H

namespace util{
    class replica_exchange_slave_2d_s_eoff_eds: public replica_exchange_base_2d_s_eoff_eds, public replica_exchange_slave_interface {
    public:
        replica_exchange_slave_2d_s_eoff_eds(io::Argument & _args,
                unsigned int cont,
                unsigned int globalThreadID,
                replica_graph_control &replicaGraphMPIControl,
                simulation::MpiControl &replica_mpi_control);
        /**
        * sends information of all replicas to master
        */
        void send_to_master() const override;

    private:
        //replica_exchange_slave_2d_s_eoff_eds(const replica_exchange_slave_2d_s_eoff_eds& orig);
        virtual ~replica_exchange_slave_2d_s_eoff_eds(){};
    };
}
#endif /* REPLICA_EXCHANGE_SLAVE_2D_S_EOFF_EDS_H */
