/* 
 * File:   hamiltonian_simulatedAnnealing_slave_eds.h
 * Author: bschroed
 *
 * Created on August 31, 2018, 10:43 AM
 * Modified June 18, 2021 - bschroed, srieder
 */

#include <util/replicaExchange/replica_exchangers/replica_exchange_slave_interface.h>
#include <util/replicaExchange/replica_exchangers/1D_S_HSA_EDS/hamiltonian_simulatedAnnealing_base_eds.h>
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

#ifndef hamiltonian_simulatedAnnealing_slave_eds_H
#define hamiltonian_simulatedAnnealing_slave_eds_H

namespace re{
    class hamiltonian_simulatedAnnealing_slave_eds: public hamiltonian_simulatedAnnealing_base_eds, public replica_exchange_slave_interface {
    public:
        hamiltonian_simulatedAnnealing_slave_eds(io::Argument & _args,
                unsigned int cont,
                unsigned int globalThreadID,
                replica_graph_control & replicaGraphMPIControl,
                simulation::MpiControl & replica_mpi_control);
        /**
        * sends information of all replicas to master
        */
        void send_to_master() const override;

    private:
        //hamiltonian_simulatedAnnealing_slave_eds(const hamiltonian_simulatedAnnealing_slave_eds& orig);
        virtual ~hamiltonian_simulatedAnnealing_slave_eds(){};

        //give all information of this node to Master.
        hamiltonian_simulatedAnnealing_slave_eds(const hamiltonian_simulatedAnnealing_slave_eds& orig); //Todo: Messy method, bschroed
        
        
            };
}
#endif /* hamiltonian_simulatedAnnealing_slave_eds_H */

