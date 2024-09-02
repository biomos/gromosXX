/*
 * This file is part of GROMOS.
 * 
 * Copyright (c) 2011, 2012, 2016, 2018, 2021, 2023 Biomos b.v.
 * See <https://www.gromos.net> for details.
 * 
 * GROMOS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */


/*
 * File:   re_master_eds.h
 * Author: bschroed
 *
 * Created on April 18, 2018, 3:20 PM
 * Modified June 18, 2021 - bschroed, srieder
 */

#include <replicaExchange/replica_exchangers/replica_exchange_base_interface.h>
#include <replicaExchange/replica_exchangers/replica_exchange_master_interface.h>
#include <replicaExchange/replica_exchangers/1D_S_RE_EDS/re_base_eds.h>

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


#ifndef REPLICA_EXCHANGE_MASTER_EDS_H
#define REPLICA_EXCHANGE_MASTER_EDS_H

namespace re{

    class replica_exchange_master_eds :  public  replica_exchange_base_eds, public  replica_exchange_master_interface  {
    public:
        /**
         * constructor
         * @param args io::Argument, passed on to replica
         * @param rank integer, rank of node
         * @param _size size of mpi comm world
         * @param _numReplicas total number of replicas
         * @param repIDs std::vector<int>, IDs of replicas the instance has to manage
         * @param repMap std::map<int,int>, maps replica IDs to nodes; needed for communication
         */
        replica_exchange_master_eds(io::Argument _args,
                                    unsigned int cont,
                                    unsigned int globalThreeadID,
                                    replica_graph_control & replicaGraphMPIControl,
                                    simulation::MpiControl & replica_mpi_control);

        /**
         * @override
         * init MD simulation for all replicas; one by one
         */
        //void init();
        void receive_from_all_slaves() override;

        void write() override;

        /**
         * functions, for initing the repout
         * @param repoutPath
         */
        void init_repOut_stat_file() override;
    protected:
        /**
         * destructor
         */
        ~replica_exchange_master_eds(){};
        //using re::replica_exchange_base_eds::replicas;
        /**
         * Column Size for redpat out-floating point nums
         */
        int svalPrecision= 5;
        int potEPrecision = 2;
        int generalPrecision = 2;

         /**
         * determines the digits needed to represent the smalles S-value in eds sim
         */
        int getSValPrecision(double minLambda);
        /*
        * global Parameters for replica exchange simulation
        * int num_T;
        * int num_l;
        * std::vector<double> temperature;
        * bool scale;
        * std::vector<double> lambda;
        * std::vector<double> dt;
        * int trials;
        * int equilibrate;
        * int slave_runs;
        * int write;
        */
        /**
        * information of all replicas
        */
        std::vector<re::reeds_replica_data> replicaData;


        /**
         * functions, for initing the repout
         * @param repoutPath
         */
        //void init_repOut_stat_file(std::string repoutPath) override;

    };
}
#endif /* REPLICA_EXCHANGE_MASTER_EDS_H */
