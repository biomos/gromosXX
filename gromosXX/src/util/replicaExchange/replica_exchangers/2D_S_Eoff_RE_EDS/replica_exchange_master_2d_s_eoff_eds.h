/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * File:   replica_exchange_master_2d_s_eoff_eds.h
 * Author: theosm
 *
 * Created on March 29, 2020, 11:00 AM
 */

#include <util/replicaExchange/replica_exchangers/replica_exchange_base_interface.h>
#include <util/replicaExchange/replica_exchangers/replica_exchange_master_interface.h>
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


#ifndef REPLICA_EXCHANGE_MASTER_EDS_H
#define REPLICA_EXCHANGE_MASTER_EDS_H

namespace util{

    class replica_exchange_master_2d_s_eoff_eds :  public  replica_exchange_base_2d_s_eoff_eds, public  replica_exchange_master_interface  {
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
        replica_exchange_master_2d_s_eoff_eds(io::Argument _args,
                                    unsigned int cont,
                                    unsigned int globalThreeadID,
                                    replica_graph_mpi_control replicaGraphMPIControl,
                                    simulation::mpi_control_struct replica_mpi_control);

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
        ~replica_exchange_master_2d_s_eoff_eds(){};
        //using util::replica_exchange_base_2d_s_eoff_eds::replicas;
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
        std::vector<util::reeds_replica_data> replicaData;


        /**
         * functions, for initing the repout
         * @param repoutPath
         */
        //void init_repOut_stat_file(std::string repoutPath) override;

    };
}
#endif /* REPLICA_EXCHANGE_MASTER_EDS_H */
