/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   replica_exchange_base_eds.h
 * Author: bschroed
 *
 * Created on April 18, 2018, 3:38 PM
 */
#include <util/replicaExchange/replica_exchangers/replica_exchange_base.h>
#include <util/replicaExchange/replica_data.h>
#include <util/replicaExchange/repex_mpi.h>
#include <util/replicaExchange/replica/replica_reeds.h>

#ifndef REPLICA_EXCHANGE_BASE_EDS_H
#define REPLICA_EXCHANGE_BASE_EDS_H

//const

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

#include <io/configuration/out_configuration.h>


#include <string>
#include <math/random.h>

#ifdef XXMPI
#include <mpi.h>
#endif

namespace util
{
     class replica_exchange_base_eds : public virtual replica_exchange_base {
    public:
        replica_exchange_base_eds(io::Argument _args, int cont, int rank, int simulationRank, int simulationID, int simulationThreads, 
                std::vector<int> repIDs, std::map<ID_t, rank_t>& _repMap);
         
        /**
         * @override
         * Tries a swapping of configuration if possible. Calculates energies, probabilities
         * and sends information via MPI communication if necessary. Added also reeds information
         */
        //using util::replica_exchange_base::swap;
        void swap() override;
        /**
        * runs MD simulation for all replicas; one by one
        */
        void run_MD() override;
        /**
         * inits replica_reeds
         */
        void init() override;
        
         /** 
         * calculate and write out acceptance probabilities for non-neighbour
         * switching attempts. Useful to optimize replica EDS parameter choice. 
         */
        void eds_stat();
          /** 
         * Initialize file and data for eds_stat output.
         */
        void init_eds_stat();
        
    protected:
        /**
        * all replicas on this node
        */
        typedef std::vector< util::replica_reeds* >::iterator repIterator;   //iterator for replicas
        std::vector<util::replica_reeds*> replicas;
        
        //void switch_coords_on_node(repIterator it1, const unsigned int partner) ;
        //void swap_on_node(repIterator it1, const unsigned int partner) ;
        virtual ~replica_exchange_base_eds();

       /**
       * stat information of all replicas for exchange optimisation 
       */
       std::map<ID_t, util::reeds_replica_stat_data > replicaStatData;
       /**
       *  output file stream for output file
       */
       std::map<ID_t, std::ofstream *> eds_stat_out;
       
        /**
        *  Other Functions:
        */
        //init Replicas - used in contstructor, initialises the replica objs.
        void createReplicas(int cont, std::vector<int>  repIDs, int rank) override;
    };
}
#endif /* REPLICA_EXCHANGE_BASE_EDS_H */

