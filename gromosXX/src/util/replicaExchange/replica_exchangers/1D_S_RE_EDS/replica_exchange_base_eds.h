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
#include <util/replicaExchange/replica_exchangers/replica_exchange_base_interface.h>
#include <util/replicaExchange/replica_data.h>
#include <util/replicaExchange/repex_mpi.h>

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
     class replica_exchange_base_eds : public virtual replica_exchange_base_interface {
    public:
        replica_exchange_base_eds(io::Argument _args, 
                unsigned int cont, 
                unsigned int globalThreadID, 
                replica_graph_control & replicaGraphMPIControl,
                simulation::MpiControl & replica_mpi_control);
        /**
         * inits replica_reeds
         */
        void init() override;
        /** 
         * Initialize file and data for eds_stat output.
         */
        void init_eds_stat();
        
        
    protected:       
        virtual ~replica_exchange_base_eds();

       /*
       * Attributes
       */

        /**
         * for exchanging params easily
         */
        simulation::Parameter::reeds_struct& reedsParam;

       /**
       * stat information of all replicas for exchange optimisation 
       */
       std::map<ID_t, util::reeds_replica_stat_data > replicaStatData;
       /**
       *  output file stream for output file
       */
       std::map<ID_t, std::ofstream *> eds_stat_out;
       /**
        *  Exchange parameters
        */
       simulation::Parameter::eds_struct eds_para;
        /**
        * contains original forces which have to be reset after RE-EDS exchange energy calculation
        */
       math::VArray force_orig;
       /**
        * contains original virial which have to be reset after RE-EDS exchange energy calculation
        */
       math::Matrix virial_tensor_orig;
    
       
       /*
        * Functions
        */
       /**
        * Set RE - params
        */
       void setParams() override;

       void set_s();
       
        // RE-Exchange functions
        /**
        * Sets eds_struct() parameters to original value of replica
        */
        void reset_eds();
        /**
         * Sets  eds_struct() parameters to value of partner (for energy calculations)
         */
        void change_eds(const unsigned int partner);
        
        //Build exchange probabilities
        void determine_switch_probabilities() ;
        void swap_replicas_1D(const unsigned int partnerReplicaID);
        
        /*
        * energy calculation for statistical purposes of eds_stat() in replica_exchange_base.cc
        * for given configuration with new smoothing parameter s.
        */
        double calc_energy_eds_stat(double s);
        double calculate_energy_core();
        double calculate_energy(const unsigned int partner);
    };
}
#endif /* REPLICA_EXCHANGE_BASE_EDS_H */

