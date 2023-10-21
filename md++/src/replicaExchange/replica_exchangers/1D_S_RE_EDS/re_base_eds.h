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
 * File:   replica_exchange_base_eds.h
 * Author: bschroed
 *
 * Created on April 18, 2018, 3:38 PM
 * Modified June 18, 2021 - bschroed, srieder
 */
#include <replicaExchange/replica_exchangers/replica_exchange_base_interface.h>
#include <replicaExchange/replica_data.h>
#include <replicaExchange/repex_mpi.h>

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

namespace re
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
       std::map<ID_t, re::reeds_replica_stat_data > replicaStatData;
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

        int find_partner() const override;

        /*
        * energy calculation for statistical purposes of eds_stat() in replica_exchange_base_2d_l_T_HREMD.cc
        * for given configuration with new smoothing parameter s.
        */
        double calc_energy_eds_stat(double s);
        double calculate_energy(const unsigned int partner) override;
        double calculate_energy_core();
        void calculate_energy_helper(const unsigned int partnerReplicaID) override;
    };
}
#endif /* REPLICA_EXCHANGE_BASE_EDS_H */

