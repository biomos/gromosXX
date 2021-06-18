/*
 * File:   replica_exchange_base_2d_s_eoff_eds.h
 * Author: theosm
 *
 * Created on March 29, 2020, 11:00 AM
 */
#include <util/replicaExchange/replica_exchangers/replica_exchange_base_interface.h>
#include <util/replicaExchange/replica_data.h>
#include <util/replicaExchange/repex_mpi.h>

#ifndef REPLICA_EXCHANGE_BASE_2D_S_EOFF_EDS_H
#define REPLICA_EXCHANGE_BASE_2D_S_EOFF_EDS_H

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
     class replica_exchange_base_2d_s_eoff_eds : public virtual replica_exchange_base_interface {
    public:
        replica_exchange_base_2d_s_eoff_eds(io::Argument _args,
                unsigned int cont,
                unsigned int globalThreadID,
                replica_graph_control &replicaGraphMPIControl,
                simulation::MpiControl &replica_mpi_control);
        /**
         * inits replica_reeds
         */
        void init() override;
        /**
         * Initialize file and data for eds_stat output.
         */
        void init_eds_stat();

    protected:
        virtual ~replica_exchange_base_2d_s_eoff_eds();

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
        * Just a wrapper for the actual swap functions swap_s and swap_eoff. Exchanges coordinates and averages.
        */
        //void swap() override;
        void determine_switch_probabilities() override;
        /**
        * Tries a swapping in s dimension of configuration if possible. Calculates energies, probabilities
        * and sends information via MPI communication if necessary.
        */
        void swap_s(const unsigned int partnerReplicaID);
        /**
        * Tries a swapping in eoff dimension of configuration if possible. Calculates energies, probabilities
        * and sends information via MPI communication if necessary.
        */
        void swap_eoff(const unsigned int partnerReplicaID);
        /**
        * Finds partner for current switch
        * @return ID of partner, own ID if no switching in current trial
        */
        int find_partner() const override;
        /**
        * Finds partner for current switch in eoff dimension for even number of eoffs and trial % 4 == 1
        * @return ID of partner, own ID if no switching in current trial
        */
        int partner_eoffDim_numEoffeven_firstCase() const;
        /**
        * Finds partner for current switch in eoff dimension for odd number of eoffs and trial % 4 == 1
        * @return ID of partner, own ID if no switching in current trial
        */
        int partner_eoffDim_numEoffodd_firstCase() const;
        /**
        * Finds partner for current switch in eoff dimension for even number of eoffs and trial % 4 == 3
        * @return ID of partner, own ID if no switching in current trial
        */
        int partner_eoffDim_numEoffeven_secondCase() const;
        /**
        * Finds partner for current switch in eoff dimension for odd number of eoffs and trial % 4 == 3
        * @return ID of partner, own ID if no switching in current trial
        */
        int partner_eoffDim_numEoffodd_secondCase() const;
        /**
        * Finds partner for current switch in eoff dimension for odd number of eoffs in cyclic manner
        * @return ID of partner, own ID if no switching in current trial
        */
        int partner_eoffDim_numEoffodd_cyclic() const;
        /**
        * calculates probability of switch with current partner in eoff dimension, may involve MPI communication
        */
        double calc_probability_for_eoff_exchange(const unsigned int partnerReplicaID);
        /**
        * Sets eds_struct() parameters to original value of replica
        */
        void reset_eds();
        /**
         * Sets  eds_struct() parameters to value of partner (for energy calculations)
         */
        void change_eds(const unsigned int partner);

        //TODO
        /*
        * energy calculation for statistical purposes of eds_stat() in replica_exchange_base.cc
        * for given configuration with new smoothing parameter s.
        */
        double calc_energy_eds_stat(double s);
        double calculate_energy_core();
        double calculate_energy(const unsigned int partner);
    };
}
#endif /* REPLICA_EXCHANGE_BASE_2D_S_EOFF_EDS_H */
