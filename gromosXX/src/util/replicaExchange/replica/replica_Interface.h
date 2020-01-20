/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   replica_Interface.h
 * Author: bschroed
 *
 * Created on November 8, 2019, 12:53 PM
 */

#include <stdheader.h>
#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <simulation/simulation.h>
#include <configuration/configuration.h>
#include <algorithm/algorithm/algorithm_sequence.h>
#include <interaction/interaction.h>
#include <interaction/forcefield/forcefield.h>
#include <io/configuration/out_configuration.h>
#include <io/argument.h>
#include <io/read_input.h>
#include <math/random.h>
#include <util/usage.h>
#include <util/replicaExchange/repex_mpi.h>

#ifndef REPLICA_INTERFACE_H
#define REPLICA_INTERFACE_H

namespace util{
    class replica_Interface {
    /**
       * replica Class Interface
       * All the data members are held public, so Master/Slave classes have full 
       * control over it. - Replica interace provied needed base functionality!
       */

    /**
        * ATTRIBUTES   
        */
      public:
        /**
         * ID, unique identifier
         */
        const unsigned int ID;
        /**
         * rank of node the replica is on
         */
        const unsigned int rank;

        /* NOTE:
         * These cannot be references, because every node has only one copy, so
         * the repicas would then reference the same conf/topo
         */
        /**
         * io::Argument object to read in all simulation parameters
         */
        io::Argument args;
        /**
         * trajectories
         */
        io::Out_Configuration * traj;
        /**
         * output file stream for out_configuration
         */
        std::ofstream * os;
        /**
         * configuration of replica
         */
        configuration::Configuration conf;
        /**
         * a temporary VArray for switching the lattice shifts
         */
        math::VArray * latticeTMP;
        /**
         * algorithm sequence
         */
        algorithm::Algorithm_Sequence md;
        /**
         * simulation
         */
        simulation::Simulation sim;
        /**
         * topology
         */
        topology::Topology topo;
        /**
         * random number generator
         */
        math::RandomGeneratorGSL * rng;
        /**
         * Temperature of replica
         */
        double T;
        /**
         * Lambda of replica
         */
        double l;
        /**
         * Timestep for current run
         */
        double dt;
        /**
         * potential energy using current Hamiltonian
         */
        double epot;
        /**
         * potential energy of the partner Hamiltonian (for bookkeeping)
         */
        double epot_partner;
        /**
         * ID of partner
         */
        unsigned int partner;

        /**
         * run number of replica
         */
        unsigned int run;

        /**
         * total number of runs per replica
         */
        unsigned int total_runs;

        /**
         * total number of steps simulated
         */
        unsigned int steps;
        /**
         * maximal number of steps to simulate per run
         */
        unsigned int maxSteps;

        /**
         * probability of last switch
         */
        double probability;
        /**
         * switched last time?
         */
        bool switched;
        /**
         * current simulation time
         */
        double time;

        /**
         * FUNCTIONS
         */
      public:
         
          replica_Interface(int ID, int rank, io::Argument args);
          virtual ~replica_Interface();
        /**
         * runs an MD simulation
         */
        virtual void run_MD() = 0;
        /**
         * init an MD simulation
         */
        virtual void init();
        /**
         * TODO: bring into Exchanger partially!!!
         * Finds partner for current switch
         * @return ID of partner, own ID if no switching in current trial
         */
        virtual int find_partner() const;
        /**
         * TODO: bring into Exchanger 
         * calculates probability of switch with current partner, may involve MPI communication
         * @param partnerID
         * @param partnerRank
         * @return probability
         */
        virtual double calc_probability(const int partnerID, const int partnerRank);
        /**TODO: bring into Exchanger 
         * calculates potential energy of current configuration with lambda(Hamiltonian) of partner
         * @param partner ID of partner
         * @return potential energy of configuration with lambda(Hamiltonian) of partner
         */
        virtual double calculate_energy(const int partner);
        /**
         * calculates potential energy for current configuration with own lambda
         */
        virtual double calculate_energy();
        /**TODO: bring into Exchanger 
         * finds out if configurations are to be switched; sets switched to true if so
         */
        virtual void swap(const unsigned int partner, const unsigned int partnerRank);
        /**TODO: bring into Exchanger 
         * switch back the averages
         */
         virtual void exchange_averages();

        /**
         * write final cnf
         */
        virtual void write_final_conf();

        /**TODO: bring into Exchanger 
         * Initiates MPI communication to receive new configuration information
         * @param senderID
         * @param senderRank
         */
        virtual void receive_new_coord(const int senderID, const int senderRank);
        
        /**TODO: bring into Exchanger 
         * Initiates MPI communication to send new configuration information
         * @param receiverID
         * @param senderRank
         */
        virtual void send_coord(const int receiverID, const int senderRank);
        /**TODO: bring into Exchanger 
         * Prints information of replica to std::cout for debugging purposes
         */
        virtual void print_info(std::string bla)const;
        /**
         * Scales the velocities after an exchange in temperature replica exchange
         */
         virtual void velscale(int i);

         
        protected:
        /**TODO: bring into Exchanger 
         * Sets lambda parameter to original value of replica
         */
        void set_lambda();
        /**TODO: bring into Exchanger 
         * Sets temperature to original value of replica
         */
        void set_temp();
        /**TODO: bring into Exchanger 
         * Sets lambda parameter to value of partner (for energy calculations)
         */
        void change_lambda(const unsigned int partner);
        /**TODO: bring into Exchanger 
         * Sets temperature to value of partner (for energy calculations)
         */
        void change_temp(const unsigned int partner);

      };
    }//namespace util
#endif /* REPLICA_INTERFACE_H */

