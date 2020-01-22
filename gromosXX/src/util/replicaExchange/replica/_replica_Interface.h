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
        // MPI
        /**
         * ID, unique identifier
         */
        const unsigned int ID;
        /**
         * rank of node the replica is on
         */
        const unsigned int rank;
        
        unsigned int curentStepNumber=0;
        unsigned int stepsPerRun=0;
        unsigned int totalStepNumber=0;
        
        /* NOTE:
         * These cannot be references, because every node has only one copy, so
         * the replicas would then reference the same conf/topo
         * TODO: NOT TRUE ANYMORE bschroed
         */
         //SIMULATION PARAMS
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
        
        // IO
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
        
        //UTILS
        /**
         * random number generator
         */
        math::RandomGeneratorGSL * rng;
        
    /**
     * FUNCTIONS
     */
    public:
          
        //Constructor
        replica_Interface(int ID, int rank, io::Argument args);
        
        //Destructor
        virtual ~replica_Interface();
        /**
         * init an MD simulation
         */
        virtual void init();
        /**
         * runs an MD simulation
         */
        virtual void run_MD() = 0;
        /**
         * write final cnf
         */
        virtual void write_final_conf();
      };
    }//namespace util
#endif /* REPLICA_INTERFACE_H */

