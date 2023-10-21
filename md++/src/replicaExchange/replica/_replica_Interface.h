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
 * File:   replica_Interface.h
 * Author: bschroed
 *
 * Created on November 8, 2019, 12:53 PM
 * Modified June 18, 2021 - bschroed, srieder
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
#include <replicaExchange/repex_mpi.h>

#ifndef REPLICA_INTERFACE_H
#define REPLICA_INTERFACE_H

namespace re{
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
        const unsigned int globalThreadID;
        const unsigned int simulationID;
        const unsigned int simulationThreadID;
        simulation::MpiControl replica_mpi_control;
        
        /**
         * rank of node the replica is on
         */
        
        unsigned int curentStepNumber=0;
        unsigned int stepsPerRun=0;
        unsigned int totalStepNumber=0;

        double energy; // not activly used currently
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
        replica_Interface(int globalThreadID, simulation::MpiControl & _replica_mpi_control, io::Argument args);
        
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
         * update replica energies
         */
        virtual double calculateEnergies();
        /**
         * write final cnf
         */
        virtual void write_final_conf();


        /**
        * contains original forces which have to be reset after RE-EDS exchange energy calculation
        */
        math::VArray force_orig;
        /**
         * contains original virial which have to be reset after RE-EDS exchange energy calculation
         */
        math::Matrix virial_tensor_orig;

      };
    }//namespace re
#endif /* REPLICA_INTERFACE_H */

