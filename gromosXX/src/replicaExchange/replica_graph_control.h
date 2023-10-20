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
 * File:   replica_graph_control.h
 * Author: bschroed
 *
 * Created on May 4, 2020, 10:21 PM
 * Modified June 18, 2021 - bschroed, srieder
 */

#ifndef REPLICA_GRAPH_CONTROL_H
#define REPLICA_GRAPH_CONTROL_H
#ifdef XXMPI
#include "mpi.h"
#endif

#include "stdheader.h"

namespace re{
    class replica_graph_control {
        /*
         * This structure is used to facilitate the MPI-Kommunikations of the Replica Graph. 
         * it is initialized in repexMPI.a nd shall be used in all RE-Graph MPI calls.
         */
        public:

            replica_graph_control();
            
            replica_graph_control(unsigned int graphID, unsigned int masterID,unsigned int threadID,
                                  unsigned int numberOfThreads, unsigned int numberOfReplicas,
                                  unsigned int mpiColor);

        #ifdef XXMPI
            replica_graph_control(unsigned int graphID, unsigned int masterID,unsigned int threadID,
                                  unsigned int numberOfThreads, unsigned int numberOfReplicas,
                                  unsigned int mpiColor, MPI_Comm comm);
        #endif

            
            replica_graph_control(const replica_graph_control& orig);
            virtual ~replica_graph_control();


            unsigned int graphID;
            unsigned int masterID;
            unsigned int threadID;
            unsigned int numberOfThreads;
            unsigned int numberOfReplicas;
            unsigned int mpiColor;

            std::vector<unsigned int> replicaMasterIDs = std::vector<unsigned int>();
            std::vector<std::vector<unsigned int>> replicaThreads = std::vector<std::vector<unsigned int>>();
            std::map<unsigned int, unsigned int> threadReplicaMap = std::map<unsigned int, unsigned int>();

            #ifdef XXMPI
                MPI_Comm comm;
            #endif

            void print_struct();
            void print_struct(std::string stage);
        private:

    };
}
#endif /* REPLICA_GRAPH_CONTROL_H */

