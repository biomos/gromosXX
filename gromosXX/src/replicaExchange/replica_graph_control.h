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

