/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   replica_graph_control.h
 * Author: bschroed
 *
 * Created on May 4, 2020, 10:21 PM
 */

#ifndef REPLICA_GRAPH_CONTROL_H
#define REPLICA_GRAPH_CONTROL_H
#ifdef XXMPI
#include "mpi.h"
#endif

#include "stdheader.h"

namespace util{
    class replica_graph_control {
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

