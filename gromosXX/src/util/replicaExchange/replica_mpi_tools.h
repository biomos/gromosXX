/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   mpi_tools.h
 * Author: bschroed
 *
 * Created on November 14, 2019, 2:12 PM
 */
#include<stdheader.h>

#ifndef MPI_TOOLS_H
#define MPI_TOOLS_H
#ifdef XXMPI
    #include <mpi.h>
    namespace util{
          struct replica_graph_mpi_control{

            replica_graph_mpi_control() : graphID(0), masterID(0), threadID(-1), numberOfThreads(0)
            { };

            unsigned int graphID;
            unsigned int masterID;
            unsigned int threadID;
            unsigned int numberOfThreads;
            unsigned int numberOfReplicas;

            MPI_Comm comm;
            int mpiColor;

            std::vector<unsigned int> replicaMasterIDs;
            std::vector<std::vector<unsigned int>> replicaThreads;
            std::map<unsigned int, unsigned int> threadReplicaMap;
        };

        //tools for thread tracking
       std::vector<std::vector<unsigned int > > calculate_Replica_Thread_Coordination(int rank, int totalNumberOfThreads, int numReplicas);

    }//namespace util
#else
    namespace util{
      struct replica_graph_mpi_control{

        replica_graph_mpi_control() : graphID(0), masterID(0), threadID(-1), numberOfThreads(0)
        { };

        unsigned int graphID;
        unsigned int masterID;
        unsigned int threadID;
        unsigned int numberOfThreads;
        unsigned int numberOfReplicas;

        int mpiColor;

        std::vector<unsigned int> replicaMasterIDs;
        std::vector<std::vector<unsigned int>> replicaThreads;
        std::map<unsigned int, unsigned int> threadReplicaMap;
    };

    //tools for thread tracking
   std::vector<std::vector<unsigned int > > calculate_Replica_Thread_Coordination(int rank, int totalNumberOfThreads, int numReplicas);

    }//namespace util
#endif
#endif /* MPI_TOOLS_H */

