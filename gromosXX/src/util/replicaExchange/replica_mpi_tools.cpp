/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "stdheader.h"
#include "replica_mpi_tools.h"


//std::vector<std::vector<int > >  util::calculate_Replica_Thread_Coordination(int rank, int totalNumberOfThreads, int numReplicas ){
        /**calculate_Replica_Thread_Coordination
         *  This function is calculating the Replica Simulation mapping to the threads.
         * @param rank  this is the rank of the current thread
         * @param totalNumberOfThreads  total number of threads 
         * @param numReplicas   number of Replica simulations.
         * @return 
         */
/*
       std::vector<std::vector<int > > repIDs;
       repIDs.resize(numReplicas);

       int threadsPerReplicaSimulation = totalNumberOfThreads / numReplicas;
       int leftOverThreads = totalNumberOfThreads % numReplicas;

       unsigned int threadID =0;
       int replica_offset = 0;
       std::cout << "left_over "<< leftOverThreads << "\n";
       for (int replicaSimulationID = 0; replicaSimulationID < numReplicas; replicaSimulationID++) {

           for (int replicaSimulationSubThread = 0; replicaSimulationSubThread < threadsPerReplicaSimulation; replicaSimulationSubThread++) {

               threadID = replicaSimulationSubThread + replicaSimulationID*threadsPerReplicaSimulation+replica_offset;
               repIDs[replicaSimulationID].push_back(threadID);

           }
           if(leftOverThreads>0 and replicaSimulationID < leftOverThreads){    //left over threads are evenly distirbuted.
               threadID = threadsPerReplicaSimulation + replicaSimulationID*totalNumberOfThreads+replica_offset;
               repIDs[replicaSimulationID].push_back(threadID);
               replica_offset++;

           }    
       }            
    return repIDs;
}*/