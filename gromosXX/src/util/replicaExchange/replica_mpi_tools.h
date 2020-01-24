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

#ifndef MPI_TOOLS_H
#define MPI_TOOLS_H



namespace util{
    //tools for thred tracking
   std::vector<std::vector<unsigned int > > calculate_Replica_Thread_Coordination(int rank, int totalNumberOfThreads, int numReplicas);
    
}//namespace util

#endif /* MPI_TOOLS_H */

