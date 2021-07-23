/* 
 * File:   replica_graph_control.cpp
 * Author: bschroed
 * 
 * Created on May 4, 2020, 10:21 PM
 * Modified June 18, 2021 - bschroed, srieder
 */

#ifdef XXMPI
#include "mpi.h"
#endif
#include "stdheader.h"
#include "replica_graph_control.h"

namespace re{
    replica_graph_control::replica_graph_control() : graphID(0), masterID(0), threadID(-1), numberOfThreads(0) {
    }

    replica_graph_control::replica_graph_control(unsigned int graphID, unsigned int masterID,unsigned int threadID,
                                  unsigned int numberOfThreads, unsigned int numberOfReplicas,
                                  unsigned int mpiColor): graphID(graphID), masterID(masterID), threadID(threadID),
                                          numberOfThreads(numberOfThreads), numberOfReplicas(numberOfReplicas),
                                          mpiColor(mpiColor){}
#ifdef XXMPI
    replica_graph_control::replica_graph_control(unsigned int graphID, unsigned int masterID,unsigned int threadID,
                                  unsigned int numberOfThreads, unsigned int numberOfReplicas,
                                  unsigned int mpiColor, MPI_Comm comm): graphID(graphID), masterID(masterID), threadID(threadID),
                                          numberOfThreads(numberOfThreads), numberOfReplicas(numberOfReplicas),
                                          mpiColor(mpiColor), comm(comm){}
#endif
                                  
    replica_graph_control::replica_graph_control(const replica_graph_control& orig) {
    }

    replica_graph_control::~replica_graph_control() {
    }

    void replica_graph_control::print_struct(std::string stage){
        std::cerr << "\n\n REPLICA_GRAPH_CONSTRUCT\n";
        std::cerr <<stage<< "\tgraphID:\t"<<graphID <<"\n";
        std::cerr <<stage<< "\tmasterID:\t"<< masterID<<"\n";
        std::cerr <<stage<< "\tthreadID:\t"<<threadID <<"\n";
        std::cerr <<stage<< "\tnumberOfThreads:\t"<<numberOfThreads <<"\n";
        std::cerr <<stage<< "\tnumberOfReplicas:\t"<<numberOfReplicas <<"\n";
        std::cerr <<stage<< "\tmpiColor:\t"<<mpiColor <<"\n";
        
        std::cerr << "\tMasterThreads:\t";
        for (auto x: replicaMasterIDs){
            std::cerr <<stage<< " "<<x ;
        }
        std::cerr << "\n";
        
        std::cerr << "\tMasterThreads:\t\n";
        for (auto replicaThreadsVectors: replicaThreads){
            std::cerr <<stage<< "\t";
            for(auto x: replicaThreadsVectors){
                std::cerr << " "<<x;            
            }
            std::cerr << "\n";
        }
        std::cerr << "\n";         
    }
    
        void replica_graph_control::print_struct(){
        std::cerr << "\n\n REPLICA_GRAPH_CONSTRUCT\n";
        std::cerr << "\tgraphID:\t"<<graphID <<"\n";
        std::cerr << "\tmasterID:\t"<< masterID<<"\n";
        std::cerr << "\tthreadID:\t"<<threadID <<"\n";
        std::cerr << "\tnumberOfThreads:\t"<<numberOfThreads <<"\n";
        std::cerr << "\tnumberOfReplicas:\t"<<numberOfReplicas <<"\n";
        std::cerr << "\tmpiColor:\t"<<mpiColor <<"\n";
        
        std::cerr << "\tMasterThreads:\t";
        for (auto x: replicaMasterIDs){
            std::cerr << " "<<x ;
        }
        std::cerr << "\n";
        
        std::cerr << "\tMasterThreads:\t\n";
        for (auto replicaThreadsVectors: replicaThreads){
            std::cerr << "\t";
            for(auto x: replicaThreadsVectors){
                std::cerr << " "<<x;            
            }
            std::cerr << "\n";
        }
        std::cerr << "\n";         
    }
}