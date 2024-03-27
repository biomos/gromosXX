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

#undef MODULE
#undef SUBMODULE
#define MODULE re
#define SUBMODULE replica_exchanger

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