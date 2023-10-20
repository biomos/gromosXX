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
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   mpiControl.cpp
 * Author: bschroed
 * 
 * Created on May 4, 2020, 11:04 AM
 */
#ifdef XXMPI
    #include<mpi.h>
#endif
#include "stdheader.h"
#include "mpiControl.h"

namespace simulation {

    MpiControl::MpiControl(): mpi(false), simulationID(0), numberOfThreads(1), masterID(0), threadID(0), mpiColor(1) {

    }

    MpiControl::MpiControl(bool mpi): 
    mpi(mpi),  simulationID(0), numberOfThreads(1), masterID(0), threadID(0), mpiColor(0)
    {
        simulationOwnedThreads = std::vector<unsigned int>();
        simulationOwnedThreads.resize(numberOfThreads);
    }
    
    MpiControl::MpiControl(bool mpi, int simulationID, int numberOfThreads, int masterID, int threadID, int mpiColor): 
    mpi(mpi), simulationID(simulationID), numberOfThreads(numberOfThreads), masterID(masterID), threadID(threadID), mpiColor(mpiColor) 
    {
        simulationOwnedThreads = std::vector<unsigned int>();
        simulationOwnedThreads.resize(numberOfThreads);
        
    }
    
        MpiControl::MpiControl(bool mpi, int simulationID, int numberOfThreads, int masterID, int threadID, int mpiColor, std::vector<unsigned int> simulationOwnedThreads): 
    mpi(mpi), simulationID(simulationID), numberOfThreads(numberOfThreads), masterID(masterID), threadID(threadID), mpiColor(mpiColor), simulationOwnedThreads(simulationOwnedThreads)
    {
    }
        
#ifdef XXMPI
    MpiControl::MpiControl(bool mpi, int simulationID, int numberOfThreads, int masterID, int threadID, int mpiColor, MPI_Comm comm): 
    mpi(mpi), simulationID(simulationID), numberOfThreads(numberOfThreads), masterID(masterID), threadID(threadID), mpiColor(mpiColor), comm(comm)
    {
        simulationOwnedThreads = std::vector<unsigned int>();
        simulationOwnedThreads.resize(numberOfThreads);
        
    }
    MpiControl::MpiControl(bool mpi, int simulationID, int numberOfThreads, int masterID, int threadID, int mpiColor, std::vector<unsigned int> simulationOwnedThreads, MPI_Comm comm): 
    mpi(mpi), simulationID(simulationID), numberOfThreads(numberOfThreads), masterID(masterID), threadID(threadID), mpiColor(mpiColor), simulationOwnedThreads(simulationOwnedThreads), comm(comm)
    {   
    }
#endif
        


    MpiControl::MpiControl(const MpiControl& orig) {
    }

    MpiControl::~MpiControl() {
    }

    void MpiControl::print_struct(){
        std::cerr << "\n\nMD MPI Control\n";
        std::cerr << "\tmpi: "<<mpi<<"\n";
        std::cerr << "\tsimulationID: "<<simulationID<<"\n";
        std::cerr << "\tnumberOfThreads: "<<numberOfThreads<<"\n";
        std::cerr << "\tmasterID: "<<masterID<<"\n";
        std::cerr << "\tthreadID: "<<threadID<<"\n";
        std::cerr << "\tmpiColor: "<<mpiColor<<"\n";
        
        std::cerr<< "\tsimulationOwned Threads: ";
        for(auto x: simulationOwnedThreads){
            std::cerr<< " "<<x;
        }
        std::cerr << "\n";
    }


    void MpiControl::print_struct(std::string stage){
        std::cerr << "\n\nMD MPI Control\n";
        std::cerr <<stage<< "\tmpi: "<<mpi<<"\n";
        std::cerr <<stage<< "\tsimulationID: "<<simulationID<<"\n";
        std::cerr <<stage<< "\tnumberOfThreads: "<<numberOfThreads<<"\n";
        std::cerr <<stage<< "\tmasterID: "<<masterID<<"\n";
        std::cerr <<stage<< "\tthreadID: "<<threadID<<"\n";
        std::cerr <<stage<< "\tmpiColor: "<<mpiColor<<"\n";
        
        std::cerr <<stage<< "\tsimulationOwned Threads: ";
        for(auto x: simulationOwnedThreads){
            std::cerr<< " "<<x;
        }
        std::cerr << "\n";
   
    }
}