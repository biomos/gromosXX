/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   replica_Interface.cpp
 * Author: bschroed
 * 
 * Created on November 8, 2019, 12:53 PM
 */


#include <io/argument.h>
#include <util/error.h>
#include <math/volume.h>

#ifdef XXMPI
#include <mpi.h>
#endif
#include <replicaExchange/replica/_replica_Interface.h>

#undef MODULE
#undef SUBMODULE
#define MODULE re
#define SUBMODULE replica


re::replica_Interface::replica_Interface(int globalThreadID, simulation::MpiControl & _replica_mpi_control, io::Argument args):
        globalThreadID(globalThreadID),  
        simulationID(_replica_mpi_control.simulationID), simulationThreadID(_replica_mpi_control.threadID), replica_mpi_control(_replica_mpi_control), 
        args(args){
    
    
        //Attributes:  - strange Value copy
        replica_mpi_control.mpi = _replica_mpi_control.mpi;
        replica_mpi_control.simulationID = _replica_mpi_control.simulationID; //local replica id of simulation
        replica_mpi_control.numberOfThreads = _replica_mpi_control.numberOfThreads;    //total_number_of_threads      
        replica_mpi_control.masterID = _replica_mpi_control.masterID; //local master of this 
        replica_mpi_control.threadID = _replica_mpi_control.threadID;
        replica_mpi_control.mpiColor = _replica_mpi_control.mpiColor;
        replica_mpi_control.simulationOwnedThreads = _replica_mpi_control.simulationOwnedThreads; 
    #ifdef XXMPI
        replica_mpi_control.comm = _replica_mpi_control.comm; 
    #endif
}

re::replica_Interface::~replica_Interface(){
    delete rng;
    delete traj;
    delete os;
}
void re::replica_Interface::init() {
  // init MD simulation
  md.init(topo, conf, sim, *os, true);
}

void re::replica_Interface::write_final_conf() {
    traj->write(conf, topo, sim, io::final);
    delete  traj;   //bschroed needs to be done, so the flush functions are triggered
}

