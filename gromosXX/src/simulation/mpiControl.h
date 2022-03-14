/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   mpiControl.h
 * Author: bschroed
 *
 * Created on May 4, 2020, 11:04 AM
 */

#ifndef MPICONTROL_H
#define MPICONTROL_H
#ifdef XXMPI
    #include<mpi.h>
#endif

namespace simulation {
    
    class MpiControl {
    public:
        MpiControl();
        MpiControl(bool mpi);
        MpiControl(bool mpi, int simulationID, int numberOfThreads, int masterID, int threadID, int mpiColor);
        MpiControl(bool mpi, int simulationID, int numberOfThreads, int masterID, int threadID, int mpiColor, std::vector<unsigned int> simulationOwnedThreads);
        
        #ifdef XXMPI
            MpiControl(bool mpi, int simulationID, int numberOfThreads, int masterID, int threadID, int mpiColor, MPI_Comm comm);
            MpiControl(bool mpi, int simulationID, int numberOfThreads, int masterID, int threadID, int mpiColor, std::vector<unsigned int> simulationOwnedThreads, MPI_Comm comm);
        #endif
        
        MpiControl(const MpiControl& orig);
        virtual ~MpiControl();

        //Attributes: 
        bool mpi;
        int simulationID; //local replica id of simulation
        int numberOfThreads;    //total_number_of_threads      
        int masterID; //local master of this 
        int threadID;
        int mpiColor;
        std::vector<unsigned int> simulationOwnedThreads; 

        #ifdef XXMPI
            MPI_Comm comm; 
        #endif
        
        void print_struct();
        void print_struct(std::string stage);
    private:
        


    };
}

#endif /* MPICONTROL_H */

