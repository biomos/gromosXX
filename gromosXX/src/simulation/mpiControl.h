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

