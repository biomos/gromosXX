/*
 * File:   cuda_mpi.h
 * Author: Axel && Mel
 *
 * Created on October 20, 2011, 11:34 AM
 *
 * MPI specific CUDA include file
 */

#ifndef CUDA_MPI_H
#define	CUDA_MPI_H
#include <memory>
#include <new>
#include <iostream>

// define namespace
namespace gcuda {


void MPIstrider(unsigned int elements[], unsigned int size, int rank, int stride) {
    
    int stridecounter = rank;

    for (unsigned int i = 0; i < size; i++) {

        if ((stridecounter - i) != 0) {
            elements[i] = 0;
        }else  {
            stridecounter += stride;
        } 
    }

    return;
}

}//namespace


#endif	/* CUDA_MPI_H */