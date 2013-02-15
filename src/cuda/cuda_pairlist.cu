/*
 * CUDA PAIRLIST KERNELS
 *
 */

#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <map>
#include <set>
#include <vector>
#include <stdexcept>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <typeinfo>
#include <list>
#include <stdio.h>

#include "../../config.h"

#include "../math/gmath.h"
#include "../util/debug.h"

#include "../io/message.h"
#include "../io/gzstream.h"

#include "../util/timing.h"


#include "c_pairlist.h"
#include "gromos_cuda.h"
#include </usr/local/cuda/include/cuda.h>
#include </usr/local/cuda/include/driver_types.h>
#include </usr/local/cuda/include/cuda_runtime_api.h>
#include </usr/local/cuda/include/driver_functions.h>
#include </usr/local/cuda/include/driver_types.h>


extern "C" int gcuda::has_gpu(){

    int number_gpus;
    if ( cudaSuccess != cudaGetDeviceCount(&number_gpus)) {
        return -1;
    } else {
        return number_gpus;
    }
}

extern "C" gcuda::devices gcuda::get_dev_properties(int number_gpus) {

    // initialise devices structure
    gcuda::devices gpus;
    gpus.number_devices = number_gpus;
    cudaDeviceProp device;

    // loop over devices
    gpus.error = false; // starting without an error
    for (unsigned int i = 0; i < number_gpus; i++) {
    if ( cudaSuccess != cudaGetDeviceProperties(&device, i)) {
            gpus.error = true; // error :(
        } else {
        gpus.property.push_back(device);
        }
    }
    return gpus;
}

// declaration of constant values
    /* box -- an array would be a performance problem as multiple threads would
    * access different elements and cause additional read cycles
    */
__constant__ double box_00, box_01, box_02, box_10, box_11, box_12, box_20, box_21, box_22, cs, cs2, cl, cl2;
__constant__ unsigned int na, nc, nsc, nsol;


extern "C" bool copybox(math::Box box)
{
  // copy values
  cudaMemcpyToSymbol("box_00", &box(0)(0), 8, 0, cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol("box_01", &box(0)(1), 8, 0, cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol("box_02", &box(0)(2), 8, 0, cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol("box_10", &box(1)(0), 8, 0, cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol("box_11", &box(1)(1), 8, 0, cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol("box_12", &box(1)(2), 8, 0, cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol("box_20", &box(2)(0), 8, 0, cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol("box_21", &box(2)(1), 8, 0, cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol("box_22", &box(2)(2), 8, 0, cudaMemcpyHostToDevice);
  return true;
}

extern "C" bool copycutoff(double *cutoff)
{
  // copy values
  cudaMemcpyToSymbol("cs", &cutoff[0], 8, 0, cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol("cs2", &cutoff[1], 8, 0, cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol("cl", &cutoff[2], 8, 0, cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol("cl2", &cutoff[3], 8, 0, cudaMemcpyHostToDevice);
  return true;
}

extern "C" bool copynratoms(unsigned int nratoms )
{
  // copy values
  cudaMemcpyToSymbol("na", &nratoms, 4, 0, cudaMemcpyHostToDevice);
  return true;
}

extern "C" bool copynrcg(unsigned int nrcg )
{
  // copy values
  cudaMemcpyToSymbol("nc", &nrcg, 4, 0, cudaMemcpyHostToDevice);
  return true;
}

extern "C" bool copynrsolutecg(unsigned int nrsolutecg )
{
  // copy values
  cudaMemcpyToSymbol("nsc", &nrsolutecg, 4, 0, cudaMemcpyHostToDevice);
  return true;
}

extern "C" bool copynrsolute(unsigned int nrsolute )
{
  // copy values
  cudaMemcpyToSymbol("nsol", &nrsolute, 4, 0, cudaMemcpyHostToDevice);
  return true;
}
/*
 * Pointers to the allocated memory in constant memory
 *
 */
__constant__ size_t positions;
__constant__ size_t cogs;
__constant__ size_t chargegroups, exclusions, exclusions_elements, solute_short, solute_short_elements, solute_long, solute_long_elements, solvent_short, solvent_short_elements, solvent_long, solvent_long_elements;
__constant__ size_t positions_pitch, chargegroups_pitch, exclusions_pitch, exclusions_elements_pitch, cogs_pitch, solute_short_pitch, solute_short_elements_pitch, solute_long_pitch, solute_long_elements_pitch, solvent_short_pitch, solvent_short_elements_pitch, solvent_long_pitch, solvent_long_elements_pitch;

extern "C" bool copymemorymap(gcuda::memmap &gpumemory, bool atomic) {

  // doubles
  cudaMemcpyToSymbol("positions", &gpumemory.memory[gpumemory.map["positions"]] , 8, 0, cudaMemcpyHostToDevice);

  if (!atomic) {
  cudaMemcpyToSymbol("cogs", &gpumemory.memory[gpumemory.map["cogs"]] , 8, 0, cudaMemcpyHostToDevice);
  }

  // ints
  if (!atomic) {
  cudaMemcpyToSymbol("chargegroups", &gpumemory.memory[gpumemory.map["chargegroups"]] , 8, 0, cudaMemcpyHostToDevice);
  }


  cudaMemcpyToSymbol("exclusions", &gpumemory.memory[gpumemory.map["exclusions"]] , 8, 0, cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol("exclusions_elements", &gpumemory.memory[gpumemory.map["exclusions_elements"]] , 8, 0, cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol("solute_short", &gpumemory.memory[gpumemory.map["solute_short"]] , 8, 0, cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol("solute_short_elements", &gpumemory.memory[gpumemory.map["solute_short_elements"]] , 8, 0, cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol("solute_long", &gpumemory.memory[gpumemory.map["solute_long"]] , 8, 0, cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol("solute_long_elements", &gpumemory.memory[gpumemory.map["solute_long_elements"]] , 8, 0, cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol("solvent_short", &gpumemory.memory[gpumemory.map["solvent_short"]] , 8, 0, cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol("solvent_short_elements", &gpumemory.memory[gpumemory.map["solvent_short_elements"]] , 8, 0, cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol("solvent_long", &gpumemory.memory[gpumemory.map["solvent_long"]] , 8, 0, cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol("solvent_long_elements", &gpumemory.memory[gpumemory.map["solvent_long_elements"]] , 8, 0, cudaMemcpyHostToDevice);

 
  // pitches
  cudaMemcpyToSymbol("positions_pitch", &gpumemory.pitches[gpumemory.map["positions"]] , sizeof(size_t), 0, cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol("cogs_pitch", &gpumemory.pitches[gpumemory.map["cogs"]] , sizeof(size_t), 0, cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol("chargegroups_pitch", &gpumemory.pitches[gpumemory.map["chargegroups"]] , sizeof(size_t), 0, cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol("exclusions_pitch", &gpumemory.pitches[gpumemory.map["exclusions"]] , sizeof(size_t), 0, cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol("exclusions_elements_pitch", &gpumemory.pitches[gpumemory.map["exclusions_elements"]] , sizeof(size_t), 0, cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol("solute_short_pitch", &gpumemory.pitches[gpumemory.map["solute_short"]] , sizeof(size_t), 0, cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol("solute_short_elements_pitch", &gpumemory.pitches[gpumemory.map["solute_short_elements"]] , sizeof(size_t), 0, cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol("solute_long_pitch", &gpumemory.pitches[gpumemory.map["solute_long"]] , sizeof(size_t), 0, cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol("solute_long_elements_pitch", &gpumemory.pitches[gpumemory.map["solute_long_elements"]] , sizeof(size_t), 0, cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol("solvent_short_pitch", &gpumemory.pitches[gpumemory.map["solvent_short"]] , sizeof(size_t), 0, cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol("solvent_short_elements_pitch", &gpumemory.pitches[gpumemory.map["solvent_short_elements"]] , sizeof(size_t), 0, cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol("solvent_long_pitch", &gpumemory.pitches[gpumemory.map["solvent_long"]] , sizeof(size_t), 0, cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol("solvent_long_elements_pitch", &gpumemory.pitches[gpumemory.map["solvent_long_elements"]] , sizeof(size_t), 0, cudaMemcpyHostToDevice);


  return true;

}

// device functions-------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

__device__ bool is_excluded_solute (unsigned int atom_a, unsigned int atom_b) {

    
    unsigned int* size = (unsigned int*)(exclusions_elements + atom_a * exclusions_elements_pitch);
    for (unsigned int i = 0; i < *size; i++) {
        unsigned int* ex = (unsigned int*)(exclusions + atom_a * exclusions_pitch) + i;
            if (*ex == atom_b) {
              return true;
            } 
    }
    return false;
}

__device__ unsigned int is_excluded_solvent (unsigned int me) {

    if ((na-(me-1))%3 == 0) {
            // first H
        return 1;
    }

    else if ((na-(me-1))%3 == 2) {
            // second H
        return 2;
    }

    else if ((na-(me-1))%3 == 1) {
            // oxygen H
        return 3;
    }

    else {
        // muh
        return 0;
    }

}

// kernels ---------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

__global__ void calc_cogs () {

   //calculate cogs
    //get my number
    unsigned int me = threadIdx.x+(blockDim.x*blockIdx.x);

    if (me < nc) {
        // cog vector
        double cog_x = 0.0;
        double cog_y = 0.0;
        double cog_z = 0.0;
        double* cog_x_ptr;
        double* cog_y_ptr;
        double* cog_z_ptr;
        double* rcx;
        double* rcy;
        double* rcz;

        unsigned int* mystart = (unsigned int*)((char*)chargegroups + me * chargegroups_pitch) + 0;
        unsigned int* myend = (unsigned int*)((char*)chargegroups + me * chargegroups_pitch) + 1;

        if (me < nsc) {
        //add
        for (unsigned int i = *mystart; i < *myend; i++) {
            cog_x_ptr = (double*)((char*)positions + i * positions_pitch) + 0;
            cog_x += *cog_x_ptr;
            cog_y_ptr = (double*)((char*)positions + i * positions_pitch) + 1;
            cog_y += *cog_y_ptr;
            cog_z_ptr = (double*)((char*)positions + i * positions_pitch) + 2;
            cog_z += *cog_z_ptr;
        }

        //devide
        cog_x = (cog_x / (*myend - *mystart));
        cog_y = (cog_y / (*myend - *mystart));
        cog_z = (cog_z / (*myend - *mystart));

        } else {
            cog_x_ptr = (double*)((char*)positions + *mystart * positions_pitch) + 0;
            cog_x += *cog_x_ptr;
            cog_y_ptr = (double*)((char*)positions + *mystart * positions_pitch) + 1;
            cog_y += *cog_y_ptr;
            cog_z_ptr = (double*)((char*)positions + *mystart * positions_pitch) + 2;
            cog_z += *cog_z_ptr;
        }

        //save
         rcx = (double*)((char*)cogs + me * cogs_pitch) + 0;
         *rcx = cog_x;
         rcy = (double*)((char*)cogs + me * cogs_pitch) + 1;
         *rcy = cog_y;
         rcz = (double*)((char*)cogs + me * cogs_pitch) + 2;
         *rcz = cog_z;
    }
}



__global__ void calc_pairlist_cg_self () {

   //calculate pairlist
    //get my number
    unsigned int me = threadIdx.x+(blockDim.x*blockIdx.x);
     if (me < nsc) {
    // get my cg

    unsigned int* mystart = (unsigned int*)((char*)chargegroups + me * chargegroups_pitch) + 0;
    unsigned int* myend = (unsigned int*)((char*)chargegroups + me * chargegroups_pitch) + 1;


    // add self
       
           
            for (unsigned int i = *mystart; i < *myend; i++) {
                for (unsigned int j = i+1; j < *myend; j++) {
                    if(!is_excluded_solute(i, j)) {
                       unsigned int* re = (unsigned int*)((char*)solute_short_elements + i * solute_short_elements_pitch);
                       unsigned int* rp = (unsigned int*)((char*)solute_short + i * solute_short_pitch) + *re;
                       *rp = j;
                       *re = *re + 1;
                    }
                }
            }
        }
}

__global__ void calc_cg_pairlist_solute() {

    unsigned long me = (blockIdx.x + (blockIdx.y*gridDim.x));
    if (me < nsc) {

       __shared__ unsigned int mystart;
       __shared__ unsigned int myend;
       __shared__ double mypos[3];
       __shared__ double pos[128][3];
       __shared__ unsigned int res_short[128][2];
       __shared__ unsigned int res_long[128][2];
       __shared__ int res_counter_short;
       __shared__ int res_counter_long;
       __shared__ unsigned int pos_counter;
       __shared__ unsigned int counters[40][2];

       //offset
       unsigned int offset = threadIdx.x+1;
       // distances
        double dis_x;
        double dis_y;
        double dis_z;
        // local result counter
        unsigned int myres;


       if (threadIdx.x == 0) {
        // save mycg
        mystart = *((unsigned int*)((char*)chargegroups + me * chargegroups_pitch) + 0);
        myend = *((unsigned int*)((char*)chargegroups + me * chargegroups_pitch) + 1);
        //save my position
        mypos[0] = *((double*)((char*)cogs + me * cogs_pitch) + 0);
        mypos[1] = *((double*)((char*)cogs + me * cogs_pitch) + 1);
        mypos[2] = *((double*)((char*)cogs + me * cogs_pitch) + 2);
        // initialize
       res_counter_short = 0;
       res_counter_long = 0;
       pos_counter = 0;
       }
       // cache
         while ((me+offset) < nc) {
       //find my neighbours
       pos[threadIdx.x][0] = *((double*)((char*)cogs + (me+offset) * cogs_pitch) + 0);
       pos[threadIdx.x][1] = *((double*)((char*)cogs + (me+offset) * cogs_pitch) + 1);
       pos[threadIdx.x][2] = *((double*)((char*)cogs + (me+offset) * cogs_pitch) + 2);
       atomicAdd(&pos_counter, 1);
       // sync us
       __syncthreads();

            if (threadIdx.x < pos_counter) {
                dis_x = mypos[0] - pos[threadIdx.x][0];
                dis_x -= box_00 * rint(dis_x * (1.0 / box_00));
                dis_y = mypos[1] - pos[threadIdx.x][1];
                dis_y -= box_11 * rint(dis_y * (1.0 / box_11));
                dis_z = mypos[2] - pos[threadIdx.x][2];
                dis_z -= box_22 * rint(dis_z * (1.0 / box_22));

            double dis_2 = (dis_x * dis_x) + (dis_y * dis_y) + (dis_z * dis_z);

            //compare dis_2 with cutoff
            if (dis_2 < cs2) {
                //short
                myres = atomicAdd(&res_counter_short, 1);
                res_short[myres][0] = *((unsigned int*)((char*)chargegroups + (me+offset) * chargegroups_pitch) + 0);
                res_short[myres][1] = *((unsigned int*)((char*)chargegroups + (me+offset) * chargegroups_pitch) + 1);
             } else if (dis_2 > cs2 && dis_2 < cl2) {
             //long
                myres = atomicAdd(&res_counter_long, 1);
                res_long[myres][0] = *((unsigned int*)((char*)chargegroups + (me+offset) * chargegroups_pitch) + 0);
                res_long[myres][1] = *((unsigned int*)((char*)chargegroups + (me+offset) * chargegroups_pitch) + 1);
             }
            //done comparing
            // sync us
            __syncthreads();
            } // if thread < pos
           // get couters
           if (threadIdx.x < ((myend-mystart)+1)) {
              counters[threadIdx.x][0] = *((unsigned int*)((char*)solute_short_elements + (mystart+threadIdx.x) * solute_short_elements_pitch));
              counters[threadIdx.x][1] = *((unsigned int*)((char*)solute_long_elements + (mystart+threadIdx.x) * solute_long_elements_pitch));
           }

             if (threadIdx.x == 0) {
            //////////////////
                 int mycounter_short = 0;
                 int mycounter_long = 0;
                 int ccount;
                    for (unsigned int j = mystart; j < myend; j++) {
                        ccount = j-mystart;
                         for (unsigned int k = res_short[ccount][0]; k < res_short[ccount][1]; k++) {
                              if(!is_excluded_solute(j, k)) {
                                    *((unsigned int*)((char*)solute_short + j * solute_short_pitch) + (counters[ccount][0]+mycounter_short)) = k;
                                    mycounter_short++;
                              }
                          }
                         for (unsigned int k = res_long[ccount][0]; k < res_long[ccount][1]; k++) {
                              *((unsigned int*)((char*)solute_long + j * solute_long_pitch) + (counters[ccount][1]+mycounter_long)) = k;
                              mycounter_long++;
                          }
                         *((unsigned int*)((char*)solute_short_elements + j * solute_short_elements_pitch)) = (counters[ccount][0]+mycounter_short);
                         *((unsigned int*)((char*)solute_long_elements + j * solute_long_elements_pitch)) = (counters[ccount][1]+mycounter_long);
                         mycounter_short = 0;
                         mycounter_long = 0;
                   }
                ////
            // reset
            res_counter_short = 0;
            res_counter_long = 0;
            pos_counter = 0;
            }
            // done
            // sync us
            __syncthreads();

    offset = offset+128;
    }// while have atoms

  } // if end
}

__global__ void calc_cg_pairlist_solvent() {

    unsigned long me = (blockIdx.x + (blockIdx.y*gridDim.x)) + nsc;
    if (me < nc) {

       __shared__ unsigned int mystart;
       __shared__ unsigned int myend;
       __shared__ double mypos[3];
       __shared__ double pos[128][3];
       __shared__ unsigned int res_short[128][2];
       __shared__ unsigned int res_long[128][2];
       __shared__ int res_counter_short;
       __shared__ int res_counter_long;
       __shared__ unsigned int pos_counter;
       __shared__ unsigned int counters[40][2];

       //offset
       unsigned int offset = threadIdx.x+1;
       // distances
        double dis_x;
        double dis_y;
        double dis_z;
        // local result counter
        unsigned int myres;


       if (threadIdx.x == 0) {
        // save mycg
        mystart = *((unsigned int*)((char*)chargegroups + me * chargegroups_pitch) + 0);
        myend = *((unsigned int*)((char*)chargegroups + me * chargegroups_pitch) + 1);
        //save my position
        mypos[0] = *((double*)((char*)cogs + me * cogs_pitch) + 0);
        mypos[1] = *((double*)((char*)cogs + me * cogs_pitch) + 1);
        mypos[2] = *((double*)((char*)cogs + me * cogs_pitch) + 2);
        // initialize
       res_counter_short = 0;
       res_counter_long = 0;
       pos_counter = 0;
       }
       // cache
         while ((me+offset) < nc) {
       //find my neighbours
       pos[threadIdx.x][0] = *((double*)((char*)cogs + (me+offset) * cogs_pitch) + 0);
       pos[threadIdx.x][1] = *((double*)((char*)cogs + (me+offset) * cogs_pitch) + 1);
       pos[threadIdx.x][2] = *((double*)((char*)cogs + (me+offset) * cogs_pitch) + 2);
       atomicAdd(&pos_counter, 1);
       // sync us
       __syncthreads();

            if (threadIdx.x < pos_counter) {
                dis_x = mypos[0] - pos[threadIdx.x][0];
                dis_x -= box_00 * rint(dis_x * (1.0 / box_00));
                dis_y = mypos[1] - pos[threadIdx.x][1];
                dis_y -= box_11 * rint(dis_y * (1.0 / box_11));
                dis_z = mypos[2] - pos[threadIdx.x][2];
                dis_z -= box_22 * rint(dis_z * (1.0 / box_22));

            double dis_2 = (dis_x * dis_x) + (dis_y * dis_y) + (dis_z * dis_z);

            //compare dis_2 with cutoff
            if (dis_2 < cs2) {
                //short
                myres = atomicAdd(&res_counter_short, 1);
                res_short[myres][0] = *((unsigned int*)((char*)chargegroups + (me+offset) * chargegroups_pitch) + 0);
                res_short[myres][1] = *((unsigned int*)((char*)chargegroups + (me+offset) * chargegroups_pitch) + 1);
             } else if (dis_2 > cs2 && dis_2 < cl2) {
             //long
                myres = atomicAdd(&res_counter_long, 1);
                res_long[myres][0] = *((unsigned int*)((char*)chargegroups + (me+offset) * chargegroups_pitch) + 0);
                res_long[myres][1] = *((unsigned int*)((char*)chargegroups + (me+offset) * chargegroups_pitch) + 1);
             } 
            //done comparing
            // sync us
            __syncthreads();
            } // if thread < pos
           // get couters
           if (threadIdx.x < ((myend-mystart)+1)) {
              counters[threadIdx.x][0] = *((unsigned int*)((char*)solvent_short_elements + (mystart+threadIdx.x) * solvent_short_elements_pitch));  
              counters[threadIdx.x][1] = *((unsigned int*)((char*)solvent_long_elements + (mystart+threadIdx.x) * solvent_long_elements_pitch)); 
           }
       
             if (threadIdx.x == 0) {
            //////////////////
                 int mycounter_short = 0;
                 int mycounter_long = 0;
                 int ccount;
                 int rcounts = 0;
                 int rcountl = 0;

                    for (unsigned int j = mystart; j < myend; j++) {
                        ccount = j - mystart;
                        if (rcounts < res_counter_short) {
                         for (unsigned int k = res_short[0][0]; k < res_short[0][1]; k++) {
                              *((unsigned int*)((char*)solvent_short + j * solvent_short_pitch) + (counters[ccount][0]+mycounter_short)) = k;
                              mycounter_short++;
                          }
                         rcounts++;
                        }
                        if (rcountl < res_counter_long) {
                         for (unsigned int k = res_long[rcountl][0]; k < res_long[rcountl][1]; k++) {
                              *((unsigned int*)((char*)solvent_long + j * solvent_long_pitch) + (counters[ccount][0]+mycounter_long)) = k;
                              mycounter_long++;
                          }
                         rcountl++;
                        }
                         *((unsigned int*)((char*)solvent_short_elements + j * solvent_short_elements_pitch)) = (counters[ccount][0]+mycounter_short);
                         *((unsigned int*)((char*)solvent_long_elements + j * solvent_long_elements_pitch)) = (counters[ccount][1]+mycounter_long);
                         mycounter_short = 0;
                         mycounter_long = 0;
                   }
                ////
            // reset
            res_counter_short = 0;
            res_counter_long = 0;
            pos_counter = 0;
            }
            // done
            // sync us
            __syncthreads();

    offset = offset+128;
    }// while have atoms

  } // if end
}

__global__ void calc_cg_pairlist() {

    unsigned int me = threadIdx.x+(blockDim.x*blockIdx.x);
    if (me < nc) {
        //save my position
        double* my_x = (double*)((char*)cogs + me * cogs_pitch) + 0;
        double* my_y = (double*)((char*)cogs + me * cogs_pitch) + 1;
        double* my_z = (double*)((char*)cogs + me * cogs_pitch) + 2;
        //find my neighbours
        double dis_x;
        double dis_y;
        double dis_z;

        unsigned int* mystart = (unsigned int*)((char*)chargegroups + me * chargegroups_pitch) + 0;
        unsigned int* myend = (unsigned int*)((char*)chargegroups + me * chargegroups_pitch) + 1;

        // add self
         if (me < nsc) {
            //im a solute atom
            for (unsigned int i = *mystart; i < *myend; i++) {
                for (unsigned int j = i+1; j < *myend; j++) {
                    if(!is_excluded_solute(i, j)) {
                       unsigned int* re = (unsigned int*)((char*)solute_short_elements + i * solute_short_elements_pitch);
                       unsigned int* rp = (unsigned int*)((char*)solute_short + i * solute_short_pitch) + *re;
                       *rp = j;
                       *re = *re + 1;
                    }
                }
            }
        } 

        for (unsigned int i = me + 1; i < nc; i++) {
            double* ox = (double*)((char*)cogs + i * cogs_pitch) + 0;
            dis_x = *my_x - *ox;
            dis_x -= box_00 * rint(dis_x * (1.0 / box_00));

            double* oy = (double*)((char*)cogs + i * cogs_pitch) + 1;
            dis_y = *my_y - *oy;
            dis_y -= box_11 * rint(dis_y * (1.0 / box_11));

            double* oz = (double*)((char*)cogs + i * cogs_pitch) + 2;
            dis_z = *my_z - *oz;
            dis_z -= box_22 * rint(dis_z * (1.0 / box_22));

            double dis_2 = (dis_x * dis_x) + (dis_y * dis_y) + (dis_z * dis_z);

            //compare dis_2 with cutoff
            if (dis_2 < cs2) {
                //short
                if (me < nsc) {
                    //im a solute atom
                    unsigned int* start = (unsigned int*)((char*)chargegroups + i * chargegroups_pitch) + 0;
                    unsigned int* end = (unsigned int*)((char*)chargegroups + i * chargegroups_pitch) + 1;

                    for (unsigned int j = *mystart; j < *myend; j++) {
                              for (unsigned int k = *start; k < *end; k++) {
                                     if(!is_excluded_solute(j, k)) {  // short range solvent?
                                        unsigned int* re = (unsigned int*)((char*)solute_short_elements + j * solute_short_elements_pitch);
                                         unsigned int* rp = (unsigned int*)((char*)solute_short + j * solute_short_pitch) + *re;
                                         *rp = k;
                                         *re = *re + 1;
                                    }
                             }
                    }
                } else {
                    // im a solvent atom
                    unsigned int* start = (unsigned int*)((char*)chargegroups + i * chargegroups_pitch) + 0;
                    unsigned int* end = (unsigned int*)((char*)chargegroups + i * chargegroups_pitch) + 1;
                    
                    for (unsigned int j = *mystart; j < *myend; j++) {
                         for (unsigned int k = *start; k < *end; k++) {
                              unsigned int* re = (unsigned int*)((char*)solvent_short_elements + j * solvent_short_elements_pitch);
                              unsigned int* rp = (unsigned int*)((char*)solvent_short + j * solvent_short_pitch) + *re;
                              *rp = k;
                              *re = *re + 1;
                          }
                   }
            }

         } else if (dis_2 > cs2 && dis_2 < cl2) {
                //long
                if (me < nsc) {
                    //im a solute atom
                    unsigned int* start = (unsigned int*)((char*)chargegroups + i * chargegroups_pitch) + 0;
                    unsigned int* end = (unsigned int*)((char*)chargegroups + i * chargegroups_pitch) + 1;

                    for (unsigned int j = *mystart; j < *myend; j++) {
                         for (unsigned int k = *start; k < *end; k++) {
                              unsigned int* re = (unsigned int*)((char*)solute_long_elements + j * solute_long_elements_pitch);
                              unsigned int* rp = (unsigned int*)((char*)solute_long + j * solute_long_pitch) + *re;
                              *rp = k;
                              *re = *re + 1;
                         }
                    }
                } else {
                    // im a solvent atom
                    unsigned int* start = (unsigned int*)((char*)chargegroups + i * chargegroups_pitch) + 0;
                    unsigned int* end = (unsigned int*)((char*)chargegroups + i * chargegroups_pitch) + 1;

                    for (unsigned int j = *mystart; j < *myend; j++) {
                         for (unsigned int k = *start; k < *end; k++) {
                              unsigned int* re = (unsigned int*)((char*)solvent_long_elements + j * solvent_long_elements_pitch);
                              unsigned int* rp = (unsigned int*)((char*)solvent_long + j * solvent_long_pitch) + *re;
                              *rp = k;
                              *re = *re + 1;
                         }
                    }
                }
            }
    } // for loop

    } // if end
}


__global__ void calc_pairlist_solute () {

   //calculate pairlist
    //get my number
    unsigned long me = blockIdx.x + (blockIdx.y*gridDim.x);
    if (me < nsol) {
    // cache
    __shared__ double pos[128][3];
    __shared__ unsigned int res_short[128];
    __shared__ unsigned int res_long[128];
    __shared__ int res_counter_short;
    __shared__ int res_counter_long;
    __shared__ double mypos[3];
    __shared__ unsigned int pos_counter;
    __shared__ unsigned int elements_short ;
    __shared__ unsigned int elements_long ;
    __shared__ unsigned int myexcl[40];
    __shared__ unsigned int num_excl;

    // distances
    double dis_x;
    double dis_y;
    double dis_z;

    // exculsion
    bool isex;

    // local result counter
    unsigned int myres;
    // loop konstant
    unsigned int offset = threadIdx.x + 1;
    
    //save my position
    if (threadIdx.x == 0) {
    // get positions once
    mypos[0] = *((double*)((char*)positions + me * positions_pitch) + 0);
    mypos[1] = *((double*)((char*)positions + me * positions_pitch) + 1);
    mypos[2] = *((double*)((char*)positions + me * positions_pitch) + 2);

    // initialize
    res_counter_short = 0;
    res_counter_long = 0;
    pos_counter = 0;
    elements_short = 0;
    elements_long = 0;
    num_excl = *((unsigned int*)(exclusions_elements + me * exclusions_elements_pitch));
    } else if (threadIdx.x > 0 && threadIdx.x < 20) { //FALSCH!!!!!! -Mel-
    // chache exclusions
     myexcl[threadIdx.x-1] = *((unsigned int*)(exclusions + me * exclusions_pitch) + (threadIdx.x-1));
    }

    //  __syncthreads(); ??? -Mel-
    while ((me+offset) < na) {
    //find my neighbours
    pos[threadIdx.x][0] = *((double*)((char*)positions + (me+offset) * positions_pitch) + 0);
    pos[threadIdx.x][1] = *((double*)((char*)positions + (me+offset) * positions_pitch) + 1);
    pos[threadIdx.x][2] = *((double*)((char*)positions + (me+offset) * positions_pitch) + 2);
    atomicAdd(&pos_counter, 1);
    // sync us
    __syncthreads();

    if (threadIdx.x < pos_counter) {
        dis_x = mypos[0] - pos[threadIdx.x][0];
        dis_x -= box_00 * rint(dis_x * (1.0 / box_00));
        dis_y = mypos[1] - pos[threadIdx.x][1];
        dis_y -= box_11 * rint(dis_y * (1.0 / box_11));
        dis_z = mypos[2] - pos[threadIdx.x][2];
        dis_z -= box_22 * rint(dis_z * (1.0 / box_22));

        double dis_2 = (dis_x * dis_x) + (dis_y * dis_y) + (dis_z * dis_z);
 
            //im a solute atom
            //compare dis_2 with cutoff
            if (dis_2 < cs2) {
                // check exclusion
                isex = false;
                for (unsigned int i = 0; i < num_excl; i++) {
                    if ((me+offset) == myexcl[i]) {
                        isex = true;
                    }
                }
                if (!isex) {
                    myres = atomicAdd(&res_counter_short, 1);
                    res_short[myres] = me+offset;    
                }
                isex = false; 
            } else if (dis_2 > cs2 && dis_2 < cl2) {
                myres = atomicAdd(&res_counter_long, 1);
                res_long[myres] = me+offset;
                
            } 
            //done comparing
            // get number of elements in pairlist
            if (threadIdx.x == 0) {
                elements_short = *((unsigned int*)((char*)solute_short_elements + me * solute_short_elements_pitch));
                elements_long = *((unsigned int*)((char*)solute_long_elements + me * solute_long_elements_pitch));                
            }
            // sync us
            __syncthreads();

         
       // write to memory
            if (threadIdx.x < res_counter_short) {
                *((unsigned int*)((char*)solute_short + me * solute_short_pitch) + (elements_short+threadIdx.x)) = res_short[threadIdx.x];
            }
            if (threadIdx.x < res_counter_long) {
                *((unsigned int*)((char*)solute_long + me * solute_long_pitch) + (elements_long+threadIdx.x)) = res_long[threadIdx.x];
            }
            // sync us
            __syncthreads();
          } // if thread < pos
            // add counter
            if (threadIdx.x == 0) {
                *((unsigned int*)((char*)solute_short_elements + me * solute_short_elements_pitch)) = elements_short + res_counter_short;
                *((unsigned int*)((char*)solute_long_elements + me * solute_long_elements_pitch)) = elements_long + res_counter_long ;

            // reset
            res_counter_short = 0;
            res_counter_long = 0;
            pos_counter = 0;
            elements_short = 0;
            elements_long = 0;
            }
            // done
            // sync us
            __syncthreads();
            // done

    // next block
    offset = offset+128;
    }// while have atoms
  } // me nsol
}

__global__ void calc_pairlist_solvent () {

   //calculate pairlist
    //get my number (solvent)
    unsigned long me = (blockIdx.x + (blockIdx.y*gridDim.x)) + nsol;
    if (me < na) {
    // cache
    __shared__ double pos[128][3];
    __shared__ unsigned int res_short[128];
    __shared__ unsigned int res_long[128];
    __shared__ int res_counter_short;
    __shared__ int res_counter_long;
    __shared__ double mypos[3];
    __shared__ unsigned int pos_counter;
    __shared__ unsigned int elements_short ;
    __shared__ unsigned int elements_long ;

    // distances
    double dis_x;
    double dis_y;
    double dis_z;

    // local result counter
    unsigned int myres;
    // loop konstant
    unsigned int offset;
    if ((na-me)%3 == 0) {
        // ox
        offset = threadIdx.x + 3;
    } else if ((na-me)%3 == 2) {
        // first H
         offset = threadIdx.x + 2;
    } else if ((na-me)%3 == 1) {
        // second H
        offset = threadIdx.x + 1;
    }

    //save my position
    if (threadIdx.x == 0) {
    // get positions once
    mypos[0] = *((double*)((char*)positions + me * positions_pitch) + 0);
    mypos[1] = *((double*)((char*)positions + me * positions_pitch) + 1);
    mypos[2] = *((double*)((char*)positions + me * positions_pitch) + 2);

    // initialize
    res_counter_short = 0;
    res_counter_long = 0;
    pos_counter = 0;
    elements_short = 0;
    elements_long = 0;
    }

    while ((me+offset) < na) {
    //find my neighbours
    pos[threadIdx.x][0] = *((double*)((char*)positions + (me+offset) * positions_pitch) + 0);
    pos[threadIdx.x][1] = *((double*)((char*)positions + (me+offset) * positions_pitch) + 1);
    pos[threadIdx.x][2] = *((double*)((char*)positions + (me+offset) * positions_pitch) + 2);
    atomicAdd(&pos_counter, 1);
    // sync us
    __syncthreads();

    if (threadIdx.x < pos_counter) {
        dis_x = mypos[0] - pos[threadIdx.x][0];
        dis_x -= box_00 * rint(dis_x * (1.0 / box_00));
        dis_y = mypos[1] - pos[threadIdx.x][1];
        dis_y -= box_11 * rint(dis_y * (1.0 / box_11));
        dis_z = mypos[2] - pos[threadIdx.x][2];
        dis_z -= box_22 * rint(dis_z * (1.0 / box_22));

        double dis_2 = (dis_x * dis_x) + (dis_y * dis_y) + (dis_z * dis_z);

            // im a solvent atom
            //compare dis_2 with cutoff
            if (dis_2 < cs2) {
                myres = atomicAdd(&res_counter_short, 1);
                res_short[myres] = me+offset;
            } else if (dis_2 > cs2 && dis_2 < cl2) {
                myres = atomicAdd(&res_counter_long, 1);
                res_long[myres] = me+offset;
            }
            //done comparing
            // get number of elements in pairlist
            if (threadIdx.x == 0) {
                elements_short = *((unsigned int*)((char*)solvent_short_elements + me * solvent_short_elements_pitch));
                elements_long = *((unsigned int*)((char*)solvent_long_elements + me * solvent_long_elements_pitch));
            }
            // sync us
            __syncthreads();
            // write to memory
            if (threadIdx.x < res_counter_short) {
                 *((unsigned int*)((char*)solvent_short + me * solvent_short_pitch) + (elements_short+threadIdx.x)) = res_short[threadIdx.x];
            }
            if (threadIdx.x < res_counter_long) {
                 *((unsigned int*)((char*)solvent_long + me * solvent_long_pitch) + (elements_long+threadIdx.x)) = res_long[threadIdx.x];
            }
             // sync us
             __syncthreads();
            } // if thread < pos
            // add counter
            if (threadIdx.x == 0) {
                *((unsigned int*)((char*)solvent_short_elements + me * solvent_short_elements_pitch)) = elements_short + res_counter_short;
                *((unsigned int*)((char*)solvent_long_elements + me * solvent_long_elements_pitch)) = elements_long + res_counter_long ;

            // reset
            res_counter_short = 0;
            res_counter_long = 0;
            pos_counter = 0;
            elements_short = 0;
            elements_long = 0;
            }
            // done
            // sync us
            __syncthreads();

   
    // next block
    offset = offset+128;   
    }// while have atoms
  }// me na
}

/*
//DEBUG
__global__ void tester() {

   unsigned int me = threadIdx.x+(blockDim.x*blockIdx.x);

    //DEBUG ////////////////////////////////////
    if (me == 0) {
    // box
    printf ("box_00: %f \n", box_00);
    printf ("box_01: %f \n", box_01);
    printf ("box_02: %f \n", box_02);
    printf ("box_10: %f \n", box_10);
    printf ("box_11: %f \n", box_11);
    printf ("box_12: %f \n", box_12);
    printf ("box_20: %f \n", box_20);
    printf ("box_21: %f \n", box_21);
    printf ("box_22: %f \n", box_22);
    // cutoff
    printf ("cs: %f \n", cs);
    printf ("cs2: %f \n", cs2);
    printf ("cl: %f \n", cl);
    printf ("cl2: %f \n", cl2);
    //atom numbers
    printf ("na: %i \n", na);
    printf ("nc: %i \n", nc);
    printf ("nsc: %i \n", nsc);
    //positions
    printf ("positions pointer: %p \n", positions);
    printf ("position pitch: %i \n", positions_pitch);
    double* testposition;
    testposition = (double*)((char*)positions + me * positions_pitch) + 0;          
    printf ("first position: %f \n", *testposition);
    //cogs
    printf ("cogs pointer: %p \n", cogs);
    printf ("cog pitch: %i \n", cogs_pitch);
    double* testcogs;
    testcogs = (double*)((char*)cogs + me * cogs_pitch) + 0;
    printf ("first cog: %f \n", *testcogs);
    //chargegroups
    printf ("chargegroup pointer: %p \n", chargegroups);
    printf ("chargegroup pitch: %i \n", chargegroups_pitch);
    unsigned int* testcg;
    testcg = (unsigned int*)((char*)chargegroups + me * chargegroups_pitch) + 0;
    printf ("first chargegroup: %i \n", *testcg);
    //exclusions
    printf ("exclusions pointer: %p \n", exclusions);
    printf ("exclusions pitch: %i \n", exclusions_pitch);
    printf ("exclusions elements pointer: %p \n", exclusions_elements);
    printf ("exclusions elements pitch: %i \n", exclusions_elements_pitch);
    unsigned int* testex;
    unsigned int* testexe;
    testex = (unsigned int*)((char*)exclusions + me * exclusions_pitch) + 0;
    testexe = (unsigned int*)((char*)exclusions_elements + me * exclusions_elements_pitch) + 0;
    printf ("first exclusions: %i \n", *testex);
    printf ("first exclusions element: %i \n", *testexe);
    //results
        //solute short
        printf ("solute short pointer: %p \n", solute_short);
        printf ("solute short pitch: %i \n", solute_short_pitch);
        printf ("solute short elements pointer: %p \n", solute_short_elements);
        printf ("solute short elements pitch: %i \n", solute_short_elements_pitch);
        //solute long
        printf ("solute long pointer: %p \n", solute_long);
        printf ("solute long pitch: %i \n", solute_long_pitch);
        printf ("solute long elements pointer: %p \n", solute_long_elements);
        printf ("solute long elements pitch: %i \n", solute_long_elements_pitch);
        //solvent short
        printf ("solvent short pointer: %p \n", solvent_short);
        printf ("solvent short pitch: %i \n", solvent_short_pitch);
        printf ("solvent short elements pointer: %p \n", solvent_short_elements);
        printf ("solvent short elements pitch: %i \n", solvent_short_elements_pitch);
        //solute long
        printf ("solvent long pointer: %p \n", solvent_long);
        printf ("solvent long pitch: %i \n", solvent_long_pitch);
        printf ("solvent long elements pointer: %p \n", solvent_long_elements);
        printf ("solvent long elements pitch: %i \n", solvent_long_elements_pitch);
   }
    ////////////////////////////////////
    //printf("%s \n", "TESTER IN ACTION");
   // printf("ME: %i \n", me);
    

}
// */

__global__ void clear_results() {

    unsigned int me = threadIdx.x+(blockDim.x*blockIdx.x);

    if (me < na) {

        if (me < nsol) {
            unsigned int* re_ss = (unsigned int*)((char*)solute_short_elements + me * solute_short_elements_pitch) + 0;
            unsigned int* re_sl = (unsigned int*)((char*)solute_long_elements + me * solute_long_elements_pitch) + 0;
            *re_ss = 0;
            *re_sl = 0;
        }
        unsigned int* re_sls = (unsigned int*)((char*)solvent_short_elements + me * solvent_short_elements_pitch) + 0;
        unsigned int* re_sll = (unsigned int*)((char*)solvent_long_elements + me * solvent_long_elements_pitch) + 0;
       *re_sls = 0;
       *re_sll = 0;
    }
 }
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

//DEBUG FUNCTION
/*
void debug_memmap(gcuda::memmap &gpumemory) {
    
    unsigned int map_counter = gpumemory.counter;
    unsigned int map_size = gpumemory.size;

    std::cout << "Current counter: " << map_counter << std::endl;
     std::cout << "Size: " << map_size << std::endl;
     
    std::map<std::string, unsigned int>::iterator map_start;
    std::cout << "Memory Map: " << std::endl;
    for (map_start = gpumemory.map.begin(); map_start != gpumemory.map.end(); map_start++) {
        std::cout << "\tKey: " << (*map_start).first << " Value: " << (*map_start).second << std::endl;
    }
    std::cout << "Memory layout: " << std::endl;
    for (unsigned int c=0; c < map_counter; c++) {
        std::cout << "\tCounter: " << c << " Pointer: " <<  gpumemory.memory[c] << " Pitch: " << gpumemory.pitches[c] << std::endl;
    }

}*/
//functions
extern "C" void update_atomic(gcuda::memmap &gpumemory, double* h_positions, math::Box box, unsigned int nratoms, unsigned int nr_solute_atoms,
                          gcuda::TwoDArray<unsigned int> *h_results_solute_short,
                          gcuda::TwoDArray<unsigned int> *h_results_solute_long,
                          gcuda::TwoDArray<unsigned int> *h_results_solvent_short,
                          gcuda::TwoDArray<unsigned int> *h_results_solvent_long) {

    //copy box
    bool cp = copybox(box);

    //copy positions
    cudaMemcpy2D((void*)gpumemory.memory[gpumemory.map["positions"]], gpumemory.pitches[gpumemory.map["positions"]], h_positions, 24, 24, nratoms, cudaMemcpyHostToDevice);


    // kernel calls
   
    int threadsPerBlock = 128;
    int numBlocks = (nratoms / threadsPerBlock) + 1;
    dim3 numBlocksPsolute;
    numBlocksPsolute.x = sqrt(nr_solute_atoms) + 1;
    numBlocksPsolute.y = sqrt(nr_solute_atoms) + 1;
    numBlocksPsolute.z = 1;

    dim3 numBlocksPsolvent;
    numBlocksPsolvent.x = sqrt((nratoms - nr_solute_atoms)) + 1;
    numBlocksPsolvent.y = sqrt((nratoms - nr_solute_atoms)) + 1;
    numBlocksPsolvent.z = 1;

   /* int numBlocksPsolute = nr_solute_atoms;
    int numBlocksPsolvent = nratoms - nr_solute_atoms; */
     // clear
    clear_results<<<numBlocks, threadsPerBlock>>>();

    //calc
    calc_pairlist_solute<<<numBlocksPsolute, threadsPerBlock>>>();
    calc_pairlist_solvent<<<numBlocksPsolvent, threadsPerBlock>>>();
    //copy back results
    size_t epitch = 4;

     cudaMemcpy2D (*h_results_solute_short->ptr, h_results_solute_short->pitch, (void*)gpumemory.memory[gpumemory.map["solute_short"]], gpumemory.pitches[gpumemory.map["solute_short"]], h_results_solute_short->pitch, nr_solute_atoms, cudaMemcpyDeviceToHost);
     cudaMemcpy2D (h_results_solute_short->elements, epitch, (void*)gpumemory.memory[gpumemory.map["solute_short_elements"]], gpumemory.pitches[gpumemory.map["solute_short_elements"]], epitch, nr_solute_atoms, cudaMemcpyDeviceToHost);
     cudaMemcpy2D (*h_results_solute_long->ptr, h_results_solute_long->pitch, (void*)gpumemory.memory[gpumemory.map["solute_long"]], gpumemory.pitches[gpumemory.map["solute_long"]], h_results_solute_long->pitch, nr_solute_atoms, cudaMemcpyDeviceToHost);
     cudaMemcpy2D (h_results_solute_long->elements, epitch, (void*)gpumemory.memory[gpumemory.map["solute_long_elements"]], gpumemory.pitches[gpumemory.map["solute_long_elements"]], epitch, nr_solute_atoms, cudaMemcpyDeviceToHost);
     cudaMemcpy2D (*h_results_solvent_short->ptr, h_results_solvent_short->pitch, (void*)gpumemory.memory[gpumemory.map["solvent_short"]], gpumemory.pitches[gpumemory.map["solvent_short"]], h_results_solvent_short->pitch, nratoms, cudaMemcpyDeviceToHost);
     cudaMemcpy2D (h_results_solvent_short->elements, epitch, (void*)gpumemory.memory[gpumemory.map["solvent_short_elements"]], gpumemory.pitches[gpumemory.map["solvent_short_elements"]], epitch, nratoms, cudaMemcpyDeviceToHost);
     cudaMemcpy2D (*h_results_solvent_long->ptr, h_results_solvent_long->pitch, (void*)gpumemory.memory[gpumemory.map["solvent_long"]], gpumemory.pitches[gpumemory.map["solvent_long"]], h_results_solvent_long->pitch, nratoms, cudaMemcpyDeviceToHost);
     cudaMemcpy2D (h_results_solvent_long->elements, epitch, (void*)gpumemory.memory[gpumemory.map["solvent_long_elements"]], gpumemory.pitches[gpumemory.map["solvent_long_elements"]], epitch, nratoms, cudaMemcpyDeviceToHost);



     //TEST
  /*          cudaFuncAttributes details;
            cudaError_t e6;
            std::cout << "############ CLEAR RESULTS ATOMIC ##################" << std::endl;
            e6 = cudaFuncGetAttributes(&details, clear_results);
            std::cout << "e6 " << e6 << std::endl;
            std::cout << "Compute capability " << details.binaryVersion << std::endl;
            std::cout << "Size of constant memory " << details.constSizeBytes / 1024 << " kb" << std::endl;
            std::cout << "Size of local memory " << details.localSizeBytes / 1024 << " kb" << std::endl;
            std::cout << "Max num threads per block " << details.maxThreadsPerBlock << std::endl;
            std::cout << "Number registers  " << details.numRegs << std::endl;
            std::cout << "Shared memory " << details.sharedSizeBytes / 1024 << " kb" << std::endl;
            std::cout << "####################################################" << std::endl;

            std::cout << "############ CALC PAIRLIST ATOMIC ##################" << std::endl;
            e6 = cudaFuncGetAttributes(&details, calc_pairlist_solute);
            std::cout << "e6 " << e6 << std::endl;
            std::cout << "Compute capability " << details.binaryVersion << std::endl;
            std::cout << "Size of constant memory " << details.constSizeBytes / 1024 << " kb" << std::endl;
            std::cout << "Size of local memory " << details.localSizeBytes / 1024 << " kb" << std::endl;
            std::cout << "Max num threads per block " << details.maxThreadsPerBlock << std::endl;
            std::cout << "Number registers  " << details.numRegs << std::endl;
            std::cout << "Shared memory " << details.sharedSizeBytes / 1024 << " kb" << std::endl;
            std::cout << "####################################################" << std::endl;
            
             std::cout << "############ CALC PAIRLIST ATOMIC ##################" << std::endl;
            e6 = cudaFuncGetAttributes(&details, calc_pairlist_solvent);
            std::cout << "e6 " << e6 << std::endl;
            std::cout << "Compute capability " << details.binaryVersion << std::endl;
            std::cout << "Size of constant memory " << details.constSizeBytes / 1024 << " kb" << std::endl;
            std::cout << "Size of local memory " << details.localSizeBytes / 1024 << " kb" << std::endl;
            std::cout << "Max num threads per block " << details.maxThreadsPerBlock << std::endl;
            std::cout << "Number registers  " << details.numRegs << std::endl;
            std::cout << "Shared memory " << details.sharedSizeBytes / 1024 << " kb" << std::endl;
            std::cout << "####################################################" << std::endl;
     //
*/
      return ;
}

extern "C" void update_cg(gcuda::memmap &gpumemory, double* h_positions, math::Box box, unsigned int nratoms, unsigned int nr_solute_atoms, unsigned int nrcg,
                          gcuda::TwoDArray<unsigned int> *h_results_solute_short,
                          gcuda::TwoDArray<unsigned int> *h_results_solute_long,
                          gcuda::TwoDArray<unsigned int> *h_results_solvent_short,
                          gcuda::TwoDArray<unsigned int> *h_results_solvent_long) {

    //copy constants
    bool cp = copybox(box);

    //pos

    //DEBUG
    /*double* d_positions;
    size_t ppitch;
    cudaMallocPitch((void**)&d_positions, &ppitch, 24, nratoms);
      cudaMemcpyToSymbol("positions", &d_positions , 8, 0, cudaMemcpyHostToDevice);
     cudaMemcpyToSymbol("positions_pitch", &ppitch , sizeof(size_t), 0, cudaMemcpyHostToDevice);
     std::cout << "TEST position pointer => " << d_positions << std::endl;
      std::cout << "TEST position pitch => " << ppitch << std::endl;
      std::cout << "TEST position host => " << *h_positions << std::endl;
     // */
    cudaMemcpy2D((void*)gpumemory.memory[gpumemory.map["positions"]], gpumemory.pitches[gpumemory.map["positions"]], h_positions, 24, 24, nratoms, cudaMemcpyHostToDevice);
    

   // debug_memmap(gpumemory);

     // kernel calls

    int threadsPerBlock = 128;
    int numBlocks = (nratoms / threadsPerBlock) + 1;

    dim3 numBlocksPsolute;
    numBlocksPsolute.x = sqrt(nrcg) + 1;
    numBlocksPsolute.y = sqrt(nrcg) + 1;
    numBlocksPsolute.z = 1;

    dim3 numBlocksPsolvent;
    numBlocksPsolvent.x = sqrt(nrcg) + 1;
    numBlocksPsolvent.y = sqrt(nrcg) + 1;
    numBlocksPsolvent.z = 1;
    //DEBUG
    //tester<<<numBlocks, threadsPerBlock>>>();

    // clear
    clear_results<<<numBlocks, threadsPerBlock>>>();

    // calculate cogs
    calc_cogs<<<numBlocks, threadsPerBlock>>>();

    // calculate pairlist
    calc_pairlist_cg_self<<<numBlocks, threadsPerBlock>>>();
    calc_cg_pairlist_solute<<<numBlocksPsolute, threadsPerBlock>>>();
    calc_cg_pairlist_solvent<<<numBlocksPsolvent, threadsPerBlock>>>();
    //calc_cg_pairlist<<<numBlocks, threadsPerBlock>>>();
    


    //copy back results
    size_t epitch = 4;

     cudaMemcpy2D (*h_results_solute_short->ptr, h_results_solute_short->pitch, (void*)gpumemory.memory[gpumemory.map["solute_short"]], gpumemory.pitches[gpumemory.map["solute_short"]], h_results_solute_short->pitch, nr_solute_atoms, cudaMemcpyDeviceToHost);
     cudaMemcpy2D (h_results_solute_short->elements, epitch, (void*)gpumemory.memory[gpumemory.map["solute_short_elements"]], gpumemory.pitches[gpumemory.map["solute_short_elements"]], epitch, nr_solute_atoms, cudaMemcpyDeviceToHost);
     cudaMemcpy2D (*h_results_solute_long->ptr, h_results_solute_long->pitch, (void*)gpumemory.memory[gpumemory.map["solute_long"]], gpumemory.pitches[gpumemory.map["solute_long"]], h_results_solute_long->pitch, nr_solute_atoms, cudaMemcpyDeviceToHost);
     cudaMemcpy2D (h_results_solute_long->elements, epitch, (void*)gpumemory.memory[gpumemory.map["solute_long_elements"]], gpumemory.pitches[gpumemory.map["solute_long_elements"]], epitch, nr_solute_atoms, cudaMemcpyDeviceToHost);
     cudaMemcpy2D (*h_results_solvent_short->ptr, h_results_solvent_short->pitch, (void*)gpumemory.memory[gpumemory.map["solvent_short"]], gpumemory.pitches[gpumemory.map["solvent_short"]], h_results_solvent_short->pitch, nratoms, cudaMemcpyDeviceToHost);
     cudaMemcpy2D (h_results_solvent_short->elements, epitch, (void*)gpumemory.memory[gpumemory.map["solvent_short_elements"]], gpumemory.pitches[gpumemory.map["solvent_short_elements"]], epitch, nratoms, cudaMemcpyDeviceToHost);
     cudaMemcpy2D (*h_results_solvent_long->ptr, h_results_solvent_long->pitch, (void*)gpumemory.memory[gpumemory.map["solvent_long"]], gpumemory.pitches[gpumemory.map["solvent_long"]], h_results_solvent_long->pitch, nratoms, cudaMemcpyDeviceToHost);
     cudaMemcpy2D (h_results_solvent_long->elements, epitch, (void*)gpumemory.memory[gpumemory.map["solvent_long_elements"]], gpumemory.pitches[gpumemory.map["solvent_long_elements"]], epitch, nratoms, cudaMemcpyDeviceToHost);


      /*    //TEST
            cudaFuncAttributes details;
            cudaError_t e6;
            std::cout << "############ CLEAR RESULTS CG ##################" << std::endl;
            e6 = cudaFuncGetAttributes(&details, clear_results);
            std::cout << "e6 " << e6 << std::endl;
            std::cout << "Compute capability " << details.binaryVersion << std::endl;
            std::cout << "Size of constant memory " << details.constSizeBytes / 1024 << " kb" << std::endl;
            std::cout << "Size of local memory " << details.localSizeBytes / 1024 << " kb" << std::endl;
            std::cout << "Max num threads per block " << details.maxThreadsPerBlock << std::endl;
            std::cout << "Number registers  " << details.numRegs << std::endl;
            std::cout << "Shared memory " << details.sharedSizeBytes / 1024 << " kb" << std::endl;
            std::cout << "####################################################" << std::endl;

            std::cout << "############ CALC COGS CG ##################" << std::endl;
            e6 = cudaFuncGetAttributes(&details, calc_cogs);
            std::cout << "e6 " << e6 << std::endl;
            std::cout << "Compute capability " << details.binaryVersion << std::endl;
            std::cout << "Size of constant memory " << details.constSizeBytes / 1024 << " kb" << std::endl;
            std::cout << "Size of local memory " << details.localSizeBytes / 1024 << " kb" << std::endl;
            std::cout << "Max num threads per block " << details.maxThreadsPerBlock << std::endl;
            std::cout << "Number registers  " << details.numRegs << std::endl;
            std::cout << "Shared memory " << details.sharedSizeBytes / 1024 << " kb" << std::endl;
            std::cout << "####################################################" << std::endl;

            std::cout << "############ CALC PAIRLIST CG ##################" << std::endl;
            e6 = cudaFuncGetAttributes(&details, calc_pairlist_cg_self);
            std::cout << "e6 " << e6 << std::endl;
            std::cout << "Compute capability " << details.binaryVersion << std::endl;
            std::cout << "Size of constant memory " << details.constSizeBytes / 1024 << " kb" << std::endl;
            std::cout << "Size of local memory " << details.localSizeBytes / 1024 << " kb" << std::endl;
            std::cout << "Max num threads per block " << details.maxThreadsPerBlock << std::endl;
            std::cout << "Number registers  " << details.numRegs << std::endl;
            std::cout << "Shared memory " << details.sharedSizeBytes / 1024 << " kb" << std::endl;
            std::cout << "####################################################" << std::endl;
            std::cout << "############ CALC PAIRLIST CG ##################" << std::endl;
            e6 = cudaFuncGetAttributes(&details, calc_cg_pairlist_solute);
            std::cout << "e6 " << e6 << std::endl;
            std::cout << "Compute capability " << details.binaryVersion << std::endl;
            std::cout << "Size of constant memory " << details.constSizeBytes / 1024 << " kb" << std::endl;
            std::cout << "Size of local memory " << details.localSizeBytes / 1024 << " kb" << std::endl;
            std::cout << "Max num threads per block " << details.maxThreadsPerBlock << std::endl;
            std::cout << "Number registers  " << details.numRegs << std::endl;
            std::cout << "Shared memory " << details.sharedSizeBytes / 1024 << " kb" << std::endl;
            std::cout << "####################################################" << std::endl;
            std::cout << "############ CALC PAIRLIST CG ##################" << std::endl;
            e6 = cudaFuncGetAttributes(&details, calc_cg_pairlist_solvent);
            std::cout << "e6 " << e6 << std::endl;
            std::cout << "Compute capability " << details.binaryVersion << std::endl;
            std::cout << "Size of constant memory " << details.constSizeBytes / 1024 << " kb" << std::endl;
            std::cout << "Size of local memory " << details.localSizeBytes / 1024 << " kb" << std::endl;
            std::cout << "Max num threads per block " << details.maxThreadsPerBlock << std::endl;
            std::cout << "Number registers  " << details.numRegs << std::endl;
            std::cout << "Shared memory " << details.sharedSizeBytes / 1024 << " kb" << std::endl;
            std::cout << "####################################################" << std::endl;
     // */

          return ;
}

// initialise the GPU memory - chargegroup
extern "C" void initcuda_cg (gcuda::memmap &gpumemory, unsigned int nratoms, unsigned int nr_solute_atoms, unsigned int nrcg, unsigned int nrsolutecg, double *cutoff, size_t host_results_pitch) {

    //definition of variables
    double *d_positions;
    unsigned int *d_chargegroups;
    unsigned int *d_exclusions;
    unsigned int *d_e_exclusions;

    // device only
    double *d_cogs;

    //results
    unsigned int *d_results_solute_short;
    unsigned int *d_results_solute_long;
    unsigned int *d_results_solvent_short;
    unsigned int *d_results_solvent_long;
    // results element counter array
    unsigned int *d_e_results_solute_short;
    unsigned int *d_e_results_solute_long;
    unsigned int *d_e_results_solvent_short;
    unsigned int *d_e_results_solvent_long;

    //pitches
    //pos
    size_t  ppitch;
    //results
    size_t  r_pitch_solute_short;
    size_t  r_e_pitch_solute_short;
    size_t  r_pitch_solute_long;
    size_t  r_e_pitch_solute_long;
    size_t  r_pitch_solvent_short;
    size_t  r_e_pitch_solvent_short;
    size_t  r_pitch_solvent_long;
    size_t  r_e_pitch_solvent_long;

    //other
    size_t  epitch = 4;  // elements array pitch on the host
    size_t  cgpitch;
    size_t  cogpitch;
    size_t  expitch;
    size_t  ex_e_pitch;
   
    // allocate memory
    //pos

    cudaMallocPitch((void**)&d_positions, &ppitch, 24, nratoms);
    //cgs
    cudaMallocPitch((void**)&d_chargegroups, &cgpitch, 8, nrcg);
    //exclusions
    cudaMallocPitch((void**)&d_exclusions, &expitch, gpumemory.pitches[gpumemory.map["host_exclusions"]] , nr_solute_atoms);
    cudaMallocPitch((void**)&d_e_exclusions, &ex_e_pitch, 4, nr_solute_atoms);
    //cogs
    cudaMallocPitch((void**)&d_cogs, &cogpitch, 24, nrcg);
    //results
    cudaMallocPitch((void**)&d_results_solute_short, &r_pitch_solute_short, host_results_pitch, nr_solute_atoms);
    cudaMallocPitch((void**)&d_e_results_solute_short, &r_e_pitch_solute_short, epitch, nr_solute_atoms);
    cudaMallocPitch((void**)&d_results_solute_long, &r_pitch_solute_long, host_results_pitch, nr_solute_atoms);
    cudaMallocPitch((void**)&d_e_results_solute_long, &r_e_pitch_solute_long, epitch, nr_solute_atoms);
    cudaMallocPitch((void**)&d_results_solvent_short, &r_pitch_solvent_short, host_results_pitch, nratoms);
    cudaMallocPitch((void**)&d_e_results_solvent_short, &r_e_pitch_solvent_short, epitch, nratoms);
    cudaMallocPitch((void**)&d_results_solvent_long, &r_pitch_solvent_long, host_results_pitch, nratoms);
    cudaMallocPitch((void**)&d_e_results_solvent_long, &r_e_pitch_solvent_long, epitch, nratoms);

    // copy one time stuff
    cudaMemcpy2D(d_chargegroups, cgpitch, (void*)gpumemory.memory[gpumemory.map["host_chargegroups"]], gpumemory.pitches[gpumemory.map["host_chargegroups"]], 8, nrcg, cudaMemcpyHostToDevice);
    cudaMemcpy2D(d_exclusions, expitch, (void*)gpumemory.memory[gpumemory.map["host_exclusions"]], gpumemory.pitches[gpumemory.map["host_exclusions"]],gpumemory.pitches[gpumemory.map["host_exclusions"]], nr_solute_atoms, cudaMemcpyHostToDevice);
    cudaMemcpy2D(d_e_exclusions, ex_e_pitch, (void*)gpumemory.memory[gpumemory.map["host_exclusions_elements"]], epitch, 4, nr_solute_atoms, cudaMemcpyHostToDevice);
    
    // draw map
    gpumemory.add("positions", d_positions, ppitch);
    gpumemory.add("chargegroups", d_chargegroups, cgpitch);
    gpumemory.add("exclusions", d_exclusions, expitch);
    gpumemory.add("exclusions_elements", d_e_exclusions, ex_e_pitch);
    gpumemory.add("cogs", d_cogs, cogpitch);
    gpumemory.add("solute_short", d_results_solute_short, r_pitch_solute_short);
    gpumemory.add("solute_short_elements", d_e_results_solute_short, r_e_pitch_solute_short);
    gpumemory.add("solute_long", d_results_solute_long, r_pitch_solute_long);
    gpumemory.add("solute_long_elements", d_e_results_solute_long, r_e_pitch_solute_long);
    gpumemory.add("solvent_short", d_results_solvent_short, r_pitch_solvent_short);
    gpumemory.add("solvent_short_elements", d_e_results_solvent_short, r_e_pitch_solvent_short);
    gpumemory.add("solvent_long", d_results_solvent_long, r_pitch_solvent_long);
    gpumemory.add("solvent_long_elements", d_e_results_solvent_long, r_e_pitch_solvent_long);

    //DEBUG
    //debug_memmap(gpumemory);
    //
    //copy constants
    bool cc = copycutoff(cutoff);
    bool nr = copynratoms(nratoms);
    bool nc = copynrcg(nrcg);
    bool ns = copynrsolutecg(nrsolutecg);
    bool nsol = copynrsolute(nr_solute_atoms);
    bool memmap = copymemorymap(gpumemory, false);
    return;
}

// initialise the GPU memory - atomic
extern "C" void initcuda_atomic (gcuda::memmap &gpumemory, unsigned int nratoms, unsigned int nr_solute_atoms, double *cutoff, size_t host_results_pitch) {

    //definition of variables
    double *d_positions;
    unsigned int *d_exclusions;
    unsigned int *d_e_exclusions;

    //results
    unsigned int *d_results_solute_short;
    unsigned int *d_results_solute_long;
    unsigned int *d_results_solvent_short;
    unsigned int *d_results_solvent_long;
    // results element counter array
    unsigned int *d_e_results_solute_short;
    unsigned int *d_e_results_solute_long;
    unsigned int *d_e_results_solvent_short;
    unsigned int *d_e_results_solvent_long;

    //pitches
    //pos
    size_t  ppitch;
    //results
    size_t  r_pitch_solute_short;
    size_t  r_e_pitch_solute_short;
    size_t  r_pitch_solute_long;
    size_t  r_e_pitch_solute_long;
    size_t  r_pitch_solvent_short;
    size_t  r_e_pitch_solvent_short;
    size_t  r_pitch_solvent_long;
    size_t  r_e_pitch_solvent_long;

    //other
    size_t  epitch = 4;  // elements array pitch on the host
    size_t  expitch;
    size_t  ex_e_pitch;

    // allocate memory
    //pos
    cudaMallocPitch((void**)&d_positions, &ppitch, 24, nratoms);
    //exclusions
    cudaMallocPitch((void**)&d_exclusions, &expitch, gpumemory.pitches[gpumemory.map["host_exclusions"]] , nr_solute_atoms);
    cudaMallocPitch((void**)&d_e_exclusions, &ex_e_pitch, 4, nr_solute_atoms);
    //results
    cudaMallocPitch((void**)&d_results_solute_short, &r_pitch_solute_short, host_results_pitch, nr_solute_atoms);
    cudaMallocPitch((void**)&d_e_results_solute_short, &r_e_pitch_solute_short, epitch, nr_solute_atoms);
    cudaMallocPitch((void**)&d_results_solute_long, &r_pitch_solute_long, host_results_pitch, nr_solute_atoms);
    cudaMallocPitch((void**)&d_e_results_solute_long, &r_e_pitch_solute_long, epitch, nr_solute_atoms);
    cudaMallocPitch((void**)&d_results_solvent_short, &r_pitch_solvent_short, host_results_pitch, nratoms);
    cudaMallocPitch((void**)&d_e_results_solvent_short, &r_e_pitch_solvent_short, epitch, nratoms);
    cudaMallocPitch((void**)&d_results_solvent_long, &r_pitch_solvent_long, host_results_pitch, nratoms);
    cudaMallocPitch((void**)&d_e_results_solvent_long, &r_e_pitch_solvent_long, epitch, nratoms);

    // copy one time stuff
    cudaMemcpy2D(d_exclusions, expitch, (void*)gpumemory.memory[gpumemory.map["host_exclusions"]], gpumemory.pitches[gpumemory.map["host_exclusions"]],gpumemory.pitches[gpumemory.map["host_exclusions"]], nr_solute_atoms, cudaMemcpyHostToDevice);
    cudaMemcpy2D(d_e_exclusions, ex_e_pitch, (void*)gpumemory.memory[gpumemory.map["host_exclusions_elements"]], epitch, 4, nr_solute_atoms, cudaMemcpyHostToDevice);

    // draw map
    gpumemory.add("positions", d_positions, ppitch);
    gpumemory.add("exclusions", d_exclusions, expitch);
    gpumemory.add("exclusions_elements", d_e_exclusions, ex_e_pitch);
    gpumemory.add("solute_short", d_results_solute_short, r_pitch_solute_short);
    gpumemory.add("solute_short_elements", d_e_results_solute_short, r_e_pitch_solute_short);
    gpumemory.add("solute_long", d_results_solute_long, r_pitch_solute_long);
    gpumemory.add("solute_long_elements", d_e_results_solute_long, r_e_pitch_solute_long);
    gpumemory.add("solvent_short", d_results_solvent_short, r_pitch_solvent_short);
    gpumemory.add("solvent_short_elements", d_e_results_solvent_short, r_e_pitch_solvent_short);
    gpumemory.add("solvent_long", d_results_solvent_long, r_pitch_solvent_long);
    gpumemory.add("solvent_long_elements", d_e_results_solvent_long, r_e_pitch_solvent_long);


    //copy constants
    bool cc = copycutoff(cutoff);
    bool nr = copynratoms(nratoms);
    bool nsol = copynrsolute(nr_solute_atoms);
    bool memmap = copymemorymap(gpumemory, true);


    return;
}

extern "C" void freecuda (gcuda::memmap &gpumemory, bool atomic) {

    // positions
    cudaFree((void*)gpumemory.memory[gpumemory.map["positions"]]);

    // chargegroups
    if (!atomic) {
    cudaFree((void*)gpumemory.memory[gpumemory.map["chargegroups"]]);
    cudaFree((void*)gpumemory.memory[gpumemory.map["cogs"]]);
    }

    // exclusions
    cudaFree((void*)gpumemory.memory[gpumemory.map["exclusions"]]);
    cudaFree((void*)gpumemory.memory[gpumemory.map["exclusions_elements"]]);

    // resulsts
    cudaFree((void*)gpumemory.memory[gpumemory.map["solute_short"]]);
    cudaFree((void*)gpumemory.memory[gpumemory.map["solute_short_elements"]]);
    cudaFree((void*)gpumemory.memory[gpumemory.map["solute_long"]]);
    cudaFree((void*)gpumemory.memory[gpumemory.map["solute_long_elements"]]);
    cudaFree((void*)gpumemory.memory[gpumemory.map["solvent_short"]]);
    cudaFree((void*)gpumemory.memory[gpumemory.map["solvent_short_elements"]]);
    cudaFree((void*)gpumemory.memory[gpumemory.map["solvent_long"]]);
    cudaFree((void*)gpumemory.memory[gpumemory.map["solvent_long_elements"]]);


    return;
}
