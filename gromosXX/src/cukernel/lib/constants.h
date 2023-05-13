/**
 * @file constants.h
 * device constant memory
 */
#ifndef INCLUDED_CUKERNEL_CONSTANTS_H
#define INCLUDED_CUKERNEL_CONSTANTS_H

/**
 * @struct cutoffs
 */
namespace device {
    __constant__ __device__ float cutoff_short;
    __constant__ __device__ float cutoff_long;
}

#endif
