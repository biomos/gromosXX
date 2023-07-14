/**
 * @file pairlist.h
 * @author Poliak (peter.poliak@boku.ac.at)
 * @brief pairlist functions
 * @version 0.1
 * @date 2023-06-17
 * 
 * @copyright Copyright (c) 2023
 * 
 */
#ifndef INCLUDED_CUPAIRLIST_H
#define INCLUDED_CUPAIRLIST_H

#include "lib/container.h"

namespace cukernel {
    /**
     * Store exclusions in a jagged array
     */
    typedef Container<unsigned> Pairlist;

    /**
     * Update the pairlist
     */
    __global__ void update_pairlist(const float3 *pos, Pairlist pairlist) {};
}
#endif

