/**
 * @file pairlist.h
 * @author Poliak (peter.poliak@boku.ac.at)
 * @brief pairlist storage in a jagged array, using cukernel::Container
 * @version 0.1
 * @date 2023-06-17
 * 
 * @copyright Copyright (c) 2023
 * 
 */

#pragma once

#include "lib/container.h"

namespace cukernel {
    /**
     * Store pairlist in a jagged array
     */
    typedef Container<unsigned> Pairlist;

    struct PairlistContainer:
        /**
         * shortrange pairlists that holds: solute-, solute-solute, solute-solvent pairs
         */
        Pairlist solute_short;
        /**
         * longrange pairlists that holds: solute-, solute-solute, solute-solvent pairs
         */
        Pairlist solute_long;
        /**
         * shortrange pairlists that holds: solvent-solvent pairs
         */
        Pairlist solvent_short;
        /**
         * longrange pairlists that holds: solvent-solvent pairs
         */
        Pairlist solvent_long; 


    /**
     * Update the pairlist
     */
    __global__ void update_pairlist(const float3 *pos, Pairlist pairlist) {};
}

