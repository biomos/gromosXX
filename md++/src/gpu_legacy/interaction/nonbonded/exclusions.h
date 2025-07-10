/**
 * @file exclusions.h
 * @author Poliak (peter.poliak@boku.ac.at)
 * @brief nonbonded exclusions to be stored on device
 * @version 0.1
 * @date 2023-06-17
 * 
 * @copyright Copyright (c) 2023
 * 
 */
#pragma once

#include "container.h"

namespace cukernel {
    /**
     * Store nonbonded exclusions in a jagged array
     */
    typedef Container<unsigned> Exclusions;
}

