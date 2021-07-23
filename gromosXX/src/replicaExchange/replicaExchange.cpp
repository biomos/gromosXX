//
// Created by benja on 23/07/2021.
//
/*
 * This file contains general module definitions.
 *
 */

#include "../stdheader.h"

double replicaExchange_ver = 1.0;

namespace re {
#ifdef NDEBUG
    int debug_level=0;
    int replica_exchanger_debug_level=0;
    int replica_debug_level=0;
#endif
}