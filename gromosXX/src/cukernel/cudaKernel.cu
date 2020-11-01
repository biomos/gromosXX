/**
 * @file cudaKernel.cu
 * contains the implementation of the CUDA solvent loop
 */

#include <iostream>
#include <cuda.h>
#include "lib/math.h"
#include "parameter.h"
#include "pairlist.h"
#include "interaction.h"
#include "constraints.h"
#include "lib/utils.h"
#include "cudaShake.h"

using namespace cudakernel;

#include "gpu_status.h"

extern "C" void test() {
}

#include "io.cu"
#include "pairlist.cu"
#include "interaction.cu"
#include "constraints.cu"
#include "cudaShake.cu"


