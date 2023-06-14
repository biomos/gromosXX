/**
 * @file cuda_kernel.cu
 * implementation of CUDA_Kernel for the CUDA-enabled code
 */


#include "stdheader.h"
//#include "io/message.h"

#undef MODULE
#undef SUBMODULE
#define MODULE cukernel
#define SUBMODULE kernel

#include "cuda.cc"

extern "C" cukernel::CUDA_Kernel::CUDA_Kernel(){disabled();};

extern "C" cukernel::CUDA_Kernel::~CUDA_Kernel(){disabled();};

extern "C" void cukernel::CUDA_Kernel::init(topology::Topology &topo,
                        configuration::Configuration &conf,
                        simulation::Simulation &sim){disabled();};

extern "C" gpu_status * cukernel::CUDA_Kernel::cudaInitConstraints(unsigned int num_of_gpus, unsigned int gpu_id, unsigned int num_atoms,
                    unsigned int num_solvent_mol) {return nullptr;};
