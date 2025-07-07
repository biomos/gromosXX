/**
 * @file cuda_manager.cu
 * implementation of CudaManager for the CUDA-enabled code
 */


#include "stdheader.h"
//#include "io/message.h"

#undef MODULE
#undef SUBMODULE
#define MODULE cuda
#define SUBMODULE manager

#include "cukernel.cc"

extern "C" cukernel::CudaManager::CudaManager(){disabled();};

extern "C" cukernel::CudaManager::~CudaManager(){disabled();};

extern "C" void cukernel::CudaManager::init(topology::Topology &topo,
                        configuration::Configuration &conf,
                        simulation::Simulation &sim){disabled();};

extern "C" gpu_status * cukernel::CudaManager::cudaInitConstraints(unsigned int num_of_gpus, unsigned int gpu_id, unsigned int num_atoms,
                    unsigned int num_solvent_mol) {return nullptr;};
