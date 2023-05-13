/**
 * @file cudaKernel.cu
 * implementation of CUDA kernel
 */

#include "stdheader.h"
#include "cuda_kernel.h"
#include "util/debug.h"

#include "gpu_status.h"
//#include "cudaKernel.h"
#include "lib/constants.h"

#include "algorithm/algorithm.h"
#include "topology/topology.h"
#include "simulation/simulation.h"
#include "configuration/configuration.h"

#undef MODULE
#undef SUBMODULE
#define MODULE cuda
#define SUBMODULE kernel

#include "cuda.cc"

extern "C" cudakernel::CUDA_Kernel::CUDA_Kernel()
: mytopo(nullptr),
  myconf(nullptr),
  mysim(nullptr),
  device_number(-1)
{};

extern "C" cudakernel::CUDA_Kernel::~CUDA_Kernel(){};

extern "C" void cudakernel::CUDA_Kernel::init(topology::Topology &topo,
                        configuration::Configuration &conf,
                        simulation::Simulation &sim){
    this->mytopo = &topo;
    this->myconf = &conf;
    this->mysim = &sim;

    // initialize the devices
    // for now, we only support only single GPU
    DEBUG(0, "Device number:" << this->device_number);
    this->device_number = mysim->param().cuda.device_number.at(0);
    DEBUG(0, "Device number:" << this->device_number);
    if (this->device_number != -1) {
      cudaSetDevice(this->device_number);
    } else {
      // determine the device number
      cudaGetDevice(&this->device_number);
    }
    DEBUG(0, "Device number:" << this->device_number);
    
    float cutoff_short = mysim->param().pairlist.cutoff_short;
    cudaMemcpyToSymbol(device::cutoff_short, &cutoff_short, sizeof(float));
    DEBUG(0, "CUDA copied to symbol");
    std::cout << "CUDA copied to symbol" << std::endl;
    /*



    // let's first query the device properties and print them out
    cudaGetDeviceProperties(&this->device_properties, this->device_number);
    cudaGetDeviceFlags(&this->flags);
    if (this->flags == 0 && cudaSetDeviceFlags(cudaDeviceScheduleYield) == cudaErrorSetOnActiveProcess)
    std::cerr << "Cannot set flags\n";

    std::cout << "CUDA" << std::endl;
    std::cout << "\tDeviceproperties:" << std::endl;
    std::cout << "\t\tNumber: " << device_number << std::endl;
    std::cout << "\t\tName: " << this->device_properties.name << std::endl;
    std::cout << "\t\tTotal memory: " << this->device_properties.totalGlobalMem << std::endl;
    std::cout << "\t\tShared memory per block: " << this->device_properties.sharedMemPerBlock << std::endl;
    std::cout << "\t\tTotal constant memory: " << this->device_properties.totalConstMem << std::endl;

    // Create a gpu_status structure
    gpu_status * gpu_stat;
    gpu_stat = (gpu_status*) malloc (sizeof(gpu_status));
    gpu_stat->device = device_number;*/
};

extern "C" gpu_status * cudakernel::CUDA_Kernel::cudaInitConstraints(unsigned int num_of_gpus, unsigned int gpu_id, unsigned int num_atoms,
                    unsigned int num_solvent_mol) {return nullptr;};
