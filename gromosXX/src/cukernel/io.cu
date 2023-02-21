/**
 * @file io.cu
 * I/O, initalization and memory management
 */
#include <iostream>

#include "parameter.h"
#include "gpu_status.h"
#include "lib/utils.h"


extern "C" gpu_status * cudaInit(int & device_number,
            unsigned int num_atoms,
            double cutoff_short,
            double cutoff_long,
            double box_x,
            double box_y,
            double box_z,
            unsigned int num_atoms_per_mol,
            unsigned int estimated_neighbors_short,
            unsigned int estimated_neighbors_long,
            double crf_2cut3i, double crf_cut, double crf_cut3i,
            cudakernel::lj_crf_parameter * lj_crf_params,
            unsigned int num_of_gpus,
            unsigned int gpu_id,
            int * error) {
  if (device_number != -1) {
    cudaSetDevice(device_number);
  } else {
    // determine the device number
    cudaGetDevice(&device_number);
  }
  // let's first query the device properties and print them out
  cudaDeviceProp devProp;
  cudaGetDeviceProperties(&devProp, device_number);
  unsigned flags = 0;
  cudaGetDeviceFlags(&flags);
  if (flags == 0 && cudaSetDeviceFlags(cudaDeviceScheduleYield) == cudaErrorSetOnActiveProcess)
  std::cerr << "Cannot set flags\n";

  std::cout << "CUDA" << std::endl;
  std::cout << "\tDeviceproperties:" << std::endl;
  std::cout << "\t\tNumber: " << device_number << std::endl;
  std::cout << "\t\tName: " << devProp.name << std::endl;
  std::cout << "\t\tTotal memory: " << devProp.totalGlobalMem << std::endl;
  std::cout << "\t\tShared memory per block: " << devProp.sharedMemPerBlock << std::endl;
  std::cout << "\t\tTotal constant memory: " << devProp.totalConstMem << std::endl;

  // Create a gpu_status structure
  gpu_status * gpu_stat;
  gpu_stat = (gpu_status*) malloc (sizeof(gpu_status));
  gpu_stat->device = device_number;
  
  // assign all those simulation parameters to the struct
  gpu_stat->host_parameter.cutoff_long = (float) cutoff_long;
  gpu_stat->host_parameter.cutoff_long_2 = (float) cutoff_long*cutoff_long;
  gpu_stat->host_parameter.cutoff_short = (float) cutoff_short;
  gpu_stat->host_parameter.cutoff_short_2 = (float) cutoff_short*cutoff_short;
  gpu_stat->host_parameter.box_x = (float) box_x;
  gpu_stat->host_parameter.box_y = (float) box_y;
  gpu_stat->host_parameter.box_z = (float) box_z;
  gpu_stat->host_parameter.box_inv_x = 1.0 / box_x;
  gpu_stat->host_parameter.box_inv_y = 1.0 / box_y;
  gpu_stat->host_parameter.box_inv_z = 1.0 / box_z;
  gpu_stat->host_parameter.box_half_x = box_x / 2.0;
  gpu_stat->host_parameter.box_half_y = box_y / 2.0;
  gpu_stat->host_parameter.box_half_z = box_z / 2.0;
  gpu_stat->host_parameter.crf_2cut3i = crf_2cut3i;
  gpu_stat->host_parameter.crf_cut = crf_cut;
  gpu_stat->host_parameter.crf_cut3i = crf_cut3i;
  gpu_stat->host_parameter.num_atoms = num_atoms;
  gpu_stat->host_parameter.estimated_neighbors_short = estimated_neighbors_short;
  gpu_stat->host_parameter.estimated_neighbors_long = estimated_neighbors_long;
  gpu_stat->host_parameter.num_atoms_per_mol = num_atoms_per_mol;
  gpu_stat->host_parameter.num_solvent_mol = num_atoms / num_atoms_per_mol;

  // Assign the numbers of gpus and the gpu_id to the struct
  gpu_stat->host_parameter.num_of_gpus = num_of_gpus;
  gpu_stat->host_parameter.gpu_id = gpu_id;


  // Allocate memory for the simulation parameters
  unsigned int mem = 0;
  #define GPUMALLOC(x, size) cudaMalloc((void**)(x), (size)); mem += (size)
  GPUMALLOC(& gpu_stat->dev_parameter, sizeof (cudakernel::simulation_parameter));

  // Allocate memory for positions, pairlists, parameters, forces and energies
  GPUMALLOC(& gpu_stat->dev_pos, gpu_stat->host_parameter.num_atoms * sizeof(float3));
  GPUMALLOC(& gpu_stat->dev_new_pos, gpu_stat->host_parameter.num_atoms * sizeof(float3));

  // assign two dimensional arrays for the pairlists.
  unsigned int numThreads = (gpu_stat->host_parameter.num_solvent_mol / num_of_gpus + 1) * gpu_stat->host_parameter.num_atoms_per_mol;
  allocate_pairlist(gpu_stat->dev_pl_short, gpu_stat->host_parameter.num_solvent_mol / num_of_gpus + 1, gpu_stat->host_parameter.estimated_neighbors_short);
  mem += gpu_stat->host_parameter.num_solvent_mol * gpu_stat->host_parameter.estimated_neighbors_short / num_of_gpus * sizeof (unsigned int);
  allocate_pairlist(gpu_stat->dev_pl_long, gpu_stat->host_parameter.num_solvent_mol / num_of_gpus + 1, gpu_stat->host_parameter.estimated_neighbors_long);
  mem += gpu_stat->host_parameter.num_solvent_mol * gpu_stat->host_parameter.estimated_neighbors_long / num_of_gpus * sizeof (unsigned int);

  // allocate space for forces, virial and energies
  // and set it to zero
  GPUMALLOC(& gpu_stat->dev_forces, numThreads * sizeof (float3));
  cudaMemset(gpu_stat->dev_forces, 0, numThreads * sizeof (float3));
  GPUMALLOC(& gpu_stat->dev_virial, numThreads * sizeof (float9));
  cudaMemset(gpu_stat->dev_virial, 0, numThreads * sizeof (float9));
  GPUMALLOC(& gpu_stat->dev_energy, numThreads * sizeof (float2));
  cudaMemset(gpu_stat->dev_energy, 0, numThreads * sizeof (float2));
  // allocate space for parameters
  GPUMALLOC(& gpu_stat->dev_lj_crf_parameter, gpu_stat->host_parameter.num_atoms_per_mol * gpu_stat->host_parameter.num_atoms_per_mol * sizeof (cudakernel::lj_crf_parameter));

  std::cout << "\t\tMemory used: " << mem << " bytes or "
          << mem / (1024.0*1024.0) << " MB" << std::endl;

  // allocate data structure suited for copying on the host
  gpu_stat->host_energy = (float2 *) malloc(numThreads * sizeof (float2));
  gpu_stat->host_pos = (float3 *) malloc(gpu_stat->host_parameter.num_atoms * sizeof (float3));
  gpu_stat->host_forces = (float3 *) malloc(numThreads * sizeof (float3));
  gpu_stat->host_virial = (float9 *) malloc(numThreads * sizeof (float9));

  // copy over the parameters
  cudaMemcpy(gpu_stat->dev_lj_crf_parameter, lj_crf_params, gpu_stat->host_parameter.num_atoms_per_mol * gpu_stat->host_parameter.num_atoms_per_mol * sizeof (cudakernel::lj_crf_parameter), cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_stat->dev_parameter, &gpu_stat->host_parameter, sizeof (cudakernel::simulation_parameter), cudaMemcpyHostToDevice);

  std::cout << "END" << std::endl;
  *error = cudakernel::checkError("after init");

  return gpu_stat;
}

extern "C" int cudaCopyBox(gpu_status * gpu_stat, double box_x, double box_y, double box_z) {
  gpu_stat->host_parameter.box_x = (float) box_x;
  gpu_stat->host_parameter.box_y = (float) box_y;
  gpu_stat->host_parameter.box_z = (float) box_z;
  gpu_stat->host_parameter.box_inv_x = 1.0 / box_x;
  gpu_stat->host_parameter.box_inv_y = 1.0 / box_y;
  gpu_stat->host_parameter.box_inv_z = 1.0 / box_z;
  gpu_stat->host_parameter.box_half_x = box_x / 2.0;
  gpu_stat->host_parameter.box_half_y = box_y / 2.0;
  gpu_stat->host_parameter.box_half_z = box_z / 2.0;
  cudaMemcpy(gpu_stat->dev_parameter, &gpu_stat->host_parameter, sizeof (cudakernel::simulation_parameter), cudaMemcpyHostToDevice);
  return cudakernel::checkError("after copying the box");
}

extern "C" int cudaCopyPositions(double * pos, gpu_status * gpu_stat) {
  double3 * positions = (double3*) pos;
  for (unsigned int i = 0; i < gpu_stat->host_parameter.num_atoms; i++) {
    gpu_stat->host_pos[i].x = positions[i].x;
    gpu_stat->host_pos[i].y = positions[i].y;
    gpu_stat->host_pos[i].z = positions[i].z;
  }

  cudaMemcpy(gpu_stat->dev_pos, gpu_stat->host_pos, gpu_stat->host_parameter.num_atoms * sizeof (float3), cudaMemcpyHostToDevice);
  return cudakernel::checkError("after copying the positions");
}

extern "C" int CleanUp(gpu_status * gpu_stat) {
  free(gpu_stat->host_energy);
  free(gpu_stat->host_forces);
  free(gpu_stat->host_virial);
  free(gpu_stat->host_pos);

  cudaFree(gpu_stat->dev_pos);
  cudaFree(gpu_stat->dev_forces);
  cudaFree(gpu_stat->dev_virial);
  free_pairlist(gpu_stat->dev_pl_short);
  free_pairlist(gpu_stat->dev_pl_long);
  cudaFree(gpu_stat->dev_lj_crf_parameter);
  cudaFree(gpu_stat->dev_energy);
  cudaFree(gpu_stat->dev_parameter);

  return cudakernel::checkError("after clean-up");
}

