/**
 * @file interaction.cu
 * interaction compuation
 */

#include "gpu_status.h"
#include "interaction.h"
#include "lib/utils.h"
#include "lib/math.h"
#include "../util/debug.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE cuda


#define NUM_THREADS_PER_BLOCK_FORCES 96



extern "C" int cudaCalcForces(double * forces, double * virial, double * lj_energy,
        double * crf_energy, bool longrange, gpu_status * gpu_stat) {
  int error = 0;
  cudakernel::pairlist *dev_pl;
  int num_of_gpus = gpu_stat->host_parameter.num_of_gpus;
  int gpu_id = gpu_stat->host_parameter.gpu_id;
  DEBUG(10,"Num solvent atoms: " << gpu_stat->host_parameter.num_atoms);
  unsigned int numThreads = (gpu_stat->host_parameter.num_solvent_mol / num_of_gpus + 1) * gpu_stat->host_parameter.num_atoms_per_mol;
  unsigned int numBlocks = numThreads / NUM_THREADS_PER_BLOCK_FORCES + 1;

  dim3 dimGrid(numBlocks, 1);
  dim3 dimBlock(NUM_THREADS_PER_BLOCK_FORCES, 1);


  DEBUG(10,"Interactions: GPU: " << gpu_id << " of " << num_of_gpus << " Threads: " << numThreads << " Blocks: " << numBlocks)
          
  // decide which pairlist to take.
  dev_pl = longrange ? &gpu_stat->dev_pl_long : &gpu_stat->dev_pl_short;

  error += cudakernel::checkError("before calculating Forces");
  // make sure to allocate __shared__ memory for ever thread
  // we need: 1. for every thread (i.e. NUM_THREADS.. * num_atoms_per_mol
  //          2. for every atom in a solvent molecule
  //          3. the size of the LJ CRF parameter struct
  const size_t dynamic_shared_memory = NUM_THREADS_PER_BLOCK_FORCES * gpu_stat->host_parameter.num_atoms_per_mol * sizeof(cudakernel::lj_crf_parameter);

  cudakernel::kernel_CalcForces_Solvent <<<dimGrid, dimBlock, dynamic_shared_memory >>>
          (*dev_pl, gpu_stat->dev_pos, gpu_stat->dev_parameter,
          gpu_stat->dev_lj_crf_parameter, gpu_stat->dev_forces,
          gpu_stat->dev_virial, gpu_stat->dev_energy,
          num_of_gpus,
          gpu_id);

  cudaDeviceSynchronize();

  DEBUG(10,"Interactions: Executed kernel and synchronized threads")

  error += cudakernel::checkError("after calculating Forces");
  //copy forces/energies back to CPU
  // TODO: These could be reduced on GPU - numThreads-fold less communication
  cudaMemcpy(gpu_stat->host_forces, gpu_stat->dev_forces, numThreads * sizeof (float3), cudaMemcpyDeviceToHost);
  error += cudakernel::checkError("after copying Forces");
  cudaMemcpy(gpu_stat->host_virial, gpu_stat->dev_virial, numThreads * sizeof (float9), cudaMemcpyDeviceToHost);
  error += cudakernel::checkError("after copying virial");
  cudaMemcpy(gpu_stat->host_energy, gpu_stat->dev_energy, numThreads * sizeof (float2), cudaMemcpyDeviceToHost);
  error += cudakernel::checkError("after copying energy");

  
  //sum up energies and write them back
  for (unsigned int i = 0; i < numThreads; i++) {
    (*lj_energy) += gpu_stat->host_energy[i].x;
    (*crf_energy) += gpu_stat->host_energy[i].y;
  }
 
  //sum up forces -> include effects of solute
  double3 * pforces = (double3 *) forces;
  int index;
  for (unsigned int i = 0; i < numThreads; i++) {
    // convert from float to double
    index = i % gpu_stat->host_parameter.num_atoms_per_mol + ( (i / gpu_stat->host_parameter.num_atoms_per_mol) * num_of_gpus + gpu_id)*gpu_stat->host_parameter.num_atoms_per_mol;
    if (index >= gpu_stat->host_parameter.num_atoms)
      continue;
    //std::cout << gpu_id << " : i = " << i << " , index = " << index << std::endl;
    DEBUG(15,"force " << i << " index: " << (int) gpu_stat->host_forces[i].x << " ai: " << (int) gpu_stat->host_forces[i].y << " mi: " << (int) gpu_stat->host_forces[i].z);
    //std::cout << "GPU: " << gpu_id << std::endl;
    pforces[index].x += gpu_stat->host_forces[i].x;
    pforces[index].y += gpu_stat->host_forces[i].y;
    pforces[index].z += gpu_stat->host_forces[i].z;

    // and the same for the virial
    virial[0] += gpu_stat->host_virial[i].xx;
    virial[1] += gpu_stat->host_virial[i].xy;
    virial[2] += gpu_stat->host_virial[i].xz;
    virial[3] += gpu_stat->host_virial[i].yx;
    virial[4] += gpu_stat->host_virial[i].yy;
    virial[5] += gpu_stat->host_virial[i].yz;
    virial[6] += gpu_stat->host_virial[i].zx;
    virial[7] += gpu_stat->host_virial[i].zy;
    virial[8] += gpu_stat->host_virial[i].zz;
    // DEBUGGING
    //if (i%100 < 2)
    //  DEBUG(15,"GPU ID: " << gpu_stat->host_parameter.gpu_id << " of " << gpu_stat->host_parameter.num_of_gpus
    //          << " Interactions: Summing up forces and virals. Atom: " << i)
  }

  return error;

}

__global__ void cudakernel::kernel_CalcForces_Solvent
(
        pairlist pl,
        float3 * dev_pos,
        cudakernel::simulation_parameter * dev_params,
        cudakernel::lj_crf_parameter * dev_lj_crf_params,
        float3 * dev_for,
        float9 * dev_virial,
        float2 * dev_energy,
        unsigned int num_of_gpus,
        unsigned int gpu_id
) {
  // calculate atom index and check for boundaries
  // x: molecules, y: atoms of the molecule
  const int index = blockIdx.x * NUM_THREADS_PER_BLOCK_FORCES + threadIdx.x;
  const unsigned int solvent_offset = dev_params->num_atoms_per_mol;
  const int atom_index = index % solvent_offset + ((index / solvent_offset) * num_of_gpus + gpu_id) * solvent_offset;

  if (atom_index >= dev_params->num_atoms)
    return;

  // fetch the position and some parameters
  const float3 my_pos = dev_pos[atom_index];

  const unsigned int atom_type = atom_index % solvent_offset;
  const unsigned int molecule_index = atom_index / (solvent_offset * num_of_gpus);
  //const unsigned int molecule_index = index;
  //dev_for[index] = make_float3(index, atom_index, molecule_index); return;
  const unsigned int num_neighbors = pl.num_neighbors[molecule_index];

  // initialize storage for quantities computed
  float3 force = make_float3(0.0f, 0.0f, 0.0f);
  float9 virial;
  virial.xx = 0.0f;
  virial.xy = 0.0f;
  virial.xz = 0.0f;
  virial.yx = 0.0f;
  virial.yy = 0.0f;
  virial.yz = 0.0f;
  virial.zx = 0.0f;
  virial.zy = 0.0f;
  virial.zz = 0.0f;
  float e_lj = 0.0f, e_crf = 0.0f;

  // the shared memory is allocated in the kernel call.
  extern __shared__ cudakernel::lj_crf_parameter lj_crf[];
  // now cache the parameter matrix line in the shared memory
  for (unsigned int j = 0; j < solvent_offset; ++j)
    lj_crf[threadIdx.x + NUM_THREADS_PER_BLOCK_FORCES * j] = dev_lj_crf_params[atom_type * solvent_offset + j];

  // cache reaction field parameters
  const float crf_2cut3i = dev_params->crf_2cut3i;
  const float crf_cut = dev_params->crf_cut;
  const float crf_cut3i = dev_params->crf_cut3i;

  // cache box
  const float3 box_param_x = make_float3(dev_params->box_half_x, dev_params->box_x, dev_params->box_inv_x);
  const float3 box_param_y = make_float3(dev_params->box_half_y, dev_params->box_y, dev_params->box_inv_y);
  const float3 box_param_z = make_float3(dev_params->box_half_z, dev_params->box_z, dev_params->box_inv_z); 

  unsigned active = __activemask();
  // loop over neighbor list which contains the index of the first atom of the neighboring molecule
  for (unsigned int i = 0; __any_sync(active, i < num_neighbors); i++) {
    if (i < num_neighbors) {
      const unsigned int neighbor_index = pl.list[molecule_index + pl.pitch * i];
      // do we need to calculate the energy?
      bool calculate_energy = (atom_index  < neighbor_index);

      // loop over atoms in neighboring solvent molecule
      for (unsigned int j = 0; j < solvent_offset; j++) {
        const unsigned int atom_j = neighbor_index + j;
        // calculate distance
        const float3 r = nearestImage(my_pos, dev_pos[atom_j], box_param_x, box_param_y, box_param_z);
        // get the parameters
        const cudakernel::lj_crf_parameter lj_crf_param = lj_crf[threadIdx.x + NUM_THREADS_PER_BLOCK_FORCES * j];

        // calculate the force
        const float dist2 = abs2(r);
        const float dist2i = 1.0f / dist2;
        const float dist6i = dist2i * dist2i * dist2i;
        const float disti = sqrtf(dist2i);
        const float q_eps = lj_crf_param.q;
        const float c12_dist6i = lj_crf_param.c12 * dist6i;

        const float pair_force = (c12_dist6i + c12_dist6i - lj_crf_param.c6) * 6.0f * dist6i * dist2i +
                q_eps * (disti * dist2i + crf_cut3i);

        const float3 f_vec = make_float3(pair_force * r.x, pair_force * r.y, pair_force * r.z);
        force.x += f_vec.x;
        force.y += f_vec.y;
        force.z += f_vec.z;

        // calculate the energy and virial?
        if (calculate_energy) {
          e_lj += (c12_dist6i - lj_crf_param.c6) * dist6i;
          e_crf += q_eps * (disti - crf_2cut3i * dist2 - crf_cut);

          virial.xx += f_vec.x * r.x;
          virial.xy += f_vec.y * r.x;
          virial.xz += f_vec.z * r.x;
          virial.yx += f_vec.x * r.y;
          virial.yy += f_vec.y * r.y;
          virial.yz += f_vec.z * r.y;
          virial.zx += f_vec.x * r.z;
          virial.zy += f_vec.y * r.z;
          virial.zz += f_vec.z * r.z;
        }
      } //end molecule-loop
    }
  } //end pairlist-loop

  /*
  dev_energy[atom_index] = make_float2(e_lj, e_crf);
  dev_for[atom_index] = force;
  dev_virial[atom_index] = virial;
   */
  dev_energy[index] = make_float2(e_lj, e_crf);
  dev_for[index] = force;
  dev_virial[index] = virial;

}

#undef DEBUG

