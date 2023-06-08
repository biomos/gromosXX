/**
 * @file interaction.cu
 * interaction compuation
 */

#include "gpu_status.h"
#include "interaction.h"
#include "lib/utils.h"
#include "lib/math.h"
#include "lib/reduction.h"
#include "../util/debug.h"

#undef MODULE
#undef SUBMODULE
#define MODULE cuda
#define SUBMODULE interaction

#define NUM_THREADS_PER_BLOCK_FORCES 32

extern __device__ __constant__ cudakernel::simulation_parameter device_param;
extern __device__ __constant__ cudakernel::simulation_parameter::box_struct device_box;

extern "C" int cudaCalcForces(double * forces, double * virial, double * lj_energy,
        double * crf_energy, bool longrange, gpu_status * gpu_stat) {
  int error = 0;
  cudakernel::pairlist *dev_pl;
  int num_of_gpus = gpu_stat->host_parameter.num_of_gpus;
  int gpu_id = gpu_stat->host_parameter.gpu_id;
  DEBUG(10,"Num solvent atoms: " << gpu_stat->host_parameter.num_atoms.solvent);

  cudakernel::simulation_parameter tmp_param;
  cudaMemcpyFromSymbol(&tmp_param, device_param, sizeof(cudakernel::simulation_parameter));
  DEBUG(4,"num_atoms.solvent: " << gpu_stat->host_parameter.num_atoms.solvent << " " << tmp_param.num_atoms.solvent);
  DEBUG(4,"num_atoms_per_mol: " << gpu_stat->host_parameter.num_atoms_per_mol << " " << tmp_param.num_atoms_per_mol);
  DEBUG(4,"num_of_gpus: " << gpu_stat->host_parameter.num_of_gpus << " " << tmp_param.num_of_gpus);
  unsigned int numThreads = (gpu_stat->host_parameter.num_solvent_mol / num_of_gpus + 1) * gpu_stat->host_parameter.num_atoms_per_mol;
  unsigned int numBlocks = numThreads / NUM_THREADS_PER_BLOCK_FORCES + 1;

  dim3 dimGrid(numBlocks, 1);
  dim3 dimBlock(NUM_THREADS_PER_BLOCK_FORCES, 1);


  DEBUG(10,"Interactions: GPU: " << gpu_id << " of " << num_of_gpus << " Threads: " << numThreads << " Blocks: " << numBlocks)
          
  // decide which pairlist to take.
  dev_pl = longrange ? &gpu_stat->dev_pl_long : &gpu_stat->dev_pl_short;

  error += cudakernel::check_error("before calculating Forces");
  // make sure to allocate __shared__ memory for ever thread
  // we need: 1. for every thread (i.e. NUM_THREADS.. * num_atoms_per_mol
  //          2. for every atom in a solvent molecule
  //          3. the size of the LJ CRF parameter struct
  const size_t dynamic_shared_memory = NUM_THREADS_PER_BLOCK_FORCES * gpu_stat->host_parameter.num_atoms_per_mol * sizeof(cudakernel::lj_crf_parameter);

  // check what we have there
  cudaMemcpyFromSymbol(&tmp_param, device_param, sizeof(cudakernel::simulation_parameter));
  DEBUG(0, "num_atoms.total: \t" << tmp_param.num_atoms.total);
  DEBUG(0, "num_atoms.solute: \t" << tmp_param.num_atoms.solute);
  DEBUG(0, "num_atoms.solvent: \t" << tmp_param.num_atoms.solvent);
  DEBUG(0, "box: \t" << tmp_param.box.full.x << " " << tmp_param.box.full.y << " " << tmp_param.box.full.z);
  DEBUG(0, "box.inv: \t" << tmp_param.box.inv.x << " " << tmp_param.box.inv.y << " " << tmp_param.box.inv.z);
  DEBUG(0, "box.half: \t" << tmp_param.box.half.x << " " << tmp_param.box.half.y << " " << tmp_param.box.half.z);
  DEBUG(0, "cutoff_long: \t" << tmp_param.cutoff_long);
  DEBUG(0, "cutoff_long_2: \t" << tmp_param.cutoff_long_2);
  DEBUG(0, "cutoff_short: \t" << tmp_param.cutoff_short);
  DEBUG(0, "cutoff_short_2: \t" << tmp_param.cutoff_short_2);
  DEBUG(0, "crf_2cut3i: \t" << tmp_param.crf_2cut3i);
  DEBUG(0, "crf_cut: \t" << tmp_param.crf_cut);
  DEBUG(0, "crf_cut3i: \t" << tmp_param.crf_cut3i);
  DEBUG(0, "num_atoms_per_mol: \t" << tmp_param.num_atoms_per_mol);
  DEBUG(0, "num_solvent_mol: \t" << tmp_param.num_solvent_mol);
  DEBUG(0, "estimated_neighbors_long: \t" << tmp_param.estimated_neighbors_long);
  DEBUG(0, "estimated_neighbors_short: \t" << tmp_param.estimated_neighbors_short);
  DEBUG(0, "num_of_gpus: \t" << tmp_param.num_of_gpus);
  DEBUG(0, "gpu_id: \t" << tmp_param.gpu_id);

  cudakernel::kernel_CalcForces_Solvent <<<dimGrid, dimBlock, dynamic_shared_memory >>>
          (*dev_pl, gpu_stat->dev_pos, gpu_stat->dev_parameter,
          gpu_stat->dev_lj_crf_parameter, gpu_stat->dev_forces,
          gpu_stat->dev_virial, gpu_stat->dev_energy,
          num_of_gpus,
          gpu_id);

  //cudaDeviceSynchronize();

  DEBUG(10,"Interactions: Executed kernel and synchronized threads")

  error += cudakernel::check_error("after calculating Forces");
  //copy forces/energies back to CPU
  // TODO: These could be reduced on GPU - numThreads-fold less communication
  cudaMemcpy(gpu_stat->host_forces, gpu_stat->dev_forces, numThreads * sizeof (float3), cudaMemcpyDeviceToHost);
  error += cudakernel::check_error("after copying Forces");
  #ifndef NDEBUG
  cudaMemcpy(gpu_stat->host_virial, gpu_stat->dev_virial, numThreads * sizeof (float9), cudaMemcpyDeviceToHost);
  error += cudakernel::check_error("after copying virial");
  cudaMemcpy(gpu_stat->host_energy, gpu_stat->dev_energy, numThreads * sizeof (float2), cudaMemcpyDeviceToHost);
  error += cudakernel::check_error("after copying energy");
  #endif
  
  //sum up energies on gpu
  float2 energy_sum = cudakernel::calc_sum<float2,64,1>(numThreads, gpu_stat->dev_energy, gpu_stat->dev_energy_output);

  *lj_energy = energy_sum.x;
  *crf_energy = energy_sum.y;
  
  #ifndef NDEBUG
  double l_lj_energy = 0.0;
  double l_crf_energy = 0.0;

  //sum up energies and write them back
  for (unsigned int i = 0; i < numThreads; i++) {
    //(*lj_energy) += gpu_stat->host_energy[i].x;
    //(*crf_energy) += gpu_stat->host_energy[i].y;
    l_lj_energy += gpu_stat->host_energy[i].x;
    l_crf_energy += gpu_stat->host_energy[i].y;
  }

  DEBUG(7,"GPU lj_energy: " << *lj_energy);
  DEBUG(7,"GPU crf_energy: " << *crf_energy);
  DEBUG(7,"CPU lj_energy: " << l_lj_energy);
  DEBUG(7,"CPU crf_energy: " << l_crf_energy);
  #endif


  //sum up forces -> include effects of solute
  double3 * pforces = (double3 *) forces;
  int index;
  for (unsigned int i = 0; i < numThreads; i++) {
    // convert from float to double
    index = i % gpu_stat->host_parameter.num_atoms_per_mol + ( (i / gpu_stat->host_parameter.num_atoms_per_mol) * num_of_gpus + gpu_id)*gpu_stat->host_parameter.num_atoms_per_mol;
    if (index >= gpu_stat->host_parameter.num_atoms.solvent)
      continue;
    //std::cout << gpu_id << " : i = " << i << " , index = " << index << std::endl;
    DEBUG(15,"force " << i << " index: " << (int) gpu_stat->host_forces[i].x << " ai: " << (int) gpu_stat->host_forces[i].y << " mi: " << (int) gpu_stat->host_forces[i].z);
    //std::cout << "GPU: " << gpu_id << std::endl;
    pforces[index].x += gpu_stat->host_forces[i].x;
    pforces[index].y += gpu_stat->host_forces[i].y;
    pforces[index].z += gpu_stat->host_forces[i].z;
  }

  //sum up virials on gpu
  float9 virial_sum = cudakernel::calc_sum<float9,64,1>(numThreads, gpu_stat->dev_virial, gpu_stat->dev_virial_output);

  virial[0] += virial_sum.xx;
  virial[1] += virial_sum.xy;
  virial[2] += virial_sum.xz;
  virial[3] += virial_sum.yx;
  virial[4] += virial_sum.yy;
  virial[5] += virial_sum.yz;
  virial[6] += virial_sum.zx;
  virial[7] += virial_sum.zy;
  virial[8] += virial_sum.zz;

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
  //const unsigned int solvent_offset = dev_params->num_atoms_per_mol;
  const unsigned int solvent_offset = device_param.num_atoms_per_mol;
  const int atom_index = index % solvent_offset + ((index / solvent_offset) * num_of_gpus + gpu_id) * solvent_offset;

  if (atom_index >= device_param.num_atoms.solvent)
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
  volatile float e_lj = 0.0f, e_crf = 0.0f;

  // the shared memory is allocated in the kernel call.
  extern __shared__ cudakernel::lj_crf_parameter lj_crf[];
  // now cache the parameter matrix line in the shared memory
  /*for (unsigned int j = 0; j < solvent_offset; ++j)
    lj_crf[threadIdx.x + NUM_THREADS_PER_BLOCK_FORCES * j] = dev_lj_crf_params[atom_type * solvent_offset + j];
*/
  // cache reaction field parameters
  const float crf_2cut3i = device_param.crf_2cut3i;
  const float crf_cut = device_param.crf_cut;
  const float crf_cut3i = device_param.crf_cut3i;

  // cache box
  /*const float3 box_param_x = make_float3(device_param.box.half.x, device_param.box.full.x, device_param.box.inv.x);
  const float3 box_param_y = make_float3(device_param.box.half.y, device_param.box.full.y, device_param.box.inv.y);
  const float3 box_param_z = make_float3(device_param.box.half.z, device_param.box.full.z, device_param.box.inv.z); */

  unsigned active = __activemask();
  // loop over neighbor list which contains the index of the first atom of the neighboring molecule
  for (unsigned int i = 0; __any_sync(active, i < num_neighbors); i++) {
    if (i < num_neighbors) {
      const unsigned int neighbor_index = pl.list[molecule_index + pl.pitch * i];
      // do we need to calculate the energy?
      bool calculate_energy = (atom_index  < neighbor_index);

      // loop over atoms in neighboring solvent molecule
      for (unsigned int j = 0; j < solvent_offset; ++j) {
        const unsigned int atom_j = neighbor_index + j;
        // calculate distance
        //const float3 r = nearestImage(my_pos, dev_pos[atom_j], box_param_x, box_param_y, box_param_z);
        const float3 r = nearestImage(my_pos, dev_pos[atom_j]);
        // get the parameters
        //const cudakernel::lj_crf_parameter lj_crf_param = lj_crf[threadIdx.x + NUM_THREADS_PER_BLOCK_FORCES * j];
        const cudakernel::lj_crf_parameter lj_crf_param = device_param.solvent_lj_crf[atom_type * solvent_offset + j];
        // calculate the force
        //const float dist2 = abs2(r);
        const float disti = rnorm3df(r.x,r.y,r.z);
        const float dist2i = disti * disti;
        const float dist6i = dist2i * dist2i * dist2i;
        //const float disti = sqrtf(dist2i);
        //const float disti = rsqrtf(dist2); // lower register pressure and bit more precise
        //const float disti = rnorm3df(r.x,r.y,r.z); // lower register pressure and bit more precise
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
          const float dist2 = abs2(r);
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

