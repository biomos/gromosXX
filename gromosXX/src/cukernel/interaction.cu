/**
 * @file interaction.cu
 * interaction compuation
 */

#include "gpu_status.h"
#include "interaction.h"
#include "lib/utils.h"
#include "lib/math.h"
#include "lib/reduction.h"
#include "lib/device_functions.h"
#include "../util/debug.h"

#undef MODULE
#undef SUBMODULE
#define MODULE cukernel
#define SUBMODULE interaction


extern __device__ __constant__ cukernel::simulation_parameter device_param;

extern "C" int cudaCalcForces(double * forces, double * virial, double * lj_energy,
        double * crf_energy, bool longrange, gpu_status * gpu_stat) {
  int error = 0;
  cukernel::pairlist *dev_pl;
  int num_of_gpus = gpu_stat->host_parameter.num_of_gpus;
  int gpu_id = gpu_stat->host_parameter.gpu_id;

  const unsigned int num_threads = (gpu_stat->host_parameter.num_solvent_mol / num_of_gpus + 1) * gpu_stat->host_parameter.num_atoms_per_mol;
  const unsigned int num_blocks = num_threads / NUM_THREADS_PER_BLOCK_FORCES + 1;
  const unsigned int num_warps = num_blocks * NUM_THREADS_PER_BLOCK_FORCES / 32 + 1;

  dim3 dimGrid(num_blocks, 1);
  dim3 dimBlock(NUM_THREADS_PER_BLOCK_FORCES, 1);


  DEBUG(10,"Interactions: GPU: " << gpu_id << " of " << num_of_gpus << " Threads: " << num_threads << " Blocks: " << num_blocks)
          
  // decide which pairlist to take.
  dev_pl = longrange ? &gpu_stat->dev_pl_long : &gpu_stat->dev_pl_short;

  error += cukernel::check_error("before calculating Forces");
  // make sure to allocate __shared__ memory for ever thread
  // we need: 1. for every thread (i.e. NUM_THREADS.. * num_atoms_per_mol
  //          2. for every atom in a solvent molecule
  //          3. the size of the LJ CRF parameter struct
  //const size_t dynamic_shared_memory = NUM_THREADS_PER_BLOCK_FORCES * gpu_stat->host_parameter.num_atoms_per_mol * sizeof(cukernel::lj_crf_parameter);

  #ifndef NDEBUG
  // Debug device constant memory
  cukernel::simulation_parameter tmp_param;
  cudaMemcpyFromSymbol(&tmp_param, device_param, sizeof(cukernel::simulation_parameter));
  #endif
  DEBUG(7, "num_atoms.total: \t" << tmp_param.num_atoms.total);
  DEBUG(7, "num_atoms.solute: \t" << tmp_param.num_atoms.solute);
  DEBUG(7, "num_atoms.solvent: \t" << tmp_param.num_atoms.solvent);
  DEBUG(7, "box: \t" << tmp_param.box.full.x << " " << tmp_param.box.full.y << " " << tmp_param.box.full.z);
  DEBUG(7, "box.inv: \t" << tmp_param.box.inv.x << " " << tmp_param.box.inv.y << " " << tmp_param.box.inv.z);
  DEBUG(7, "box.half: \t" << tmp_param.box.half.x << " " << tmp_param.box.half.y << " " << tmp_param.box.half.z);
  DEBUG(7, "cutoff_long: \t" << tmp_param.cutoff_long);
  DEBUG(7, "cutoff_long_2: \t" << tmp_param.cutoff_long_2);
  DEBUG(7, "cutoff_short: \t" << tmp_param.cutoff_short);
  DEBUG(7, "cutoff_short_2: \t" << tmp_param.cutoff_short_2);
  DEBUG(7, "crf_2cut3i: \t" << tmp_param.crf_2cut3i);
  DEBUG(7, "crf_cut: \t" << tmp_param.crf_cut);
  DEBUG(7, "crf_cut3i: \t" << tmp_param.crf_cut3i);
  DEBUG(7, "num_atoms_per_mol: \t" << tmp_param.num_atoms_per_mol);
  DEBUG(7, "num_solvent_mol: \t" << tmp_param.num_solvent_mol);
  DEBUG(7, "estimated_neighbors_long: \t" << tmp_param.estimated_neighbors_long);
  DEBUG(7, "estimated_neighbors_short: \t" << tmp_param.estimated_neighbors_short);
  DEBUG(7, "num_of_gpus: \t" << tmp_param.num_of_gpus);
  DEBUG(7, "gpu_id: \t" << tmp_param.gpu_id);


  cukernel::kernel_CalcForces_Solvent <<<dimGrid, dimBlock >>>
          (*dev_pl, gpu_stat->dev_pos, gpu_stat->dev_parameter,
          gpu_stat->dev_lj_crf_parameter, gpu_stat->dev_forces,
          gpu_stat->dev_virial, gpu_stat->dev_energy,
          num_of_gpus,
          gpu_id);

  DEBUG(10,"Interactions: Executed kernel and synchronized threads")

  error += cukernel::check_error("after calculating Forces");
  //copy forces/energies back to CPU
  cudaMemcpy(gpu_stat->host_forces, gpu_stat->dev_forces, num_threads * sizeof (float3), cudaMemcpyDeviceToHost);
  error += cukernel::check_error("after copying Forces");
  // Virials and energies are reduced on GPU, here copy only for debugging
  #ifndef NDEBUG
  cudaMemcpy(gpu_stat->host_virial, gpu_stat->dev_virial, num_warps * sizeof (float9), cudaMemcpyDeviceToHost);
  error += cukernel::check_error("after copying virial");
  cudaMemcpy(gpu_stat->host_energy, gpu_stat->dev_energy, num_warps * sizeof (float2), cudaMemcpyDeviceToHost);
  error += cukernel::check_error("after copying energy");
  #endif

  //sum up energies on gpu
  float2 energy_sum = cukernel::calc_sum<float2, NUM_THREADS_PER_BLOCK_REDUCTION_FLOAT2, NUM_ELEMENTS_PER_THREAD_FLOAT2>(num_warps, gpu_stat->dev_energy, gpu_stat->dev_energy_output);

  *lj_energy = energy_sum.x;
  *crf_energy = energy_sum.y;
  
  #ifndef NDEBUG
  double l_lj_energy = 0.0;
  double l_crf_energy = 0.0;

  //sum up energies and write them back
  for (unsigned int i = 0; i < num_warps; ++i) {
    l_lj_energy += gpu_stat->host_energy[i].x;
    l_crf_energy += gpu_stat->host_energy[i].y;
  }

  DEBUG(7,"lj_energy: CPU: " << l_lj_energy << ",\tGPU: " << *lj_energy);
  DEBUG(7,"lj_energy: CPU: " << l_crf_energy << ",\tGPU: " << *crf_energy);
  #endif


  //sum up forces -> include effects of solute
  double3 * pforces = (double3 *) forces;
  int index;
  for (unsigned int i = 0; i < num_threads; ++i) {
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
  float9 virial_sum = cukernel::calc_sum<float9, NUM_THREADS_PER_BLOCK_REDUCTION_FLOAT9, NUM_ELEMENTS_PER_THREAD_FLOAT9>(num_warps, gpu_stat->dev_virial, gpu_stat->dev_virial_output);
  
  #ifndef NDEBUG
  double9 l_virial_sum {0.0};
  for (unsigned int i = 0; i < num_warps; i++) {
    l_virial_sum.xx += gpu_stat->host_virial[i].xx;
    l_virial_sum.xy += gpu_stat->host_virial[i].xy;
    l_virial_sum.xz += gpu_stat->host_virial[i].xz;
    l_virial_sum.yx += gpu_stat->host_virial[i].yx;
    l_virial_sum.yy += gpu_stat->host_virial[i].yy;
    l_virial_sum.yz += gpu_stat->host_virial[i].yz;
    l_virial_sum.zx += gpu_stat->host_virial[i].zx;
    l_virial_sum.zy += gpu_stat->host_virial[i].zy;
    l_virial_sum.zz += gpu_stat->host_virial[i].zz;
  }
  DEBUG(7,"virial_sum.xx: CPU: " << l_virial_sum.xx << ",\tGPU: " << virial_sum.xx);
  DEBUG(7,"virial_sum.xy: CPU: " << l_virial_sum.xy << ",\tGPU: " << virial_sum.xy);
  DEBUG(7,"virial_sum.xz: CPU: " << l_virial_sum.xz << ",\tGPU: " << virial_sum.xz);
  DEBUG(7,"virial_sum.yx: CPU: " << l_virial_sum.yx << ",\tGPU: " << virial_sum.yx);
  DEBUG(7,"virial_sum.yy: CPU: " << l_virial_sum.yy << ",\tGPU: " << virial_sum.yy);
  DEBUG(7,"virial_sum.yz: CPU: " << l_virial_sum.yz << ",\tGPU: " << virial_sum.yz);
  DEBUG(7,"virial_sum.zx: CPU: " << l_virial_sum.zx << ",\tGPU: " << virial_sum.zx);
  DEBUG(7,"virial_sum.zy: CPU: " << l_virial_sum.zy << ",\tGPU: " << virial_sum.zy);
  DEBUG(7,"virial_sum.zz: CPU: " << l_virial_sum.zz << ",\tGPU: " << virial_sum.zz);
  #endif

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

__global__ void cukernel::kernel_CalcForces_Solvent
(
        pairlist pl,
        float3 * dev_pos,
        cukernel::simulation_parameter * dev_params,
        cukernel::lj_crf_parameter * dev_lj_crf_params,
        float3 * dev_for,
        float9 * dev_virial,
        float2 * dev_energy,
        unsigned int num_of_gpus,
        unsigned int gpu_id
) {
  // calculate atom index and check for boundaries
  // x: molecules, y: atoms of the molecule
  const unsigned index = blockIdx.x * NUM_THREADS_PER_BLOCK_FORCES + threadIdx.x;
  //const unsigned int solvent_offset = dev_params->num_atoms_per_mol;
  const unsigned solvent_offset = device_param.num_atoms_per_mol;
  const unsigned atom_index = index % solvent_offset + ((index / solvent_offset) * num_of_gpus + gpu_id) * solvent_offset;
  //printf("index: %u, blockIdx.x: %u, atom_index: %u, solvent_offset: %u\n", index, blockIdx.x, atom_index, solvent_offset);
  //printf("num_atoms: total %u, solute %u, solvent %u, per_mol %u\n", device_param.num_atoms.total, device_param.num_atoms.solute, device_param.num_atoms.solvent, device_param.num_atoms_per_mol);
  if (atom_index >= device_param.num_atoms.solvent)
    return;
  //printf("index: %u, blockIdx.x: %u, atom_index: %u, solvent_offset: %u\n", index, blockIdx.x, atom_index, solvent_offset);

  // fetch the position and some parameters
  const float3 my_pos = dev_pos[atom_index];

  const unsigned atom_type = atom_index % solvent_offset;
  const unsigned molecule_index = atom_index / (solvent_offset * num_of_gpus);
  //const unsigned int molecule_index = index;
  //dev_for[index] = make_float3(index, atom_index, molecule_index); return;
  const unsigned num_neighbors = pl.num_neighbors[molecule_index];

  // initialize storage for quantities computed
  float3 force {0.0f};
  float9 virial {0.0f};
  volatile float e_lj = 0.0f, e_crf = 0.0f;
  // the shared memory is allocated in the kernel call.
  //extern __shared__ cukernel::lj_crf_parameter lj_crf[];
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

  /*if (active != 0xFFFFFFFF)
    printf("index: %u, warp: %u, lane: %u, active: 0x%08X\n", index, index/32, index%32, active);
  if (solvent_offset != 3)
    printf("index: %u, warp: %u, lane: %u, active: 0x%08X, solvent_offset: %u\n", index, index/32, index%32, active, solvent_offset);*/

  // loop over neighbor list which contains the index of the first atom of the neighboring molecule
  for (unsigned i = 0; __any_sync(__activemask(), i < num_neighbors); ++i) {
    
    /*if (__activemask() != 0xFFFFFFFF) {
      printf("index: %u, warp: %u, lane: %u, active: 0x%08X\n", index, index/32, index%32, active);
    }*/
    if (i < num_neighbors) {
      const unsigned neighbor_index = pl.list[molecule_index + pl.pitch * i];
      // do we need to calculate the energy?
      bool calculate_energy = (atom_index < neighbor_index);

      // loop over atoms in neighboring solvent molecule
      for (unsigned j = 0; j < solvent_offset; ++j) {
        const unsigned atom_j = neighbor_index + j;
        // calculate distance
        //const float3 r = nearestImage(my_pos, dev_pos[atom_j], box_param_x, box_param_y, box_param_z);
        const float3 r = nearestImage(my_pos, dev_pos[atom_j]);
        // get the parameters
        //const cukernel::lj_crf_parameter lj_crf_param = lj_crf[threadIdx.x + NUM_THREADS_PER_BLOCK_FORCES * j];
        const cukernel::lj_crf_parameter lj_crf_param = device_param.solvent_lj_crf[atom_type * solvent_offset + j];
        //const cukernel::lj_crf_parameter lj_crf_param {0.0f};
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


  dev_for[index] = force;
  //dev_energy[index] = make_float2(e_lj, e_crf);
  //dev_virial[index] = virial;
  // reduce virial and energy within warp
  const unsigned mask = __ballot_sync(FULL_MASK, index < device_param.num_atoms.solvent);
  for (unsigned offset = 16; offset > 0; offset /= 2) {
    // only if my partner participate
    //printf("atom_index: %d, threadIdx.x: %d, active: %0X\n", atom_index, threadIdx.x, active);
    float _e_lj, _e_crf;
    float9 _virial;
    _e_lj      = __shfl_down_sync(mask, e_lj, offset);
    _e_crf     = __shfl_down_sync(mask, e_crf, offset);
    _virial.xx = __shfl_down_sync(mask, virial.xx, offset);
    _virial.xy = __shfl_down_sync(mask, virial.xy, offset);
    _virial.xz = __shfl_down_sync(mask, virial.xz, offset);
    _virial.yx = __shfl_down_sync(mask, virial.yx, offset);
    _virial.yy = __shfl_down_sync(mask, virial.yy, offset);
    _virial.yz = __shfl_down_sync(mask, virial.yz, offset);
    _virial.zx = __shfl_down_sync(mask, virial.zx, offset);
    _virial.zy = __shfl_down_sync(mask, virial.zy, offset);
    _virial.zz = __shfl_down_sync(mask, virial.zz, offset);
    // shuffling from above num_atoms.solvent is garbage
    if (index + offset < device_param.num_atoms.solvent) {
      e_lj += _e_lj;
      e_crf += _e_crf;
      virial += _virial;
    }
  }
  // elect a leader that who will write to global memory
  int leader = __ffs(mask) - 1;
  int lane_id = threadIdx.x%32;
  if (lane_id == leader) {
    dev_virial[index/32] = virial;
    dev_energy[index/32] = make_float2(e_lj, e_crf);
  }
}
