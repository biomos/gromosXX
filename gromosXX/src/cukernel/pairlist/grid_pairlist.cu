/**
 * @file pairlist.cu
 * pairlist computation
 */

#include <iostream>
#include "gpu_status.h"

#include "lib/math.h"
#include "../util/debug.h"

#undef MODULE
#undef SUBMODULE
#define MODULE cuda
#define SUBMODULE pairlist

#define NUM_THREADS_PER_BLOCK 96

extern __device__ __constant__ cukernel::simulation_parameter device_param;


extern "C" void cuda_grid_pairlist(gpu_status * gpu_stat) {
  // prepare
   // get positions of chargegroups/atoms
  math::VArray pos;
  
  // assign chargegroups/atoms to subboxes
  assign_subbox<<<assignment_grid, assignment_block>>>(pos, subboxes, sb.sizes(), fpos.size());
  // array of size=num_chargegroups(or num_atoms)
  // chargegroups/atoms are sorted by their affiliation to the subbox
  // calculate pointers to subboxes

  // generate masks of contact shifts

  // generate pairlist
}

__global__ void assign_subbox(float3 *pos, unsigned *subbox, unsigned *sb_sizes, const unsigned N) {
  unsigned idx = threadIdx.x + blockIdx.x * blockDim.x;
  if (idx < N) {
    // Move to first quadrant box
    while (pos[idx].x < 0.0) {
      pos[idx].x += device_param.box.full.x;
    }
    while (pos[idx].x >= device_param.box.full.x) {
      pos[idx].x -= device_param.box.full.x;
    }

    while (pos[idx].y < 0.0) {
      pos[idx].y += device_param.box.full.y;
    }
    while (pos[idx].y >= device_param.box.full.y) {
      pos[idx].y -= device_param.box.full.y;
    }
    
    while (pos[idx].z < 0.0) {
      pos[idx].z += device_param.box.full.z;
    }
    while (pos[idx].z >= device_param.box.full.z) {
      pos[idx].z -= device_param.box.full.z;
    }

    uint3 sb;
    sb.x = pos[idx].x / device_param.box.full.x;
    sb.y = pos[idx].y / device_param.box.full.y;
    sb.z = pos[idx].z / device_param.box.full.z;
    unsigned subbox_idx = sb.x * device_param.num_subbox.y * device_param.num_subbox.z + sb.y * device_param.num_subbox.z + sb.z;
    subbox[idx] = subbox_idx;
    unsigned s = atomicAdd(&sb_sizes[subbox_idx], 1);
  }
};




void cukernel::free_pairlist(pairlist &pl) {
  cudaFree(pl.list); cudaFree(pl.num_neighbors); cudaFree(pl.overflow);
  DEBUG(4,"Pairlist: freed pairlist")
}
void cukernel::allocate_pairlist(pairlist &pl, unsigned int size, unsigned int max_neighbors) {
  size_t pitch;
  // allocate the number of neighbors
  cudaMalloc((void**) & pl.num_neighbors, size * sizeof (unsigned int));
  // Set the memory to 0
  cudaMemset(pl.num_neighbors, 0, size * sizeof (unsigned int));
  // allocate the neighbor list
  cudaMallocPitch((void**)((void*)&pl.list), &pitch, size * sizeof(unsigned int), max_neighbors);
  pl.max_size = max_neighbors;
  pl.pitch = (int)pitch / sizeof(int);
  // allocate the overflow flag and set it to false
  cudaMalloc((void**) &pl.overflow, sizeof(bool));
  // Set the memory to 0
  cudaMemset(pl.overflow, 0, sizeof(bool));
  DEBUG(10,"Pairlist: allocated memory")
}

extern "C" void cudaCalcPairlist(gpu_status * gpu_stat) {

  unsigned int num_blocks = (unsigned int) gpu_stat->host_parameter.num_solvent_mol / ( NUM_THREADS_PER_BLOCK * gpu_stat->host_parameter.num_of_gpus ) + 1;
  dim3 dimGrid(num_blocks, 1);
  dim3 dimBlock(NUM_THREADS_PER_BLOCK, 1);

  DEBUG(10,"Pairlist: GPU ID: " << gpu_stat->host_parameter.gpu_id << " of " << gpu_stat->host_parameter.num_of_gpus
            <<  ". Blocks: " << num_blocks)
  bool overflow;
  do {
    overflow = false;
    // calculate the pairlist
    kernel_CalcPairlist <<<dimGrid, dimBlock >>>(gpu_stat->dev_parameter,
            gpu_stat->dev_pos, 
            gpu_stat->dev_pl_short, gpu_stat->dev_pl_long,
            gpu_stat->host_parameter.num_of_gpus, gpu_stat->host_parameter.gpu_id);

    
    //cudaDeviceSynchronize();
    DEBUG(10,"Pairlist: Executed kernel and synchronized Threads")

    bool overflow_short, overflow_long;
    // get the overflow flags
    cudaMemcpy(&overflow_short, gpu_stat->dev_pl_short.overflow, sizeof (bool), cudaMemcpyDeviceToHost);
    cudaMemcpy(&overflow_long, gpu_stat->dev_pl_long.overflow, sizeof (bool), cudaMemcpyDeviceToHost);

    //unsigned int num_neighbors[gpu_stat->host_parameter.num_solvent_mol];
    //cudaMemcpy(num_neighbors, gpu_stat->dev_pl_short.num_neighbors, gpu_stat->host_parameter.num_solvent_mol*sizeof (unsigned int), cudaMemcpyDeviceToHost);

    // DEBUGGING
    //for(unsigned int i = 0; i < gpu_stat->host_parameter.num_solvent_mol; ++i) {
    //  if (gpu_stat->host_parameter.gpu_id==0)
    //  DEBUG(15,"nn " << i << ":" << num_neighbors[i]);
    //}

    // guard for overflow
    if (overflow_short) {
      overflow = true;
      // add 20% more space.
      DEBUG(1,"short overflow");
      DEBUG(1,"max_size = " << gpu_stat->dev_pl_short.max_size);
      const unsigned int new_estimate = gpu_stat->dev_pl_short.max_size + int(gpu_stat->dev_pl_short.max_size * 0.2f);
      free_pairlist(gpu_stat->dev_pl_short);
      allocate_pairlist(gpu_stat->dev_pl_short, gpu_stat->host_parameter.num_solvent_mol, new_estimate);
    }
    if (overflow_long) {
      overflow = true;
      // add 20% more space.
      DEBUG(1,"long overflow");
      DEBUG(1,"max_size = " << gpu_stat->dev_pl_long.max_size);
      const unsigned int new_estimate = gpu_stat->dev_pl_long.max_size + int(gpu_stat->dev_pl_long.max_size * 0.2f);
      free_pairlist(gpu_stat->dev_pl_long);
      allocate_pairlist(gpu_stat->dev_pl_long, gpu_stat->host_parameter.num_solvent_mol, new_estimate);
    }

    // warn the user that this is an issue
    if (overflow) {
      std::cout << "CUDA: Overflow. Recalculating pairlist. This is a performance issue "
              "increase size estimate." << std::endl;
    }
  } while (overflow); // recalculate the pairlist
  
}


__global__ void cukernel::kernel_CalcPairlist
(
        cukernel::simulation_parameter * dev_params,
        float3 * dev_pos,
        pairlist pl_short,
        pairlist pl_long,
        unsigned int num_of_gpus,
        unsigned int gpu_id
) {

  unsigned int num_neighbors_long = 0, num_neighbors_short = 0;
  __shared__ float shared_pos[NUM_THREADS_PER_BLOCK * 3];

  // take host_parameter local
  const unsigned int &N = device_param.num_atoms.solvent;
  const unsigned int &num_solvent_mol = device_param.num_solvent_mol;
  const float &cutoff_long_2 = device_param.cutoff_long_2;
  const float &cutoff_short_2 = device_param.cutoff_short_2;
  //box edges
  const float &box_x = device_param.box.full.x;
  const float &box_y = device_param.box.full.y;
  const float &box_z = device_param.box.full.z;
  
  const float &box_inv_x = device_param.box.inv.x;
  const float &box_inv_y = device_param.box.inv.y;
  const float &box_inv_z = device_param.box.inv.z;
  const unsigned int solvent_offset = device_param.num_atoms_per_mol;

  // calculate indices
  const unsigned int index = blockIdx.x * NUM_THREADS_PER_BLOCK + threadIdx.x;
  const unsigned int molecule_index = index * num_of_gpus + gpu_id;

  const unsigned int first_atom_index = molecule_index*solvent_offset;
  const unsigned int myThreadOffset = threadIdx.x*solvent_offset;

  float3 first_atom_pos;
  if (first_atom_index < N)
    first_atom_pos = dev_pos[first_atom_index];

  for (unsigned int i = 0; i < N; i += (NUM_THREADS_PER_BLOCK * solvent_offset)) {
    float3 neighbor_pos;
    if (i + myThreadOffset < N)
      neighbor_pos = dev_pos[i + myThreadOffset];

    // cache a block of positions
    __syncthreads();
    shared_pos[threadIdx.x] = neighbor_pos.x;
    shared_pos[threadIdx.x + NUM_THREADS_PER_BLOCK] = neighbor_pos.y;
    shared_pos[threadIdx.x + 2 * NUM_THREADS_PER_BLOCK] = neighbor_pos.z;
    __syncthreads();

    unsigned int end_i_loop = NUM_THREADS_PER_BLOCK;
    if (end_i_loop > (N - i) / solvent_offset)
      end_i_loop = (N - i) / solvent_offset;

    if (first_atom_index < N) {
      for (unsigned int start_i_loop = 0; start_i_loop < end_i_loop; start_i_loop++) {
        const unsigned int current_first_atom_index = i + start_i_loop*solvent_offset;
        if (current_first_atom_index != first_atom_index && current_first_atom_index < N) {
          //{ calculate distance and d^2
          float3 dist;
          dist.x = first_atom_pos.x - shared_pos[start_i_loop];
          dist.x -= box_x * rintf(dist.x * box_inv_x);
          dist.y = first_atom_pos.y - shared_pos[start_i_loop + NUM_THREADS_PER_BLOCK];
          dist.y -= box_y * rintf(dist.y * box_inv_y);
          dist.z = first_atom_pos.z - shared_pos[start_i_loop + 2 * NUM_THREADS_PER_BLOCK];
          dist.z -= box_z * rintf(dist.z * box_inv_z);
          const float d2 = abs2(dist);
          //} calculate distance and d^2
       // are they interacting?
       if (d2 < cutoff_long_2) {
            // longrange?
            if (d2 > cutoff_short_2) {
              if (num_neighbors_long < pl_long.max_size) {
                pl_long.list[index + pl_long.pitch * num_neighbors_long] = current_first_atom_index;
                num_neighbors_long++;
              } else {
                *pl_long.overflow = true;
              } // overflow
            } else { // shortrange then
              if (num_neighbors_short < pl_short.max_size) {
                pl_short.list[index + pl_short.pitch * num_neighbors_short] = current_first_atom_index;
                num_neighbors_short++;
              } else {
                *pl_short.overflow = true;
              } // overflow
            } // if shortrange / longrange
          } // if cutoff
        } // if atom in valid range
      } // for atoms j
    } // if atom in valid range
  } // for atoms i
  if (molecule_index < num_solvent_mol) {
    pl_long.num_neighbors[index] = num_neighbors_long;
    pl_short.num_neighbors[index] = num_neighbors_short;
  }
}


