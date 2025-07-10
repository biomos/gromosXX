/**
 * @file grid_pairlist.cu
 * grid pairlist computation
 */

/*#include <iostream>
#include "gpu_status.h"

#include "lib/math.h"
#include "../util/debug.h"*/

#undef MODULE
#undef SUBMODULE
#define MODULE cukernel
#define SUBMODULE pairlist

template <int N, bool ordered = true>
__global__ void cukernel::interaction::find_pairs_neighbour(
                              const float3 *pos
                            , const cukernel::Grid cellgrid
                            , Pairlist pairlist_long
                            , Pairlist pairlist_short
                            , unsigned cell_i
                            , unsigned cell_j
                            , float3 fpshift) {
  // we need to store positions and sizes of atoms 1 and 2
  extern __shared__ char array[];
  cooperative_groups::thread_block tb = cooperative_groups::this_thread_block();
  dim3 tidx = tb.thread_index();
  dim3 bidx = tb.group_index();
  dim3 bdim = tb.group_dim();
  tidx.x += N * bdim.x * bidx.x;
  tidx.y += N * bdim.y * bidx.y;
  const unsigned lidx = threadIdx.x + threadIdx.y * blockDim.x;
  const unsigned bsize = blockDim.x * blockDim.y;
  // number of atoms processed by this block
  const unsigned num_atoms1 = N * blockDim.x;
  const unsigned num_atoms2 = N * blockDim.y;
  float3 *spos1 = reinterpret_cast<float3*>(&array[0]); // size = num_atoms1 * sizeof(float3);
  float3 *spos2 = reinterpret_cast<float3*>(&spos1[num_atoms1]); // size = num_atoms2 * sizeof(float3);
  unsigned *indices1 = reinterpret_cast<unsigned*>(&spos2[num_atoms2]); // size = num_atoms1 * sizeof(unsigned);
  unsigned *indices2 = reinterpret_cast<unsigned*>(&indices1[num_atoms1]); // size = num_atoms2 * sizeof(unsigned);
  
  // sizes of the long pairlist
  int *sizes_long1 = reinterpret_cast<int*>(&indices2[num_atoms2]); // size = num_atoms1 * sizeof(int);
  int *sizes_long2 = reinterpret_cast<int*>(&sizes_long1[num_atoms1]); // size = num_atoms2 * sizeof(int);
  // sizes of the short pairlist
  int *sizes_short1 = reinterpret_cast<int*>(&sizes_long2[num_atoms2]); // size = num_atoms1 * sizeof(int);
  int *sizes_short2 = reinterpret_cast<int*>(&sizes_short1[num_atoms1]); // size = num_atoms2 * sizeof(int);
  // offset to get from atoms1 to atoms2
  const unsigned offset2 = num_atoms1;
  
  /*ushort2 *pairs = reinterpret_cast<ushort2*>(&sizes2[num_atoms2]); // size = 4*blocksize
  __shared__ unsigned pairs_size;
  if (lidx == 0) pairs_size = 0;*/

  // total number of atoms in the cells
  const unsigned total_num_atoms1 = cellgrid.cell_size(cell_i);
  const unsigned total_num_atoms2 = cellgrid.cell_size(cell_j);

  // initialize sizes arrays and copy positions to shared memory - we should load only data for this block
  for (unsigned i = lidx; i < num_atoms1; i += bsize) {
    const unsigned gidx = N * bidx.x * bdim.x + i;
    if (gidx < total_num_atoms1) {
      const unsigned a1 = cellgrid.get_atom_in_cell(gidx, cell_i);
      indices1[i] = a1;
      spos1[i] = pos[a1];
      sizes_long1[i] = 0;
      sizes_short1[i] = 0;
      //printf("bidx: %u, lidx: %u, i: %u, gidx: %u, a1: %u\n", bidx.x, lidx, i, gidx, a1);
    }
  }
  for (unsigned j = lidx; j < num_atoms2; j += bsize) {
    const unsigned gidx = N * bidx.y * bdim.y + j;
    if (gidx < total_num_atoms2) {
      const unsigned a2 = cellgrid.get_atom_in_cell(gidx, cell_j);
      indices2[j] = a2;
      spos2[j] = pos[a2]/* + fpshift*/;
      sizes_long2[j] = 0;
      sizes_short2[j] = 0;
      //printf("bidx: %u, lidx: %u, j: %u, gidx: %u, a2: %u\n", bidx.y, lidx, j, gidx, a2);
    }
  }

  tb.sync();
  
  unsigned *target[N][N];
  int *sizes_ptr[N][N];
  unsigned aidx[N][N];
  unsigned value[N][N];
  for (unsigned si = 0; si < N; ++si) {
    for (unsigned sj = 0; sj < N; ++sj) {
      target[si][sj] = nullptr;
      sizes_ptr[si][sj] = nullptr;
      aidx[si][sj] = 0;
      value[si][sj] = 0;
      const unsigned tx = tidx.x + si * bdim.x;
      const unsigned ty = tidx.y + sj * bdim.y;
      if (tx < total_num_atoms1 && ty < total_num_atoms2) {
        const unsigned i = threadIdx.x + si * bdim.x;
        const unsigned j = threadIdx.y + sj * bdim.y;
        const unsigned a1 = indices1[i];
        const unsigned a2 = indices2[j];
        //const float3 v_nim = spos1[i] - spos2[j];
        const float3 v_nim = nim(spos1[i], spos2[j]);
        //printf("%u %u : v_nim: %f %f %f\n", a1, a2, v_nim.x, v_nim.y, v_nim.z);
        float dist2 = abs2(v_nim);
        if (dist2 < dev_cutoff2_long) {
          const bool shortrange = dist2 < dev_cutoff2_short;
          int *sptr = nullptr;
          unsigned ai = 0;
          unsigned a = 0;
          unsigned val = 0;
          sptr = shortrange ? sizes_short1 : sizes_long1;
          ai = i;
          a = a1;
          val = a2;
          if (ordered && a1 > a2) {
            sptr += offset2;
            ai = j;
            a = a2;
            val = a1;
          }
          // To do warp-wide reduction here, we need to give up ordered pairlist, i.e. a1 < a2
          unsigned my_offset = 0;
          // if (!ordered) {
          //   /*** START WARP WIDE ***/ //
          //   // no improvement over shared writes observed
            
          //   // mask of active threads within warp
          //   int active = __activemask();
          //   // elect a leader
          //   int leader = __ffs(active) - 1;
            
          //   const unsigned warp_atom = j;
          //   //todo: do this warp wide
          //   unsigned my_offset;
          //   if (threadIdx.x == leader) {
          //     // leader writes total to shared memory
          //     my_offset = sizes1[warp_atom];
          //     sizes1[warp_atom] += __popc(active);
          //   }
          //   // broadcast offset and calculate my offset from the mask - shift away higher threads and count '1' bits
          //   my_offset = __popc(active << (32 - threadIdx.x)) + __shfl_sync(active, my_offset, leader);
          // /*** END WARP WIDE ***/
          // } else {
            my_offset = atomicAdd(&sptr[ai], 1);
          //}
          Pairlist *pairlist = shortrange ? &pairlist_short : &pairlist_long;

          target[si][sj] = (*pairlist)(a) + my_offset;
          aidx[si][sj] = ai;
          sizes_ptr[si][sj] = sptr;
          value[si][sj] = val;
        }
      }
    }
  }
  tb.sync();
  // from now we do not need spos

  // load global offsets from global to shared memory - longrange
  for (unsigned i = lidx; i < num_atoms1; i += bsize) {
    const unsigned gidx = N * bidx.x * bdim.x + i;
    if (gidx < total_num_atoms1) {
      const unsigned a1 = indices1[i];
      //printf("gidx: %u, a1: %u, sizes_long1[%u]: %u\n", gidx, a1, i, sizes_long1[i]);
        const int s = sizes_long1[i];
      if (s > 0) { // save some global reads
        const int global_offset_long1 = pairlist_long.reserve_strip<false>(a1, s);
        sizes_long1[i] = global_offset_long1;
      }
    }
  }
  if (ordered) { // but here is the speedup of unordered ca. 15%, as we do not need to write to both cells !
    for (unsigned j = lidx; j < num_atoms2; j += bsize) {
      const unsigned gidx = N * bidx.y * bdim.y + j;
      if (gidx < total_num_atoms2) {
        const unsigned a2 = indices2[j];
        //printf("gidx: %u, a2: %u, sizes_long2[%u]: %u\n", gidx, a2, j, sizes_long2[j]);
        const int s = sizes_long2[j];
        if (s > 0) {
          const int global_offset_long2 = pairlist_long.reserve_strip<false>(a2, s);
          sizes_long2[j] = global_offset_long2;
        }
      }
    }
  }

  // load global offsets from global to shared memory - shortrange
  for (unsigned i = lidx; i < num_atoms1; i += bsize) {
    const unsigned gidx = N * bidx.x * bdim.x + i;
    if (gidx < total_num_atoms1) {
      const unsigned a1 = indices1[i];
      //printf("gidx: %u, a1: %u, sizes_short1[%u]: %u\n", gidx, a1, i, sizes_short1[i]);
      const int s = sizes_short1[i];
      if (s > 0) { // save some global reads
        const int global_offset_short1 = pairlist_short.reserve_strip<false>(a1, s);
        sizes_short1[i] = global_offset_short1;
      }
    }
  }
  if (ordered) { // but here is the speedup of unordered ca. 15% !
    for (unsigned j = lidx; j < num_atoms2; j += bsize) {
      const unsigned gidx = N * bidx.y * bdim.y + j;
      if (gidx < total_num_atoms2) {
        const unsigned a2 = indices2[j];
        //printf("gidx: %u, a2: %u, sizes_short2[%u]: %u\n", gidx, a2, j, sizes_short2[j]);
        const int s = sizes_short2[j];
        if (s > 0) {
          const int global_offset2 = pairlist_short.reserve_strip<false>(a2, s);
          sizes_short2[j] = global_offset2;
        }
      }
    }
  }
  tb.sync();
  //assert(!pairlist.overflown());
  for (unsigned si = 0; si < N; ++si) {
    for (unsigned sj = 0; sj < N; ++sj) {
      if (target[si][sj] != nullptr) {
        // add global offset
        const int *const sptr = sizes_ptr[si][sj];
        unsigned* tar = target[si][sj];
        const unsigned ai = aidx[si][sj];
        if (sptr[ai] != -1) { // if not overflown
            tar += sptr[ai];
            *(tar) = value[si][sj];
        }
      }
    }
  }
}


__global__ void cukernel::interaction::find_pairs_self(
                              const float3 *pos
                            , const cukernel::Grid cellgrid
                            , Pairlist pairlist_long
                            , Pairlist pairlist_short
                            , unsigned cell_i
                            , unsigned num_pairs // this we know in advance since we do not iterate over cells
                            ) {
    // here we also have to consider atom positions and load them to shared memory, as we will be reusing them
    extern __shared__ char array[];

    const unsigned lidx = threadIdx.x;
    const unsigned idx = threadIdx.x + blockDim.x*blockIdx.x;
    //const unsigned bidx = blockIdx.x; // for best reuse of shared memory, one block should iterate over the entire cell
    const unsigned blocksize = blockDim.x;
    //const unsigned gridsize = gridDim.x;

    const unsigned num_atoms = cellgrid.cell_size(cell_i);

    float3* spos = reinterpret_cast<float3*>(array); // size = num_atoms * sizeof(float3);
    unsigned* indices = reinterpret_cast<unsigned*>(&spos[num_atoms]); // size = num_atoms * sizeof(unsigned);
    int* sizes_long = reinterpret_cast<int*>(&indices[num_atoms]); // size = num_atoms * sizeof(int);
    int* sizes_short = reinterpret_cast<int*>(&sizes_long[num_atoms]); // size = num_atoms * sizeof(int);
    // should we maybe store pairlist in shared memory as bitmask? that would be still 125kB for 1000 atoms
    //unsigned* boolmatrix = reinterpret_cast<int*>(&indices[num_atoms]); // size = num_atoms * sizeof(int);
    /*__shared__ unsigned active_threads_long;
    __shared__ unsigned active_threads_short;
    if (lidx == 0) active_threads_long = 0;
    if (lidx == 0) active_threads_short = 0;*/

    //const unsigned num_pairs = num_atoms * (num_atoms - 1) / 2; // provided in function call is faster, but reading from global for every block is slower
    // load atom positions and initialize their pairlist sizes in shared memory
    for (unsigned i = lidx; i < num_atoms; i += blocksize) {
        const unsigned a1 = cellgrid.get_atom_in_cell(i, cell_i);
        spos[i] = pos[a1];
        sizes_long[i] = 0;
        sizes_short[i] = 0;
        indices[i] = a1;
        //printf("cell_i: %u, num_atoms:%u , bidx: %u, sizes[%u]: %d, indices[%u]: %u\n", cell_i, num_atoms, bidx, i, sizes[i], i, indices[i]);
    }
    __syncthreads();
    /*for (unsigned i = lidx; i < num_atoms; i += blocksize) {
        printf("cell_i: %u, num_atoms:%u , bidx: %u, sizes[%u]: %d, indices[%u]: %u\n", cell_i, num_atoms, bidx, i, sizes[i], i, indices[i]);
    }
    __syncthreads();*/
    // monolith, just get your own pair, mate
    const unsigned pair_index = idx;
    const unsigned last_id = num_pairs - 1;
    //int my_offset = 0;
    unsigned *target = nullptr;
    unsigned aidx = 0;
    unsigned a = 0;
    unsigned value = 0;
    int *sizes = nullptr;
    bool shortrange = false;
    if (pair_index < num_pairs) {
        //const unsigned i = num_atoms - 2 - ((unsigned)sqrt((float)(8 * ((unsigned)last_id - (unsigned)idx) + 1 )) - 1) / 2;
        const unsigned i = num_atoms - 2 - (fast_isqrt(8 * (last_id - pair_index) + 1 ) - 1) / 2;
        const unsigned j = (pair_index + (i * ((i + 1) + 2) / 2 + 1)) % num_atoms;
        const float3 v_nim = nim(spos[i], spos[j]);
        const float dist2 = abs2(v_nim);
        if (dist2 < dev_cutoff2_long) {
            shortrange = dist2 < dev_cutoff2_short;
            //unsigned *_active = shortrange ? &active_threads_short : &active_threads_long;
            //atomicAdd(_active, 1);
            //pairlist = &pairlist_long;
            const unsigned a1 = indices[i];
            const unsigned a2 = indices[j];
            if (a1 < a2) {
                aidx = i;
                a = a1;
                value = a2;
            } else {
                aidx = j;
                a = a2;
                value = a1;
            }
            sizes = shortrange ? sizes_short : sizes_long;
            //sizes = sizes_long;
            const unsigned my_offset = atomicAdd(&sizes[aidx], 1);
            //target = pointer2D(pl, pl_pitch, a, my_offset);
            Pairlist *pairlist = shortrange ? &pairlist_short : &pairlist_long;
            target = (*pairlist)(a) + my_offset;
            //pairlist.push_back(a,value);
        }
    }
    __syncthreads();
    /*if (lidx == 0) {
      if (active_threads_long) {
        printf("cell: %u, block: %u, long ACTIVE (active_threads = %u/%u)\n", cell_i, blockIdx.x, active_threads_long, blockDim.x);
      } else {
        printf("cell: %u, block: %u, long INACTIVE\n", cell_i, blockIdx.x);
      }
      if (active_threads_short) {
        printf("cell: %u, block: %u, short ACTIVE (active_threads = %u/%u)\n", cell_i, blockIdx.x, active_threads_short, blockDim.x);
      } else {
        printf("cell: %u, block: %u, short INACTIVE\n", cell_i, blockIdx.x);
      }
    }*/

    // the block does not write exclusively to global, other blocks operate in this cell as well
    // reserve space in global pairlist and load global offsets to shared memory
    for (unsigned i = lidx; i < num_atoms; i += blocksize) {
        int s = sizes_long[i];
        if (s > 0) {
          int global_offset = pairlist_long.reserve_strip<false>(indices[i], s);
          sizes_long[i] = global_offset;
        }
        s = sizes_short[i];
        if (s > 0) {
          int global_offset = pairlist_short.reserve_strip<false>(indices[i], s);
          sizes_short[i] = global_offset;
        }
    }
    /*for (unsigned i = lidx; i < num_atoms; i += blocksize) {
    }*/
    __syncthreads();
    if (target != nullptr) {
        // offset -1 means global memory overflow
        int offset = sizes[aidx];
        if (offset != -1) {
          // add global offset
          target += offset;
          // write value
          *(target) = value;
        }
    }
}


