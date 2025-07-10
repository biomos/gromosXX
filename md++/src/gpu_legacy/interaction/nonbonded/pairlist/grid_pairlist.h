/**
 * @file grid_pairlist.h
 * @author Poliak (peter.poliak@boku.ac.at)
 * @date 17.06.2023
 * CUDA implementation of grid pairlist
 */

namespace cukernel {
  namespace interaction {
    /**
     * Find pairs between neighbouring cells
     * ordering has no penalty on performance, but unordered we can
     * write more densely (+15% speedup)
     * and do warp-wide reduction (so far seems comparable to shared memory)
     * 
     * shared memory requirements:
     * (num_atoms1 + num_atoms2) * (float3 + unsigned + 2*int)
     */

    /**
     * @brief Find pairs between neighbouring cells
     * ordering has no penalty on performance, but unordered we can
     * write more densely (+15% speedup)
     * and do warp-wide reduction (so far seems comparable to shared memory)
     * 
     * @tparam N size of the datablock to be processed by single block
     * @tparam ordered store the pairlist in ordered manner, i.e. atom1 < atom2
     * @param pos atom coordinates
     * @param cellgrid cellgrid object holding the cells data
     * @param pairlist_long longrange pairlist
     * @param pairlist_short shortrange pairlist
     * @param cell_i index of the main cell
     * @param cell_j index of the neighbouring cell
     */
    template <int N, bool ordered = true>
    __global__ void find_pairs_neighbour(const float3 *pos
                                , const cukernel::Grid cellgrid
                                , Pairlist pairlist_long
                                , Pairlist pairlist_short
                                , unsigned cell_i
                                , unsigned cell_j);
    
    /**
     * find pairs within the grid cell - monolith version
     * rerun this kernel for every cell
     * num_threads should be tuned
     * num_blocks to cover all pairs, one pair per thread
     * 
     * shared memory requirement: num_atoms*[sizeof(float3) + sizeof(unsigned) + sizeof(short int)]
     * 
     * this performs 100s times better than the looping version
     * probably with the huge number of blocks running in parallel,
     * GPU can better hide the memory latency and sync wait?
     * and/or broadcasts more?
     */

    /**
     * @brief find pairs within the grid cell - monolith version
     * rerun this kernel for every cell
     * num_threads should be tuned
     * num_blocks to cover all pairs, one pair per thread
     * 
     * shared memory requirement: num_atoms*[sizeof(float3) + sizeof(unsigned) + sizeof(short int)]
     * 
     * this performs 100s times better than the looping version
     * probably with the huge number of blocks running in parallel,
     * GPU can better hide the memory latency and sync wait?
     * and/or broadcasts more?
     * 
     * @param pos atom coordinates
     * @param cellgrid cellgrid object holding the cells data
     * @param pairlist_long longrange pairlist
     * @param pairlist_short shortrange pairlist
     * @param cell_i index of the cell
     * @param num_pairs number of pairs = (num_atoms * (num_atoms - 1) / 2)
     * @param cell_j index of the neighbouring cell
     */
    __global__ void find_pairs_self(const float3 *pos
                                , const cukernel::Grid cellgrid
                                , Pairlist pairlist_long
                                , Pairlist pairlist_short
                                , unsigned cell_i
                                , unsigned num_pairs // this we know in advance since we do not iterate over cells
                                );
  }
}