/**
 * @file pairlist.h
 * pairlist computation
 */
#ifndef CUKERNEL_PAIRLIST
#define CUKERNEL_PAIRLIST
namespace cudakernel {
/**
 * @struct pairlist
 *  the holds the pairlist
 */
struct pairlist {
  /**
   * the elements in the pairlist [i*pitch + j]
   */
  unsigned int *list;
  /**
   * the number of neighbors in the list
   */
  unsigned int *num_neighbors;
  /**
   * the maximal number of neighbors that can be stored in the list
   */
  unsigned int max_size;
  /**
   * the pitch (from cudaMallocPitch)
   */
  unsigned int pitch;
  /**
   * bool overflow
   */
  bool *overflow;
};

#ifndef CUKERNEL_TYPES_ONLY
/**
 * free a pairlist object
 * @param pl the pairlist you want to free on the device
 */
void free_pairlist(pairlist &pl);
/**
 * allocate a pairlist object on the device
 * @param pl the pl object to hold the pointers
 * @param size the number of solvent molecules
 * @param max_neighbors the maximum number of neighbors the pairlist can hold
 */
void allocate_pairlist(pairlist &pl, unsigned int size, unsigned int max_neighbors);

/**
 * computes the pairlist on the GPU
 * @param[in] dev_param the simulation parameters
 * @param[in] dev_pos the positions
 * @param[out] pl_short the short-range pairlist
 * @param[out] pl_long the long-range pairlist
 */
__global__ void kernel_CalcPairlist(
        cudakernel::simulation_parameter * dev_params,
        float3 * dev_pos,
        pairlist pl_short,
        pairlist pl_long,
        unsigned int num_of_gpus,
        unsigned int gpu_id);
#endif
}
#endif

