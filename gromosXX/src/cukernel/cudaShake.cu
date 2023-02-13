/**
 * @file cudaShake.cu
 * implementation of m_shake algorithm
 */

#define NUM_THREADS_PER_BLOCK_SHAKE 64

#include "gpu_status.h"
#include "lib/math.h"

#include "cudaShake.h"

#include "lib/utils.h"
#include "parameter.h"

#include "macros.h"

#include "../util/debug.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE cuda


/**
 * Initialize the GPU and create a gpu_status
 */
extern "C" gpu_status * cudaInitGPU_Shake(
        int & device_number,
        double * constr,
        double * factor, double* mass,
        double tol,
        unsigned int num_gpus, unsigned int gpu_id,
        unsigned int num_atoms,
        unsigned int num_solvent_mol,
        int * error) {

  if (device_number != -1) {
    cudaSetDevice(device_number);
  } else {
    cudaGetDevice(&device_number);
    ;;
  }

  // let's first query the device properties and print them out
  DEBUG(4,"Set device properties")
  cudaDeviceProp devProp;
  cudaGetDeviceProperties(&devProp, device_number);
  unsigned flags = 0;
  cudaGetDeviceFlags(&flags);

  if (flags == 0 && cudaSetDeviceFlags(cudaDeviceScheduleYield) == cudaErrorSetOnActiveProcess)
    std::cerr << "Cannot set flags\n";

  //Infos
  std::cout << "CUDA for M_SHAKE" << std::endl;
  std::cout << "\tDeviceproperties:" << std::endl;
  std::cout << "\t\tNumber: " << device_number << std::endl;
  std::cout << "\t\tName: " << devProp.name << std::endl;
  std::cout << "\t\tTotal memory: " << devProp.totalGlobalMem << std::endl;
  std::cout << "\t\tShared memory per block: " << devProp.sharedMemPerBlock << std::endl;
  std::cout << "\t\tTotal constant memory: " << devProp.totalConstMem << std::endl;

  // Create gpu_status
  gpu_status * gpu_stat;
  gpu_stat = (gpu_status *) malloc(sizeof (gpu_status));

  // Assign the numbers of gpus and the gpu_id to the struct
  DEBUG(4,"GPU_SHAKE : " <<gpu_id << ". GPU of " << num_gpus)
  gpu_stat->host_parameter.num_of_gpus = num_gpus;
  gpu_stat->host_parameter.gpu_id = gpu_id;

  DEBUG(4,"Number of solvent mol: " << num_solvent_mol)
  gpu_stat->host_parameter.num_atoms = num_atoms;
  gpu_stat->host_parameter.num_solvent_mol = num_solvent_mol;
  
  // Allocate space for the old and new positions
  HOST_NEW_POS = (VECTOR *) malloc((num_atoms / num_gpus + 3) * sizeof(VECTOR));
  HOST_OLD_POS = (VECTOR *) malloc((num_atoms / num_gpus + 3) * sizeof(VECTOR));
  
  // Constraint lengths
  VECTOR length2;
  length2.x = (FL_PT_NUM) constr[0];
  length2.y = (FL_PT_NUM) constr[1];
  length2.z = (FL_PT_NUM) constr[2];

  cudaMalloc((void**) & DEV_CONST_LENGTH2, sizeof (VECTOR));
  cudaMemcpy(DEV_CONST_LENGTH2, &length2, sizeof (VECTOR), cudaMemcpyHostToDevice);
  DEBUG(10,"Allocated space for const_length2");

  MATRIX factorm;

  factorm.xx = (FL_PT_NUM) factor[0];
  factorm.xy = (FL_PT_NUM) factor[1];
  factorm.xz = (FL_PT_NUM) factor[2];
  factorm.yx = (FL_PT_NUM) factor[3];
  factorm.yy = (FL_PT_NUM) factor[4];
  factorm.yz = (FL_PT_NUM) factor[5];
  factorm.zx = (FL_PT_NUM) factor[6];
  factorm.zy = (FL_PT_NUM) factor[7];
  factorm.zz = (FL_PT_NUM) factor[8];

  cudaMalloc((void**) & DEV_FACTOR, sizeof(MATRIX));
  cudaMemcpy(DEV_FACTOR, & factorm, sizeof(MATRIX), cudaMemcpyHostToDevice);

  cudaMalloc((void**) & DEV_MASS, sizeof(VECTOR));
  VECTOR massv;
  massv.x = FL_PT_NUM (mass[0]);
  massv.y = FL_PT_NUM (mass[1]);
  massv.z = FL_PT_NUM (mass[2]);
  cudaMemcpy(DEV_MASS,& massv, sizeof(VECTOR), cudaMemcpyHostToDevice);

  cudaMalloc((void**) & DEV_TOL, sizeof (FL_PT_NUM));
  FL_PT_NUM tolf = (FL_PT_NUM) tol;
  cudaMemcpy(DEV_TOL, &tolf, sizeof (FL_PT_NUM), cudaMemcpyHostToDevice);
  DEBUG(10,"Allocated space for tol");
  

  cudaMalloc((void**) & gpu_stat->dev_shake_fail_mol, sizeof (int));
  cudaMemset(gpu_stat->dev_shake_fail_mol, -1, sizeof (int));
  DEBUG(10,"Allocated space for shake fail molecule");

  // This index denotes how many molecules will be treated
  // cudaMalloc((void**) & gpu_stat->dev_highest_index, sizeof(unsigned int));

  cudaMalloc((void**) & gpu_stat->dev_parameter, sizeof (cudakernel::simulation_parameter));
  cudaMemcpy(gpu_stat->dev_parameter, &gpu_stat->host_parameter, sizeof (cudakernel::simulation_parameter), cudaMemcpyHostToDevice);
  DEBUG(10,"Allocated space for parameters");
  cudaMalloc((void**) & DEV_NEW_POS, (num_atoms / num_gpus + 3) * sizeof (VECTOR));
  cudaMemset(DEV_NEW_POS, 0, (num_atoms / num_gpus + 3) * sizeof (VECTOR));
  cudaMalloc((void**) & DEV_OLD_POS, (num_atoms / num_gpus + 3) * sizeof (VECTOR));
  cudaMemset(DEV_OLD_POS, 0, (num_atoms / num_gpus + 3) * sizeof (VECTOR));
  DEBUG(10,"Allocated space for positions");

  *error += cudakernel::checkError("after allocating Memory for the old and new positions");

  DEBUG(10,"Return gpu_stat")
  return gpu_stat;
}

/**
 * Apply M_SHAKE algorithm
 */
extern "C" int cudaGPU_Shake(double *newpos, double *oldpos, int & shake_fail_mol,
        gpu_status * gpu_stat) {

  const unsigned int gpu_id = gpu_stat->host_parameter.gpu_id;
  const unsigned int num_atoms = gpu_stat->host_parameter.num_atoms;
  const unsigned int num_gpus = gpu_stat->host_parameter.num_of_gpus;
  const unsigned int num_solvent_mol = gpu_stat->host_parameter.num_solvent_mol;
  const unsigned int first = (num_solvent_mol / num_gpus) * 3 * gpu_id;
  unsigned int last;
  if (gpu_id == num_gpus - 1){
    last = num_atoms;
  }
  else {
    last = (num_solvent_mol / num_gpus) * 3 * (gpu_id + 1);
  }
  DEBUG(15,"GPU : " << gpu_id <<  " First index : " << first << ", last : " << last << " of " << num_atoms)

  // copy the new positions
  DEBUG(10,"Copy the new positions")
  double3 * positions = (double3*) newpos;
  for (unsigned int i = first, j = 0; i < last; i++, j++) {
    HOST_NEW_POS[j].x = (FL_PT_NUM) positions[i].x;
    //if (gpu_id==1)
    DEBUG(15,"Old Pos : i " << i << " x : " << HOST_NEW_POS[j].x)
    HOST_NEW_POS[j].y = (FL_PT_NUM) positions[i].y;
    HOST_NEW_POS[j].z = (FL_PT_NUM) positions[i].z;
  }


  // Copy the old positions
  DEBUG(10,"Copy the old positions")
  double3 * old_positions = (double3*) oldpos;
  for (unsigned int i = first, j = 0; i < last; i++, j++) {
    HOST_OLD_POS[j].x = (FL_PT_NUM) old_positions[i].x;
    HOST_OLD_POS[j].y = (FL_PT_NUM) old_positions[i].y;
    HOST_OLD_POS[j].z = (FL_PT_NUM) old_positions[i].z;
  }

  // Copy the old and new positions to the GPU
  DEBUG(10,"Copy the positions to the GPU")
  cudaMemcpy(DEV_NEW_POS, HOST_NEW_POS, (num_atoms / num_gpus + 3) * sizeof (VECTOR), cudaMemcpyHostToDevice);
  cudaMemcpy(DEV_OLD_POS, HOST_OLD_POS, (num_atoms / num_gpus + 3) * sizeof (VECTOR), cudaMemcpyHostToDevice);
  cudakernel::checkError("after copying the positions");

  // Copy the highest index of the molecule

  const unsigned int highest_index =  (last - first) / 3;
  DEBUG(10,"Higest molecule index : " << highest_index)
  //cudaMemset(gpu_stat->dev_highest_index,, sizeof(unsigned int));

  // Dimensions
  unsigned int numBlocks = (unsigned int) (num_solvent_mol / num_gpus + 1) / NUM_THREADS_PER_BLOCK_SHAKE + 1;
  dim3 dimGrid(numBlocks, 1);
  dim3 dimBlock(NUM_THREADS_PER_BLOCK_SHAKE, 1);
  DEBUG(10,"numBlocks: " << numBlocks)

  DEBUG(7,"Starting kernel")
  cudakernel::kernel_Calc_Shake <<<dimGrid, dimBlock >>>
          (DEV_NEW_POS, DEV_OLD_POS, gpu_stat->dev_parameter, gpu_stat->dev_shake_fail_mol,
          DEV_TOL, DEV_MASS, DEV_CONST_LENGTH2, DEV_FACTOR, highest_index);
  cudakernel::checkError("after GPU_SHAKE");


  // Copy the new positions from the GPU
  DEBUG(10,"Get the new positions")
  cudaMemcpy(HOST_NEW_POS, DEV_NEW_POS, (num_atoms / num_gpus + 3) * sizeof (VECTOR), cudaMemcpyDeviceToHost);
  cudakernel::checkError("after copying the new positions");

  for (unsigned int i = first, j = 0; i < last; ++i, j++) {
    positions[i].x = (double) HOST_NEW_POS[j].x;
    //if(gpu_id == 1)
    DEBUG(15,"GPU : " << gpu_id << " New Pos i: " << i << " x: " << positions[i].x)
    positions[i].y = (double) HOST_NEW_POS[j].y;
    positions[i].z = (double) HOST_NEW_POS[j].z;
  }

  // Check, if everything went well
  cudaMemcpy(&shake_fail_mol, gpu_stat->dev_shake_fail_mol, sizeof (int), cudaMemcpyDeviceToHost);
  DEBUG(7,"GPU : " << gpu_id << " Get fail molecule, which is " << shake_fail_mol)
  if (shake_fail_mol > 0){
      shake_fail_mol = shake_fail_mol + first / 3;
      return 1;
  }

  return 0;
}

__global__ void cudakernel::kernel_Calc_Shake
(
        VECTOR * new_pos, VECTOR * old_pos,
        cudakernel::simulation_parameter * dev_params,
        int *shake_fail_mol, FL_PT_NUM * tol,
        VECTOR * mass, VECTOR * const_length2,
        MATRIX * factor, unsigned int highest_mol_index
        ) {

  const unsigned int mol_index = blockIdx.x * NUM_THREADS_PER_BLOCK_SHAKE + threadIdx.x;
  //const unsigned int shared_index = blockIdx.x;


  if (mol_index >= highest_mol_index)
    return;

  // Shared memory for faster access
  __shared__ VECTOR pos_shr[3][NUM_THREADS_PER_BLOCK_SHAKE];
  pos_shr[0][threadIdx.x] = new_pos[mol_index * 3];
  pos_shr[1][threadIdx.x] = new_pos[mol_index * 3 + 1];
  pos_shr[2][threadIdx.x] = new_pos[mol_index * 3 + 2];

  const VECTOR old_pos_1 = old_pos[mol_index * 3];
  const VECTOR old_pos_2 = old_pos[mol_index * 3 + 1];
  const VECTOR old_pos_3 = old_pos[mol_index * 3 + 2];

  const VECTOR dist_old_1 = old_pos_1 - old_pos_2;
  const VECTOR dist_old_2 = old_pos_1 - old_pos_3;
  const VECTOR dist_old_3 = old_pos_2 - old_pos_3;

  const FL_PT_NUM ltol = *tol * PREC(2.0);
  const VECTOR cl2 = *const_length2;
  const MATRIX lfactor = *factor;

  bool convergence = false;
  for (unsigned int i = 0;  !convergence; ++i) {
    const VECTOR dist_new_1 = pos_shr[0][threadIdx.x] - pos_shr[1][threadIdx.x];
    const VECTOR dist_new_2 = pos_shr[0][threadIdx.x] - pos_shr[2][threadIdx.x];
    const VECTOR dist_new_3 = pos_shr[1][threadIdx.x] - pos_shr[2][threadIdx.x];

    const VECTOR dist = MAKE_VECTOR(abs2(dist_new_1), abs2(dist_new_2), abs2(dist_new_3));
    const VECTOR diff = cl2 - dist;
    
    convergence = true;

    if (FABS(diff.x) >= cl2.x * ltol  ||
            FABS(diff.y) >= cl2.y * ltol ||
            FABS(diff.z) >= cl2.z * ltol) {

      convergence = false;

      // matrix A
      MATRIX A;
      A.xx = dot(dist_old_1, dist_new_1) * lfactor.xx;
      A.xy = dot(dist_old_2, dist_new_1) * lfactor.xy;
      A.xz = dot(dist_old_3, dist_new_1) * lfactor.xz;
      A.yx = dot(dist_old_1, dist_new_2) * lfactor.yx;
      A.yy = dot(dist_old_2, dist_new_2) * lfactor.yy;
      A.yz = dot(dist_old_3, dist_new_2) * lfactor.yz;
      A.zx = dot(dist_old_1, dist_new_3) * lfactor.zx;
      A.zy = dot(dist_old_2, dist_new_3) * lfactor.zy;
      A.zz = dot(dist_old_3, dist_new_3) * lfactor.zz;

      // vectors orthogonal? -> SHAKE error
      
      if (A.xx < cl2.x * EPSILOND ||
              A.yy < cl2.y * EPSILOND ||
              A.zz < cl2.z * EPSILOND || i > 1000) {
        *shake_fail_mol = int(mol_index) + 1;
        return;
      }
       
      // inverse

      const FL_PT_NUM det =
                A.xy * A.yz * A.zx
              - A.xz * A.yy * A.zx
              + A.xz * A.yx * A.zy
              - A.xx * A.yz * A.zy
              - A.xy * A.yx * A.zz
              + A.xx * A.yy * A.zz;
      const FL_PT_NUM deti = PREC(1.0) / det;
      MATRIX Ai;
      Ai.xx = (A.yy * A.zz - A.yz * A.zy) * deti;
      Ai.yx = (A.yz * A.zx - A.yx * A.zz) * deti;
      Ai.zx = (A.yx * A.zy - A.yy * A.zx) * deti;
      Ai.xy = (A.xz * A.zy - A.xy * A.zz) * deti;
      Ai.yy = (A.xx * A.zz - A.xz * A.zx) * deti;
      Ai.zy = (A.xy * A.zx - A.xx * A.zy) * deti;
      Ai.xz = (A.xy * A.yz - A.xz * A.yy) * deti;
      Ai.yz = (A.xz * A.yx - A.xx * A.yz) * deti;
      Ai.zz = (A.xx * A.yy - A.xy * A.yx) * deti;

      const FL_PT_NUM f0 = (Ai.xx * diff.x + Ai.xy * diff.y + Ai.xz * diff.z) * PREC(0.5);
      const FL_PT_NUM f1 = (Ai.yx * diff.x + Ai.yy * diff.y + Ai.yz * diff.z) * PREC(0.5);
      const FL_PT_NUM f2 = (Ai.zx * diff.x + Ai.zy * diff.y + Ai.zz * diff.z) * PREC(0.5);

      const VECTOR f01 = f0 * dist_old_1 + f1 * dist_old_2;
      const VECTOR f02 = f2 * dist_old_3 - f0 * dist_old_1;
      const VECTOR f12 = f1 * dist_old_2 + f2 * dist_old_3;

      // calculate the new constrained positions
      pos_shr[0][threadIdx.x] = pos_shr[0][threadIdx.x] + f01 / mass->x;
      pos_shr[1][threadIdx.x] = pos_shr[1][threadIdx.x] + f02 / mass->y;
      pos_shr[2][threadIdx.x] = pos_shr[2][threadIdx.x] - f12 / mass->z;

    } // if we have to shake
  } // while not convergence

  new_pos[mol_index * 3] = pos_shr[0][threadIdx.x];
  new_pos[mol_index * 3 + 1] = pos_shr[1][threadIdx.x];
  new_pos[mol_index * 3 + 2] = pos_shr[2][threadIdx.x];
}

#undef DEBUG


