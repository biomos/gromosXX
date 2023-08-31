/**
 * @file constraints.cu
 * implementation of constraints algorithms
*/

#include "gpu_status.h"
#include "constraints.h"
#include "parameter.h"
#include "lib/math.h"
#include "lib/utils.h"

#include "../util/debug.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE cuda

#define NUM_THREADS_PER_BLOCK_SETTLE 128



/**
 * inititalze the gpu for the settle calculations
 */
extern "C" gpu_status * cudaInitConstraints(unsigned int num_of_gpus, unsigned int gpu_id,
        unsigned int num_atoms,
        unsigned int num_solvent_mol) {
  
  // let's first query the device properties and print them out
  DEBUG(4,"Set device properties")
  cudaDeviceProp devProp;
  cudaGetDeviceProperties(&devProp, gpu_id);
  cudaSetDevice(gpu_id);
  if (cudaSetDeviceFlags(cudaDeviceScheduleYield) == cudaErrorSetOnActiveProcess)
    std::cerr << "Cannot set flags\n";

  //Infos
  std::cout << "\tCUDA for settle" << std::endl;
  std::cout << "\tDeviceproperties:" << std::endl;
  std::cout << "\t\tName: " << devProp.name << std::endl;
  std::cout << "\t\tTotal memory: " << devProp.totalGlobalMem << std::endl;
  std::cout << "\t\tShared memory per block: " << devProp.sharedMemPerBlock << std::endl;
  std::cout << "\t\tTotal constant memory: " << devProp.totalConstMem << std::endl;

  // Create gpu_status
  gpu_status * gpu_stat;
  gpu_stat = (gpu_status *) malloc(sizeof (gpu_status));

  // Assign the numbers of gpus and the gpu_id to the struct
  gpu_stat->host_parameter.num_of_gpus = num_of_gpus;
  gpu_stat->host_parameter.gpu_id = gpu_id;

  gpu_stat->host_parameter.num_atoms = num_atoms;
  gpu_stat->host_parameter.num_solvent_mol = num_solvent_mol;

  // Allocate space for the old and new positions
  gpu_stat->host_double_new_pos = (double3 *) malloc(num_atoms * sizeof(double3));
  gpu_stat->host_double_old_pos = (double3 *) malloc(num_atoms * sizeof(double3));

  DEBUG(4,"Allocating memory on the GPU for SETTLE")
  cudaMalloc((void**) & gpu_stat->dev_shake_fail_mol, sizeof (int));
  cudaMemset(gpu_stat->dev_shake_fail_mol, 0, sizeof(int));
  cudaMalloc((void**) & gpu_stat->dev_parameter, sizeof (cudakernel::simulation_parameter));
  cudaMemset(gpu_stat->dev_parameter, 0, sizeof(cudakernel::simulation_parameter));
  cudaMalloc((void**) & gpu_stat->dev_double_new_pos, num_atoms * sizeof (double3));
  cudaMemset(gpu_stat->dev_double_new_pos, 0, num_atoms * sizeof (double3));
  cudaMalloc((void**) & gpu_stat->dev_double_old_pos, num_atoms * sizeof (double3));
  cudaMemset(gpu_stat->dev_double_old_pos, 0, num_atoms * sizeof (double3));

  cudakernel::checkError("after allocating Memory for the old and new positions");

  return gpu_stat;
}

extern "C" int cudaConstraints(double *newpos, double *oldpos,
                    int & shake_fail_mol, gpu_status * gpu_stat) {
  
  DEBUG(4,"shake_fail_mol in the beginning is " << shake_fail_mol)
  // copy the new positions
  DEBUG(10,"Copy the new positions")
  double3 * new_positions = (double3*) newpos;
  for (unsigned int i = 0; i < gpu_stat->host_parameter.num_atoms; i++) {
    gpu_stat->host_double_new_pos[i].x = new_positions[i].x;
    gpu_stat->host_double_new_pos[i].y = new_positions[i].y;
    gpu_stat->host_double_new_pos[i].z = new_positions[i].z;
  }

  // Copy the old positions
  DEBUG(10,"Copy the old positions")
  double3 * old_positions = (double3*) oldpos;
  for (unsigned int i = 0; i < gpu_stat->host_parameter.num_atoms; i++) {
    gpu_stat->host_double_old_pos[i].x = old_positions[i].x;
    gpu_stat->host_double_old_pos[i].y = old_positions[i].y;
    gpu_stat->host_double_old_pos[i].z = old_positions[i].z;
  }

  // Copy the new and old positions on the GPU
  DEBUG(10,"Copy the new positions on to the GPU")
  cudaMemcpy(gpu_stat->dev_double_new_pos, gpu_stat->host_double_new_pos, gpu_stat->host_parameter.num_atoms * sizeof (double3), cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_stat->dev_double_old_pos, gpu_stat->host_double_old_pos, gpu_stat->host_parameter.num_atoms * sizeof (double3), cudaMemcpyHostToDevice);

  cudakernel::checkError("after copying the new positions to the GPU");

  // Copy the paramter
  cudaMemcpy(gpu_stat->dev_parameter, & gpu_stat->host_parameter, sizeof(cudakernel::simulation_parameter), cudaMemcpyHostToDevice);


  // Copy the error molecule index
  cudaMemcpy(gpu_stat->dev_shake_fail_mol, &shake_fail_mol, sizeof(int), cudaMemcpyHostToDevice);

  // Dimensions
  unsigned int numBlocks = (unsigned int) gpu_stat->host_parameter.num_solvent_mol / NUM_THREADS_PER_BLOCK_SETTLE + 1;
  dim3 dimGrid(numBlocks, 1);
  dim3 dimBlock(NUM_THREADS_PER_BLOCK_SETTLE, 1);

  DEBUG(4,"Actually calculate the constraints on the GPU. (Fail Mol : " << shake_fail_mol << ")")
  kernel_CalcConstraints_Settle <<<dimGrid, dimBlock>>>
          (gpu_stat->dev_double_new_pos, gpu_stat->dev_double_old_pos,
            gpu_stat->dev_parameter,
            gpu_stat->dev_shake_fail_mol);
  
  cudakernel::checkError("after GPU_SETTLE"); 

  DEBUG(7,"Copy the Error molecule index")
  cudaMemcpy(&shake_fail_mol, gpu_stat->dev_shake_fail_mol, sizeof (int), cudaMemcpyDeviceToHost);
  DEBUG(7,"... which is : " << shake_fail_mol << "\n")
  if (shake_fail_mol >= 0)
    return 1;

    DEBUG(10,"Copy the positions")
  cudaMemcpy(gpu_stat->host_double_new_pos, gpu_stat->dev_double_new_pos, gpu_stat->host_parameter.num_atoms * sizeof (double3), cudaMemcpyDeviceToHost);
  cudakernel::checkError("after copying the new positions");
  for(unsigned int i = 0; i < gpu_stat->host_parameter.num_atoms; ++i) {
    newpos[3*i  ] = double(gpu_stat->host_double_new_pos[i].x);
    newpos[3*i+1] = double(gpu_stat->host_double_new_pos[i].y);
    newpos[3*i+2] = double(gpu_stat->host_double_new_pos[i].z);
    DEBUG(15,"i : " << i << ". x : " << newpos[3*i  ])
    //DEBUG(15,"i : " << i << ". y : " << newpos[3*i+1])
    //DEBUG(15,"i : " << i << ". z : " << newpos[3*i+2])
  }

  return 0;

}

__global__ void cudakernel::kernel_CalcConstraints_Settle
(
        double3 * new_pos, double3 * old_pos,
        cudakernel::simulation_parameter * dev_params,
        int *shake_fail_mol
) {

  const int mol_index = blockIdx.x * NUM_THREADS_PER_BLOCK_SETTLE + threadIdx.x;

  if (mol_index >= dev_params->num_solvent_mol)
    return;
/*
  new_pos[mol_index * 3 + 0] = make_double3(mol_index, 0.0, 0.0);
  new_pos[mol_index * 3 + 1] = make_double3(mol_index, 0.0, 0.0);
  new_pos[mol_index * 3 + 2] = make_double3(mol_index, 0.0, 0.0);

*/
  // these parameters are for SPC water.
  const double mass_O = 15.99940;
  const double mass_H = 1.00800;
  const double dist_OH = 0.1000000;
  const double dist_HH = 0.1632990;

  // calculate the coordinates of the canonical triangle
  // see Figure 2 (a)
  const double half_mO_div_mH = 0.5 * mass_O / mass_H;
  const double rc = 0.5 * dist_HH;
  const double ra = sqrt(dist_OH * dist_OH - rc * rc) / (1.0 + half_mO_div_mH);
  const double rb = half_mO_div_mH * ra;

  // vectors in the plane of the old positions
  double3 b0 = old_pos[mol_index * 3 + 1] - old_pos[mol_index * 3 + 0];
  double3 c0 = old_pos[mol_index * 3 + 2] - old_pos[mol_index * 3 + 0];

  // centre of mass of new positions
  const double3 d0 = (new_pos[mol_index * 3 + 0] * mass_O +
                        new_pos[mol_index * 3 + 1] * mass_H + new_pos[mol_index * 3 + 2] * mass_H) /
                       (mass_O + mass_H + mass_H);

  // move the origin to the centre of mass
  double3 a1 = new_pos[mol_index * 3 + 0] - d0;
  double3 b1 = new_pos[mol_index * 3 + 1] - d0;
  double3 c1 = new_pos[mol_index * 3 + 2] - d0;

  // vectors describing transformation from original coordinate system to
  // the centre of mass originated coordinate system
  double3 n0 = cross(b0, c0);
  double3 n1 = cross(a1, n0);
  double3 n2 = cross(n0, n1);
  n0 = n0 / abs(n0); // this can give a NaN but it is very unlikely.
  n1 = n1 / abs(n1);
  n2 = n2 / abs(n2);

  // generate new normal vectors from the transpose in order to undo
  // the transformation
  const double3 m1 = make_double3(n1.x, n2.x, n0.x);
  const double3 m2 = make_double3(n1.y, n2.y, n0.y);
  const double3 m0 = make_double3(n1.z, n2.z, n0.z);

  // do the transformation to the centre of mass originated coordinate system
  // of the old positions
  b0 = make_double3(dot(n1, b0), dot(n2, b0), dot(n0, b0));
  c0 = make_double3(dot(n1, c0), dot(n2, c0), dot(n0, c0));

  // and of the new positions
  a1 = make_double3(dot(n1, a1), dot(n2, a1), dot(n0, a1));
  b1 = make_double3(dot(n1, b1), dot(n2, b1), dot(n0, b1));
  c1 = make_double3(dot(n1, c1), dot(n2, c1), dot(n0, c1));
  // now we can compute positions of canonical water
  // the cos is calculate from the square of the sin via sin^2 + cos^2 = 1
  const double sinphi = a1.z / ra; // this is (A8)
  const double one_minus_sinphi2 = 1.0 - sinphi*sinphi;
  if (one_minus_sinphi2 < 0.0) {
    *shake_fail_mol = mol_index;
    return;
  }
  const double cosphi = sqrt(one_minus_sinphi2);

  const double sinpsi = (b1.z - c1.z) / (2.0 * rc * cosphi); // this is (A9)
  const double one_minus_sinpsi2 = 1.0 - sinpsi*sinpsi;
  if (one_minus_sinpsi2 < 0.0) {
    *shake_fail_mol = mol_index;
    return;
  }
  const double cospsi = sqrt(one_minus_sinpsi2);

  // these are just to avoid recalculations
  const double minus_rb_cosphi = -rb * cosphi;
  const double rc_cospsi = rc * cospsi;
  const double rc_sinpsi_sinphi = rc * sinpsi*sinphi;
  const double rc_sinpsi_cosphi = rc * sinpsi*cosphi;

  // do the X. this is (A3)
  const double x_a2 = 0.0;
  const double x_b2 = - rc_cospsi;
  const double x_c2 = rc_cospsi;

  // do the Y. this is (A3)
  // I think there is an error in the paper. ra was confused with rc
  const double y_a2 = ra * cosphi;
  const double y_b2 = minus_rb_cosphi - rc_sinpsi_sinphi;
  const double y_c2 = minus_rb_cosphi + rc_sinpsi_sinphi;

  // do the Z components
  const double z_a2 = ra * sinphi; // this is (A5)
  const double z_b2 = -rb * sinphi + rc_sinpsi_cosphi; // this is (A6)
  const double z_c2 = -rb * sinphi - rc_sinpsi_cosphi; // this is (A7)

  // now construct the a2, b2 and c2 vectors
  const double3 a2 = make_double3(x_a2, y_a2, z_a2);
  const double3 b2 = make_double3(x_b2, y_b2, z_b2);
  const double3 c2 = make_double3(x_c2, y_c2, z_c2);

  // there are no a0 terms because we've already subtracted the term off
  // when we first defined b0 and c0.
  // this is (A15) but the equation is not really numbered...
  const double alpha = b2.x * (b0.x - c0.x) + b0.y * b2.y + c0.y * c2.y;
  const double beta = b2.x * (c0.y - b0.y) + b0.x * b2.y + c0.x * c2.y;
  const double gamma = b0.x * b1.y - b1.x * b0.y + c0.x * c1.y - c1.x * c0.y;

  const double alpha2_beta2 = alpha * alpha + beta * beta;
  // this is (A17)
  const double sintheta = (alpha * gamma -
          beta * sqrt(alpha2_beta2 - gamma * gamma)) / alpha2_beta2;
  const double one_minus_sintheta2 = 1.0 - sintheta * sintheta;
  if (one_minus_sintheta2 < 0.0) {
    *shake_fail_mol = mol_index;
    return;
  }
  const double costheta = sqrt(one_minus_sintheta2);

  // new finally a3, b3 and c3. These are (A4)
  const double3 a3 = make_double3(-a2.y * sintheta,
          a2.y * costheta,
          a1.z);
  const double3 b3 = make_double3(b2.x * costheta - b2.y * sintheta,
          b2.x * sintheta + b2.y * costheta,
          b1.z);
  const double3 c3 = make_double3(-b2.x * costheta - c2.y * sintheta,
          -b2.x * sintheta + c2.y * costheta,
          c1.z);

  // calculate the new positions
  new_pos[mol_index * 3 + 0] = make_double3(dot(a3, m1), dot(a3, m2), dot(a3, m0)) + d0;
  new_pos[mol_index * 3 + 1] = make_double3(dot(b3, m1), dot(b3, m2), dot(b3, m0)) + d0;
  new_pos[mol_index * 3 + 2] = make_double3(dot(c3, m1), dot(c3, m2), dot(c3, m0)) + d0;

}

#undef DEBUG

