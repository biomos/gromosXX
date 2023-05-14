/**
 * @file cudaKernel.cu
 * implementation of CUDA kernel
 */

#include "stdheader.h"
#include "cuda_kernel.h"
#include "util/debug.h"

#include "gpu_status.h"
//#include "cudaKernel.h"
#include "lib/constants.h"
#include "cukernel/parameter.h"

#include "algorithm/algorithm.h"
#include "../topology/topology.h"
#include "simulation/simulation.h"
#include "configuration/configuration.h"

#include "math/volume.h"

#undef MODULE
#undef SUBMODULE
#define MODULE cuda
#define SUBMODULE kernel

#include "cuda.cc"

__device__ __constant__ cudakernel::simulation_parameter device_param;

extern "C" cudakernel::CUDA_Kernel::CUDA_Kernel()
: mytopo(nullptr),
  myconf(nullptr),
  mysim(nullptr),
  device_number(-1)
{};

extern "C" cudakernel::CUDA_Kernel::~CUDA_Kernel(){};

extern "C" void cudakernel::CUDA_Kernel::init(topology::Topology &topo,
                        configuration::Configuration &conf,
                        simulation::Simulation &sim){
    this->mytopo = &topo;
    this->myconf = &conf;
    this->mysim = &sim;

    // initialize the devices
    // for now, we only support only single GPU
    DEBUG(0, "Device number:" << this->device_number);
    this->device_number = mysim->param().cuda.device_number.at(0);
    DEBUG(0, "Device number:" << this->device_number);
    if (this->device_number != -1) {
      cudaSetDevice(this->device_number);
    } else {
      // determine the device number
      cudaGetDevice(&this->device_number);
    }
    DEBUG(0, "Device number:" << this->device_number);
    
    float cutoff_short = mysim->param().pairlist.cutoff_short;
    //cudaMemcpyToSymbol(device::cutoff.cshort, &cutoff_short, sizeof(float));
    this->copy_parameters();
    DEBUG(0, "CUDA copied to symbol");
    /*
    // let's first query the device properties and print them out
    cudaGetDeviceProperties(&this->device_properties, this->device_number);
    cudaGetDeviceFlags(&this->flags);
    if (this->flags == 0 && cudaSetDeviceFlags(cudaDeviceScheduleYield) == cudaErrorSetOnActiveProcess)
    std::cerr << "Cannot set flags\n";

    std::cout << "CUDA" << std::endl;
    std::cout << "\tDeviceproperties:" << std::endl;
    std::cout << "\t\tNumber: " << device_number << std::endl;
    std::cout << "\t\tName: " << this->device_properties.name << std::endl;
    std::cout << "\t\tTotal memory: " << this->device_properties.totalGlobalMem << std::endl;
    std::cout << "\t\tShared memory per block: " << this->device_properties.sharedMemPerBlock << std::endl;
    std::cout << "\t\tTotal constant memory: " << this->device_properties.totalConstMem << std::endl;

    // Create a gpu_status structure
    gpu_status * gpu_stat;
    gpu_stat = (gpu_status*) malloc (sizeof(gpu_status));
    gpu_stat->device = device_number;*/
};

extern "C" void cudakernel::CUDA_Kernel::copy_parameters() {
  DEBUG(0, "cutoff before copy: " << mysim->param().pairlist.cutoff_long);
  float m_cutoff;
  // long cutoff
  m_cutoff = mysim->param().pairlist.cutoff_long;
  this->param.cutoff_long = mysim->param().pairlist.cutoff_long;
  //cudaMemcpyToSymbol(device_param.cutoff_long, &m_cutoff, sizeof(float));
  m_cutoff = m_cutoff * m_cutoff;
  this->param.cutoff_long_2 = mysim->param().pairlist.cutoff_long * mysim->param().pairlist.cutoff_long;
  //cudaMemcpyToSymbol(device_param.cutoff_long_2, &m_cutoff, sizeof(float));

  // short cutoff
  m_cutoff = mysim->param().pairlist.cutoff_short;
  this->param.cutoff_short = mysim->param().pairlist.cutoff_short;
  //cudaMemcpyToSymbol(device_param.cutoff_short, &m_cutoff, sizeof(float));
  m_cutoff = m_cutoff * m_cutoff;
  this->param.cutoff_short_2 = mysim->param().pairlist.cutoff_short * mysim->param().pairlist.cutoff_short;
  //cudaMemcpyToSymbol(device_param.cutoff_short_2, &m_cutoff, sizeof(float));

  // box edges
  float3 m_box;
  m_box.x = myconf->current().box(0)(0);
  m_box.y = myconf->current().box(1)(1);
  m_box.z = myconf->current().box(2)(2);
  this->param.box.x = myconf->current().box(0)(0);
  this->param.box.y = myconf->current().box(1)(1);
  this->param.box.z = myconf->current().box(2)(2);
  //cudaMemcpyToSymbol(device_param.box, &m_box, sizeof(float3));

  // inverted box edges
  m_box.x = 1 / m_box.x;
  m_box.y = 1 / m_box.y;
  m_box.z = 1 / m_box.z;
  this->param.box_inv.x = 1 / myconf->current().box(0)(0);
  this->param.box_inv.y = 1 / myconf->current().box(1)(1);
  this->param.box_inv.z = 1 / myconf->current().box(2)(2);
  //cudaMemcpyToSymbol(device_param.box_inv, &m_box, sizeof(float3));

  // half the box edges
  m_box.x = myconf->current().box(0)(0) / 2;
  m_box.y = myconf->current().box(1)(1) / 2;
  m_box.z = myconf->current().box(2)(2) / 2;
  this->param.box_half.x = myconf->current().box(0)(0) / 2;
  this->param.box_half.y = myconf->current().box(1)(1) / 2;
  this->param.box_half.z = myconf->current().box(2)(2) / 2;
  //cudaMemcpyToSymbol(device_param.box_half, &m_box, sizeof(float3));
  
  // reaction field constants
  float m_cut3i, m_crf, m_crf_cut, m_crf_cut3i, m_crf_2cut3i;
  float rf_cutoff = mysim->param().nonbonded.rf_cutoff;
  float epsilon = mysim->param().nonbonded.epsilon;
  float rf_epsilon = mysim->param().nonbonded.rf_epsilon;
  float rf_kappa = mysim->param().nonbonded.rf_kappa;

  m_cut3i = 1.0 / (rf_cutoff * rf_cutoff * rf_cutoff);
  m_crf = 2 * (epsilon - rf_epsilon)*(1.0 + rf_kappa * rf_cutoff) - rf_epsilon *
          (rf_kappa * rf_cutoff * rf_kappa * rf_cutoff);
  m_crf /= (epsilon + 2 * rf_epsilon)*(1.0 + rf_kappa * rf_cutoff) +
          rf_epsilon * (rf_kappa * rf_cutoff * rf_kappa * rf_cutoff);
  m_crf_cut3i = m_crf*m_cut3i;

  m_crf_2cut3i = m_crf_cut3i / 2.0;
  m_crf_cut = (1 - m_crf / 2.0) / rf_cutoff;

  this->param.crf_2cut3i = m_crf_2cut3i;
  this->param.crf_cut = m_crf_cut;
  this->param.crf_cut3i = m_crf_cut3i;
  //cudaMemcpyToSymbol(device_param.crf_2cut3i, &m_crf_2cut3i, sizeof(float));
  //cudaMemcpyToSymbol(device_param.crf_cut, &m_crf_cut, sizeof(float));
  //cudaMemcpyToSymbol(device_param.crf_cut3i, &m_crf_cut3i, sizeof(float));
  
  // number of atoms
  unsigned m_num_atoms = mytopo->num_atoms();
  this->param.num_atoms.total = mytopo->num_atoms();
  //cudaMemcpyToSymbol(device_param.num_atoms.total, &m_num_atoms, sizeof(unsigned));

  // number of solute atoms
  m_num_atoms = mytopo->num_solute_atoms();
  this->param.num_atoms.solute = mytopo->num_solute_atoms();
  //cudaMemcpyToSymbol(device_param.num_atoms.solute, &m_num_atoms, sizeof(unsigned));

  // number of solvent atoms
  m_num_atoms = mytopo->num_solvent_atoms();
  this->param.num_atoms.solvent = mytopo->num_solvent_atoms();
  //cudaMemcpyToSymbol(device_param.num_atoms.solvent, &m_num_atoms, sizeof(unsigned));

  // the number of atoms per solvent molecule
  assert(mytopo->num_solvents() <= 1);
  unsigned m_num_solvent_mol = mytopo->num_solvent_molecules(0);
  this->param.num_atoms.solvent = mytopo->num_solvent_atoms(0);
  this->param.num_solvent_mol = mytopo->num_solvent_molecules(0);
  //cudaMemcpyToSymbol(device_param.num_solvent_mol, &m_num_solvent_mol, sizeof(unsigned));
  m_num_atoms = mytopo->num_solvent_atoms() / m_num_solvent_mol;
  this->param.num_atoms_per_mol = mytopo->num_solvent_atoms() / mytopo->num_solvent_molecules(0);
  //cudaMemcpyToSymbol(device_param.num_atoms_per_mol, &m_num_atoms, sizeof(unsigned));
  unsigned m_num_gpus = mysim->param().cuda.number_gpus;
  this->param.num_of_gpus = mysim->param().cuda.number_gpus;
  //cudaMemcpyToSymbol(device_param.num_of_gpus, &m_num_gpus, sizeof(unsigned));
  
  this->estimate_pairlist();
  cudaMemcpyToSymbol(device_param, &this->param, sizeof(cudakernel::simulation_parameter));
  // check what we have there
  cudakernel::simulation_parameter tmp_param;
  cudaMemcpyFromSymbol(&tmp_param, device_param, sizeof(cudakernel::simulation_parameter));
  DEBUG(0, "num_atoms.total: \t" << tmp_param.num_atoms.total);
  DEBUG(0, "num_atoms.solute: \t" << tmp_param.num_atoms.solute);
  DEBUG(0, "num_atoms.solvent: \t" << tmp_param.num_atoms.solvent);
  DEBUG(0, "box: \t" << tmp_param.box.x << " " << tmp_param.box.y << " " << tmp_param.box.z);
  DEBUG(0, "box_inv: \t" << tmp_param.box_inv.x << " " << tmp_param.box_inv.y << " " << tmp_param.box_inv.z);
  DEBUG(0, "box_half: \t" << tmp_param.box_half.x << " " << tmp_param.box_half.y << " " << tmp_param.box_half.z);
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
  return;
};

extern "C" void cudakernel::CUDA_Kernel::estimate_pairlist() {
  // calculate density and allocate memory?
  // when we overflow, double the size and start over again?
  // calculate the particle density?
  double particle_density = 0.;
  double volume = math::volume(myconf->current().box, myconf->boundary_type);
  unsigned m_estimated_neighbors_long = 500;
  unsigned m_estimated_neighbors_short = 300;
  if (volume) {
    particle_density = mytopo->num_atoms() / volume;
    DEBUG(0,"Volume current: " << math::volume(myconf->current().box, myconf->boundary_type));
    DEBUG(0,"particle_density: " << particle_density);
    double m_cutoff_long3 = mysim->param().pairlist.cutoff_long;
    m_cutoff_long3 *= m_cutoff_long3 * m_cutoff_long3;
    double m_cutoff_short3 = mysim->param().pairlist.cutoff_short;
    m_cutoff_short3 *= m_cutoff_short3 * m_cutoff_short3;
    DEBUG(0,"m_cutoff_short3: " << m_cutoff_short3);
    DEBUG(0,"m_cutoff_long3: " << m_cutoff_long3);
    double sphere_volume_short = 4. / 3. * math::Pi * m_cutoff_short3;
    double sphere_volume_long = 4. / 3. * math::Pi * m_cutoff_long3 - sphere_volume_short;
    m_estimated_neighbors_long = 1.2 * sphere_volume_long * particle_density;
    m_estimated_neighbors_short = 1.2 * sphere_volume_short * particle_density;
  }
  this->param.estimated_neighbors_long = m_estimated_neighbors_long;
  this->param.estimated_neighbors_short = m_estimated_neighbors_short;
  DEBUG(0,"long pairlist size: " << this->param.estimated_neighbors_long);
  DEBUG(0,"short pairlist size: " << this->param.estimated_neighbors_short);
  //cudaMemcpyToSymbol(device_param.estimated_neighbors_long, &m_estimated_neighbors_long, sizeof(unsigned));
  //cudaMemcpyToSymbol(device_param.estimated_neighbors_short, &m_estimated_neighbors_short, sizeof(unsigned));
  return;
};

