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
#include "parameter.h"
#include "macros.h"

#include "algorithm/algorithm.h"
#include "../topology/topology.h"
#include "simulation/simulation.h"
#include "configuration/configuration.h"
#include "interaction/nonbonded/interaction/nonbonded_parameter.h"
#include "interaction/interaction_types.h"

#include "math/volume.h"


#undef MODULE
#undef SUBMODULE
#define MODULE cukernel
#define SUBMODULE kernel

#include "cuda.cc"

__device__ __constant__ cukernel::simulation_parameter device_param;

// Initialize the CUDA_Kernel instance to nullptr
cukernel::CUDA_Kernel * cukernel::CUDA_Kernel::m_cuda_kernel = nullptr;

extern "C" cukernel::CUDA_Kernel::CUDA_Kernel()
: mytopo(nullptr),
  myconf(nullptr),
  mysim(nullptr),
  device_number(-1)
{};


extern "C" cukernel::CUDA_Kernel * cukernel::CUDA_Kernel::get_instance(
                                                topology::Topology &topo,
                                                configuration::Configuration &conf,
                                                simulation::Simulation &sim) {
  if (m_cuda_kernel == nullptr) {
    m_cuda_kernel = new CUDA_Kernel();
    m_cuda_kernel->init(topo,conf,sim);
  } else {
    assert(&topo == m_cuda_kernel->mytopo);
    assert(&conf == m_cuda_kernel->myconf);
    assert(&sim == m_cuda_kernel->mysim);
  }
  return m_cuda_kernel;
};

extern "C" void cukernel::CUDA_Kernel::update_nonbonded(interaction::Nonbonded_Parameter *np) {
  // solvent lj and crf parameters to const memory
  if (this->param.num_atoms_per_mol > simulation_parameter::max_atoms_solvent) {
    io::messages.add("Too many solvent atoms for constant memory",
                        "CUDA_Kernel", io::message::critical);
  }
  memset(&this->param.solvent_lj_crf, 0, sizeof(this->param.solvent_lj_crf));

  const unsigned solvent_index = mytopo->num_solute_atoms();
  for (unsigned i = 0; i < this->param.num_atoms_per_mol; ++i) {
    for (unsigned j = 0; j < this->param.num_atoms_per_mol; ++j) {
      lj_crf_parameter & lj_crf_pair = this->param.solvent_lj_crf[i*this->param.num_atoms_per_mol + j];
      const interaction::lj_parameter_struct & lj = np->lj_parameter(mytopo->iac(solvent_index + i), mytopo->iac(solvent_index + j));
      lj_crf_pair.q = math::four_pi_eps_i * mytopo->charge(solvent_index + i) * mytopo->charge(solvent_index + j);
      lj_crf_pair.c6 = lj.c6;
      lj_crf_pair.c12 = lj.c12;
    }
  }
  this->sync_symbol();
}

extern "C" void cukernel::CUDA_Kernel::sync_symbol() {
  cudaMemcpyToSymbol(device_param, &this->param, sizeof(cukernel::simulation_parameter));
}

//extern "C" cukernel::CUDA_Kernel::~CUDA_Kernel(){};

extern "C" void cukernel::CUDA_Kernel::init(topology::Topology &topo,
                        configuration::Configuration &conf,
                        simulation::Simulation &sim){
    this->mytopo = &topo;
    this->myconf = &conf;
    this->mysim = &sim;

    // initialize the devices
    // for now, we only support only single GPU
    DEBUG(0, "Device number:" << this->device_number);
    
    /*
     TEMPORARY OVERRIDE
    */
    {
      this->device_number = -1;
      if (mysim->param().cuda.device_number.size()) {
        this->device_number = mysim->param().cuda.device_number.at(0);
      }
    }



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

extern "C" void cukernel::CUDA_Kernel::copy_parameters() {
  DEBUG(0, "cutoff before copy: " << mysim->param().pairlist.cutoff_long);
  // long cutoff
  this->param.cutoff_long = mysim->param().pairlist.cutoff_long;
  this->param.cutoff_long_2 = (FP_CAST) mysim->param().pairlist.cutoff_long * mysim->param().pairlist.cutoff_long;
  
  #ifndef NDEBUG
    float sp_cutoff_long_2 = float(mysim->param().pairlist.cutoff_long) * float(mysim->param().pairlist.cutoff_long);
    double dp_cutoff_long_2 = mysim->param().pairlist.cutoff_long * mysim->param().pairlist.cutoff_long;
    float dpsp_cutoff_long_2 = dp_cutoff_long_2;
    float dpspsp_cutoff_long_2 = (float) mysim->param().pairlist.cutoff_long * mysim->param().pairlist.cutoff_long;
    DEBUG(0, "double cutoff_long_2: " << std::scientific << std::setprecision(18) << dp_cutoff_long_2);
    DEBUG(0, "float cutoff_long_2: " << std::scientific << std::setprecision(18) << sp_cutoff_long_2);
    double truncation_error = 100.*((double)sp_cutoff_long_2 - dp_cutoff_long_2) / dp_cutoff_long_2;
    DEBUG(0, "truncation error: " << std::scientific << std::setprecision(18) << truncation_error << " %");

    DEBUG(0, "double->float cutoff_long_2: " << std::scientific << std::setprecision(18) << dpsp_cutoff_long_2);
    truncation_error = 100.*((double)dpsp_cutoff_long_2 - dp_cutoff_long_2) / dp_cutoff_long_2;
    DEBUG(0, "truncation error: " << std::scientific << std::setprecision(18) << truncation_error << " %");

    DEBUG(0, "float*double->float cutoff_long_2: " << std::scientific << std::setprecision(18) << dpspsp_cutoff_long_2);
    truncation_error = 100.*((double)dpspsp_cutoff_long_2 - dp_cutoff_long_2) / dp_cutoff_long_2;
    DEBUG(0, "truncation error: " << std::scientific << std::setprecision(18) << truncation_error << " %");
  #endif

  // short cutoff
  this->param.cutoff_short = mysim->param().pairlist.cutoff_short;
  this->param.cutoff_short_2 = (FP_CAST) mysim->param().pairlist.cutoff_short * mysim->param().pairlist.cutoff_short;
  
  #ifndef NDEBUG
    float sp_cutoff_short_2 = float(mysim->param().pairlist.cutoff_short) * float(mysim->param().pairlist.cutoff_short);
    double dp_cutoff_short_2 = mysim->param().pairlist.cutoff_short * mysim->param().pairlist.cutoff_short;
    float dpsp_cutoff_short_2 = dp_cutoff_short_2;
    float dpspsp_cutoff_short_2 = (float) mysim->param().pairlist.cutoff_short *  mysim->param().pairlist.cutoff_short;
    DEBUG(0, "double cutoff_short_2: " << std::scientific << std::setprecision(18) << dp_cutoff_short_2);

    DEBUG(0, "float cutoff_short_2: " << std::scientific << std::setprecision(18) << sp_cutoff_short_2);
    truncation_error = 100.*((double)sp_cutoff_short_2 - dp_cutoff_short_2) / dp_cutoff_short_2;
    DEBUG(0, "truncation error: " << std::scientific << std::setprecision(18) << truncation_error << " %");

    DEBUG(0, "double->float cutoff_short_2: " << std::scientific << std::setprecision(18) << dpsp_cutoff_short_2);
    truncation_error = 100.*((double)dpsp_cutoff_short_2 - dp_cutoff_short_2) / dp_cutoff_short_2;
    DEBUG(0, "truncation error: " << std::scientific << std::setprecision(18) << truncation_error << " %");
    
    DEBUG(0, "float*double->float cutoff_short_2: " << std::scientific << std::setprecision(18) << dpspsp_cutoff_short_2);
    truncation_error = 100.*((double)dpspsp_cutoff_short_2 - dp_cutoff_short_2) / dp_cutoff_short_2;
    DEBUG(0, "truncation error: " << std::scientific << std::setprecision(18) << truncation_error << " %");
    
  #endif
  // box edges
  this->param.box.full.x = myconf->current().box(0)(0);
  this->param.box.full.y = myconf->current().box(1)(1);
  this->param.box.full.z = myconf->current().box(2)(2);

  // inverted box edges
  this->param.box.inv.x = (FP_CAST) 1. / myconf->current().box(0)(0);
  this->param.box.inv.y = (FP_CAST) 1. / myconf->current().box(1)(1);
  this->param.box.inv.z = (FP_CAST) 1. / myconf->current().box(2)(2);

  // half the box edges
  this->param.box.half.x = (FP_CAST) myconf->current().box(0)(0) / 2.;
  this->param.box.half.y = (FP_CAST) myconf->current().box(1)(1) / 2.;
  this->param.box.half.z = (FP_CAST) myconf->current().box(2)(2) / 2.;
  
  // reaction field constants
  FP_CAST m_cut3i, m_crf, m_crf_cut, m_crf_cut3i, m_crf_2cut3i;
  FP_CAST rf_cutoff = mysim->param().nonbonded.rf_cutoff;
  FP_CAST epsilon = mysim->param().nonbonded.epsilon;
  FP_CAST rf_epsilon = mysim->param().nonbonded.rf_epsilon;
  FP_CAST rf_kappa = mysim->param().nonbonded.rf_kappa;

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
  
  // number of atoms
  this->param.num_atoms.total = mytopo->num_atoms();

  // number of solute atoms
  this->param.num_atoms.solute = mytopo->num_solute_atoms();

  // number of solvent atoms
  assert(mytopo->num_solvents() <= 1);
  this->param.num_atoms.solvent = mytopo->num_solvent_atoms(0);

  // the number of solvent molecules
  this->param.num_solvent_mol = mytopo->num_solvent_molecules(0);
  // the number of atoms per solvent molecule
  this->param.num_atoms_per_mol = mytopo->num_solvent_atoms(0) / mytopo->num_solvent_molecules(0);
  // the number of gpus
  this->param.num_of_gpus = mysim->param().cuda.number_gpus;
  //cudaMemcpyToSymbol(device_param.num_of_gpus, &m_num_gpus, sizeof(unsigned));
  this->param.gpu_id = this->device_number;
  
  this->estimate_pairlist();
  this->sync_symbol();
  // check what we have there
  cukernel::simulation_parameter tmp_param;
  cudaMemcpyFromSymbol(&tmp_param, device_param, sizeof(cukernel::simulation_parameter));
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
  return;
};

extern "C" void cukernel::CUDA_Kernel::estimate_pairlist() {
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

