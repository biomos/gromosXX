/**
 * @file cuda_kernel.cu
 * @author Poliak (peter.poliak@boku.ac.at)
 * @brief Implementation of the singleton controlling the GPUs
 * @version 0.1
 * @date 17.06.2023
 * 
 * @copyright Copyright (c) 2023
 * 
 */

#include "stdheader.h"
#include "cuda_kernel.h"
#include "util/debug.h"

#include "gpu_status.h"
//#include "cudaKernel.h"
#include "interaction.h"
#include "lib/constants.h"
#include "lib/utils.h"
#include "parameter.h"
#include "macros.h"

#include "algorithm/algorithm.h"
#include "../topology/topology.h"
#include "simulation/simulation.h"
#include "configuration/configuration.h"
#include "interaction/nonbonded/interaction/nonbonded_parameter.h"
#include "interaction/interaction_types.h"

#include "math/volume.h"

#include "pairlist/pairlist.h"

#undef MODULE
#undef SUBMODULE
#define MODULE cukernel
#define SUBMODULE kernel

#include "cukernel.cc"

__device__ __constant__ cukernel::simulation_parameter device_param;

// Initialize the CUDA_Kernel instance to nullptr
cukernel::CUDA_Kernel * cukernel::CUDA_Kernel::m_cuda_kernel = nullptr;

/*extern "C"*/ cukernel::CUDA_Kernel::CUDA_Kernel()
: mytopo(nullptr),
  myconf(nullptr),
  mysim(nullptr),
  device_number(-1)
{};


/*extern "C"*/ cukernel::CUDA_Kernel * cukernel::CUDA_Kernel::get_instance(
                                                const topology::Topology &topo,
                                                const configuration::Configuration &conf,
                                                const simulation::Simulation &sim) {
  if (m_cuda_kernel == nullptr) {
    m_cuda_kernel = new CUDA_Kernel();
    m_cuda_kernel->init(topo,conf,sim);
  } else {
    assert(&topo == m_cuda_kernel->mytopo);
    assert(&conf == m_cuda_kernel->myconf);
    assert(&sim == m_cuda_kernel->mysim);
    // also check, if the topo/conf/sim correspond to the parameters uploaded to the device
    assert(m_cuda_kernel->param == m_cuda_kernel->get_device_parameter());
  }
  return m_cuda_kernel;
};

/*extern "C"*/ void cukernel::CUDA_Kernel::update_nonbonded(interaction::Nonbonded_Parameter *np) {
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
  this->set_device_parameter();
}

/*extern "C"*/ void cukernel::CUDA_Kernel::set_device_parameter() const {
  cudaMemcpyToSymbol(device_param, &this->param, sizeof(cukernel::simulation_parameter));
}

cukernel::simulation_parameter cukernel::CUDA_Kernel::get_device_parameter() const {
  cukernel::simulation_parameter tmp_param;
  cudaMemcpyFromSymbol(&tmp_param, device_param, sizeof(cukernel::simulation_parameter));
  return tmp_param;
}

///*extern "C"*/ cukernel::CUDA_Kernel::~CUDA_Kernel(){};

/*extern "C"*/ void cukernel::CUDA_Kernel::init(
                        const topology::Topology &topo,
                        const configuration::Configuration &conf,
                        const simulation::Simulation &sim){
    this->mytopo = &topo;
    this->myconf = &conf;
    this->mysim = &sim;

    m_pos.resize(topo.num_atoms());
    m_force.resize(topo.num_atoms());

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
    cudaGetDeviceProperties(&this->device_prop, this->device_number);
    /*
    // let's first query the device properties and print them out
    cudaGetDeviceProperties(&this->device_prop, this->device_number);
    cudaGetDeviceFlags(&this->flags);
    if (this->flags == 0 && cudaSetDeviceFlags(cudaDeviceScheduleYield) == cudaErrorSetOnActiveProcess)
    std::cerr << "Cannot set flags\n";

    std::cout << "CUDA" << std::endl;
    std::cout << "\tDeviceproperties:" << std::endl;
    std::cout << "\t\tNumber: " << device_number << std::endl;
    std::cout << "\t\tName: " << this->device_prop.name << std::endl;
    std::cout << "\t\tTotal memory: " << this->device_prop.totalGlobalMem << std::endl;
    std::cout << "\t\tShared memory per block: " << this->device_prop.sharedMemPerBlock << std::endl;
    std::cout << "\t\tTotal constant memory: " << this->device_prop.totalConstMem << std::endl;

    // Create a gpu_status structure
    gpu_status * gpu_stat;
    gpu_stat = (gpu_status*) malloc (sizeof(gpu_status));
    gpu_stat->device = device_number;*/
};

/*extern "C"*/ void cukernel::CUDA_Kernel::copy_parameters() {
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
  this->set_device_parameter();
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

/*extern "C"*/ void cukernel::CUDA_Kernel::estimate_pairlist() {
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

int cukernel::CUDA_Kernel::copy_box(math::Box & box) {
  this->param.box.full.x = box(0)(0);
  this->param.box.full.y = box(1)(1);
  this->param.box.full.z = box(2)(2);
  this->param.box.inv.x = 1.0 / box(0)(0);
  this->param.box.inv.y = 1.0 / box(1)(1);
  this->param.box.inv.z = 1.0 / box(2)(2);
  this->param.box.half.x = box(0)(0) / 2.0;
  this->param.box.half.y = box(1)(1) / 2.0;
  this->param.box.half.z = box(2)(2) / 2.0;
  // copy only the box
  cudaMemcpyToSymbol(device_param, &this->param.box, sizeof(cukernel::simulation_parameter::box_struct), offsetof(cukernel::simulation_parameter, box));
  return cukernel::check_error("after copying the box");
}

int cukernel::CUDA_Kernel::copy_positions(math::VArray & pos) {
  // copy only the box
  assert(pos.size() == this->m_pos.size());
  typedef decltype(this->m_pos)::value_type vec_type; // float3
  typedef decltype(vec_type::x) val_type; // float
  for (unsigned i = 0; i < pos.size(); ++i) {
    this->m_pos[i].x = static_cast<val_type>(pos(i)(0));
    this->m_pos[i].y = static_cast<val_type>(pos(i)(1));
    this->m_pos[i].z = static_cast<val_type>(pos(i)(2));
  }
  cudaMemPrefetchAsync(this->m_pos.data(), this->m_pos.size() * sizeof(vec_type), this->device_number);
  return cukernel::check_error("after copying the positions");
}

void cukernel::CUDA_Kernel::update_pairlist(topology::Topology &topo,
                                            configuration::Configuration &conf,
                                            simulation::Simulation &sim) {
  /** number of blocks should be multiples of
   * the maximum number of blocks per multiprocessor
   */
  const unsigned threads_per_block = 128;
   /** threads_per_block should be optimized for the specific
    * problem size and compute capability (occupancy calculator)
    */
  //const unsigned max_blocks_per_mp = this->device_prop.maxBlocksPerMultiProcessor;
  const unsigned max_blocks_per_mp = 2048;
  const unsigned max_threads_per_mp = this->device_prop.maxThreadsPerMultiProcessor;
  assert(max_blocks_per_mp * max_threads_per_mp % threads_per_block == 0 && "Set better threads_per_block");
  const unsigned num_blocks = max_blocks_per_mp * max_threads_per_mp / threads_per_block;
  dim3 grid(num_blocks);
  dim3 block(threads_per_block);

  DEBUG(10,"Pairlist: GPU ID: " << this->device_number << " of " << sim.param().cuda.number_gpus
            <<  ". Blocks: " << num_blocks);
  

  Pairlist pairlist(topo.num_atoms(), 2000);

  bool overflow = false;
  do {
      //pairlist.clear();
      cukernel::update_pairlist<<<grid,block>>>(m_pos.data(), pairlist);
    
      //std::cout << "cell_i: " << cell_i << std::endl;
      overflow = pairlist.overflown();
      
      if (overflow) {
          unsigned pair_size = 2*pairlist.width();
          std::cout << "Pairlist overflown, new parlist width: " << pair_size*sizeof(Pairlist::num_type) << std::endl;
          pairlist.resize(topo.num_atoms(), pair_size);
      }
  } while (overflow);
};

void cukernel::CUDA_Kernel::calculate_interactions() {
  cukernel::calculate_interactions<<<32,32>>>();
};

// __global__ void cukernel::update_pairlist(float3 * dev_pos,
//                                           pairlist pl_short,
//                                           pairlist pl_long,
//                                           bool overflow) {
//   // two tables
//   // to convert from global atom index to local index
//   __shared__ char shmem[];
//   __shared__ uint2 glob_to_loc[];
//   __shared__ uint loc_to_glob[];
//   // local thread index (within the block)
//   const unsigned lid = threadIdx.x;
//   // global thread index
//   const unsigned tid = threadIdx.x + blockDim.x * blockIdx.x;
//   // loop over chargegroups or atoms


//   // search for all the indices in our rows
//   unsigned int num_neighbors_long = 0, num_neighbors_short = 0;
//   __shared__ float shared_pos[NUM_THREADS_PER_BLOCK * 3];

//   // take host_parameter local
//   const unsigned int &N = device_param.num_atoms.solvent;
//   const unsigned int &num_solvent_mol = device_param.num_solvent_mol;
//   const float &cutoff_long_2 = device_param.cutoff_long_2;
//   const float &cutoff_short_2 = device_param.cutoff_short_2;
//   //box edges
//   const cukernel::simulation_parameter::box_struct& device_box = device_param.box;
//   const float &box_x = device_box.full.x;
//   const float &box_y = device_box.full.y;
//   const float &box_z = device_box.full.z;

//   const float &box_inv_x = device_box.inv.x;
//   const float &box_inv_y = device_box.inv.y;
//   const float &box_inv_z = device_box.inv.z;
//   const unsigned int solvent_offset = device_param.num_atoms_per_mol;

//   // calculate indices
//   const unsigned int index = blockIdx.x * NUM_THREADS_PER_BLOCK + threadIdx.x;
//   const unsigned int molecule_index = index * num_of_gpus + gpu_id;

//   const unsigned int first_atom_index = molecule_index*solvent_offset;
//   const unsigned int myThreadOffset = threadIdx.x*solvent_offset;

//   float3 first_atom_pos;
//   if (first_atom_index < N)
//     first_atom_pos = dev_pos[first_atom_index];

//   for (unsigned int i = 0; i < N; i += (NUM_THREADS_PER_BLOCK * solvent_offset)) {
//     float3 neighbor_pos;
//     if (i + myThreadOffset < N)
//       neighbor_pos = dev_pos[i + myThreadOffset];

//     // cache a block of positions
//     __syncthreads();
//     shared_pos[threadIdx.x] = neighbor_pos.x;
//     shared_pos[threadIdx.x + NUM_THREADS_PER_BLOCK] = neighbor_pos.y;
//     shared_pos[threadIdx.x + 2 * NUM_THREADS_PER_BLOCK] = neighbor_pos.z;
//     __syncthreads();

//     /*unsigned int end_i_loop = NUM_THREADS_PER_BLOCK;
//     if (end_i_loop > (N - i) / solvent_offset)
//       end_i_loop = (N - i) / solvent_offset;*/
    
//     const unsigned int end_i_loop = min(NUM_THREADS_PER_BLOCK, (N - i) / solvent_offset);
//     if (first_atom_index < N) {
//       for (unsigned int start_i_loop = 0; start_i_loop < end_i_loop; start_i_loop++) {
//         const unsigned int current_first_atom_index = i + start_i_loop*solvent_offset;
//         if (current_first_atom_index != first_atom_index && current_first_atom_index < N) {
//           //{ calculate distance and d^2
//           float3 dist;
//           dist.x = first_atom_pos.x - shared_pos[start_i_loop];
//           dist.x -= box_x * rintf(dist.x * box_inv_x);
//           dist.y = first_atom_pos.y - shared_pos[start_i_loop + NUM_THREADS_PER_BLOCK];
//           dist.y -= box_y * rintf(dist.y * box_inv_y);
//           dist.z = first_atom_pos.z - shared_pos[start_i_loop + 2 * NUM_THREADS_PER_BLOCK];
//           dist.z -= box_z * rintf(dist.z * box_inv_z);
//           const float d2 = abs2(dist);
//           //} calculate distance and d^2
//        // are they interacting?
//        if (d2 < cutoff_long_2) {
//             // longrange?
//             if (d2 > cutoff_short_2) {
//               if (num_neighbors_long < pl_long.max_size) {
//                 pl_long.list[index + pl_long.pitch * num_neighbors_long] = current_first_atom_index;
//                 num_neighbors_long++;
//               } else {
//                 *pl_long.overflow = true;
//               } // overflow
//             } else { // shortrange then
//               if (num_neighbors_short < pl_short.max_size) {
//                 pl_short.list[index + pl_short.pitch * num_neighbors_short] = current_first_atom_index;
//                 num_neighbors_short++;
//               } else {
//                 *pl_short.overflow = true;
//               } // overflow
//             } // if shortrange / longrange
//           } // if cutoff
//         } // if atom in valid range
//       } // for atoms j
//     } // if atom in valid range
//   } // for atoms i
//   if (molecule_index < num_solvent_mol) {
//     pl_long.num_neighbors[index] = num_neighbors_long;
//     pl_short.num_neighbors[index] = num_neighbors_short;
//   }
// }