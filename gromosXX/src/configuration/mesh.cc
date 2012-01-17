/**
 * @file mesh.cc
 * Mesh class
 */

#include "../stdheader.h"
#include "../configuration/mesh.h"
#include "../util/error.h"
#include "../util/debug.h"

// for memcpy
#include <cstring>

#undef MODULE
#undef SUBMODULE
#define MODULE configuration
#define SUBMODULE configuration

template<typename complex_type>
configuration::GenericMesh<complex_type>::GenericMesh() : m_x(0), m_y(0), m_z(0), m_volume(0), 
        m_mesh(NULL), plan_forward(NULL), plan_backward(NULL) {
}

// this function is only used for complex numbers which has to be aligned
// for FFTW
namespace configuration {
template<>
GenericMesh<complex_number>::GenericMesh(unsigned int x, unsigned int y, unsigned int z) :
m_x(x), m_y(y), m_z(z), m_volume(x*y*z) {
  // allocate the grid arrays. Here we have to use a dynamic array
  // because static arrays would lie on the stack and are not aligned
  // in the SIMD way. 
  // Elements can be accessed by:
  // statarr[x][y][z] == dynarr[z + Nz * (y + Ny * x)]
  // for more details see 3.2 Multidimensional Array format in section
  
  m_mesh = (complex_number*) FFTW3(malloc(m_volume * sizeof(FFTW3(complex))));
  plan_forward = FFTW3(plan_dft_3d(m_x, m_y, m_z,
              reinterpret_cast<FFTW3(complex)*> (m_mesh), reinterpret_cast<FFTW3(complex)*> (m_mesh),
              FFTW_FORWARD, FFTW_ESTIMATE));
  plan_backward = FFTW3(plan_dft_3d(m_x, m_y, m_z,
              reinterpret_cast<FFTW3(complex)*> (m_mesh), reinterpret_cast<FFTW3(complex)*> (m_mesh),
              FFTW_BACKWARD, FFTW_ESTIMATE));
}

// this function is used for all other types.
template<typename complex_type>
GenericMesh<complex_type>::GenericMesh(unsigned int x, unsigned int y, unsigned int z) :
m_x(x), m_y(y), m_z(z), m_volume(x*y*z) {
  m_mesh = new complex_type[m_volume];
}

// for complex number we also have to destroy the plans
template<>
GenericMesh<complex_number>::~GenericMesh() {
  // destory the plans
  if (plan_forward != NULL)
    FFTW3(destroy_plan(plan_forward));
  if (plan_backward != NULL)
    FFTW3(destroy_plan(plan_backward));
  
  if (m_mesh != NULL) {
    // make sure you free the mesh. It can be big.
    FFTW3(free((configuration::complex_number*) m_mesh));
  }
}

// this function is used for all other types.
template<typename complex_type>
GenericMesh<complex_type>::~GenericMesh() {
  if (m_mesh != NULL)
    delete m_mesh;
}

template<>
GenericMesh<complex_number>::GenericMesh(const GenericMesh<complex_number> & mesh) {
  m_x = mesh.x();
  m_y = mesh.y();
  m_z = mesh.z();
  m_volume = m_x * m_y * m_z;
  // reserve memory for the new mesh
  m_mesh = (configuration::complex_number*) FFTW3(malloc(m_volume * sizeof(configuration::complex_number)));
  // just copy the memory
  memcpy(m_mesh, mesh.m_mesh, m_volume * sizeof(configuration::complex_number));
}

template<typename complex_type>
GenericMesh<complex_type>::GenericMesh(const GenericMesh<complex_type> & mesh) {
  m_x = mesh.x();
  m_y = mesh.y();
  m_z = mesh.z();
  m_volume = m_x * m_y * m_z;
  // reserve memory for the new mesh
  m_mesh = new complex_type[10];
  const complex_type* old_mesh = mesh.mesh();
  for(unsigned int i = 0; i < m_volume; ++i)
    m_mesh[i] = old_mesh[i]; // use the assignment operator
}

// for the complex number we have to align
template<>
void GenericMesh<complex_number>::resize(unsigned int x, unsigned int y, unsigned int z){
  m_x = x;
  m_y = y;
  m_z = z;
  m_volume = m_x * m_y * m_z;
  
  if (plan_forward != NULL)
    FFTW3(destroy_plan(plan_forward));
  if (plan_backward != NULL)
    FFTW3(destroy_plan(plan_forward));
  
  if (m_mesh != NULL)
    FFTW3(free(m_mesh));
  
  m_mesh = (complex_number*) FFTW3(malloc(m_volume * sizeof(FFTW3(complex))));
  plan_forward = FFTW3(plan_dft_3d(m_x, m_y, m_z,
              reinterpret_cast<FFTW3(complex)*> (m_mesh), reinterpret_cast<FFTW3(complex)*> (m_mesh),
              FFTW_FORWARD, FFTW_ESTIMATE));
  plan_backward = FFTW3(plan_dft_3d(m_x, m_y, m_z,
              reinterpret_cast<FFTW3(complex)*> (m_mesh), reinterpret_cast<FFTW3(complex)*> (m_mesh),
              FFTW_BACKWARD, FFTW_ESTIMATE));
}

template<typename complex_type>
void GenericMesh<complex_type>::resize(unsigned int x, unsigned int y, unsigned int z){
  m_x = x;
  m_y = y;
  m_z = z;
  m_volume = m_x * m_y * m_z;
  
  if (m_mesh != NULL)
    delete m_mesh;
  
  m_mesh = new complex_type[m_volume];
}

// this function is only for the complex number
template<>
void GenericMesh<complex_number>::fft(GenericMesh<complex_number>::fft_type type) {
  switch (type) {
    case GenericMesh<complex_number>::fft_forward :
              FFTW3(execute(plan_forward));
      break;
    case GenericMesh<complex_number>::fft_backward :
              FFTW3(execute(plan_backward));
      break;
    default:
      io::messages.add("fft type invalid", "Lattice Sum", io::message::critical);
      return;
  }
}

template<typename complex_type>
void GenericMesh<complex_type>::fft(GenericMesh<complex_type>::fft_type type) {
  throw std::runtime_error("FFT of a complex data type not implemented");
}

// make sure you instantiate the template with all types you need here.
// if you fail to do so, linking will fail.
template class GenericMesh<double>;
template class GenericMesh<complex_number>;
template class GenericMesh<math::SymmetricMatrix>;

ParallelMesh::ParallelMesh(unsigned int size, unsigned int arank, unsigned int acache_size) :
      GenericMesh<complex_number>(), mesh_left(NULL), mesh_right(NULL), mesh_tmp(NULL),
      num_threads(size), rank(arank), cache_size(acache_size), slice_width(0),
      slice_end(0){
        DEBUG(15, "Cache size: " << cache_size);
}
      
ParallelMesh::ParallelMesh(unsigned int size, unsigned int arank, unsigned int acache_size,
      unsigned int x, unsigned int y, unsigned int z) : GenericMesh<complex_number>(),
      mesh_left(NULL), mesh_right(NULL), mesh_tmp(NULL),
      num_threads(size), rank(arank), cache_size(cache_size), slice_width(0),
      slice_start(0), slice_end(0) {
  resize(x,y,z);
}
      
ParallelMesh::~ParallelMesh() {
  // destory the plans
  if (plan_forward != NULL)
    FFTW3(destroy_plan(plan_forward));
  if (plan_backward != NULL)
    FFTW3(destroy_plan(plan_backward));
  
  if (m_mesh != NULL) {
    // make sure you free the mesh. It can be big.
    FFTW3(free((configuration::complex_number*) m_mesh));
  }
  if (mesh_left != NULL)
    FFTW3(free(mesh_left));
  if (mesh_right != NULL)
    FFTW3(free(mesh_right));
  if (mesh_tmp != NULL)
    FFTW3(free(mesh_tmp));
}
      
void ParallelMesh::resize(unsigned int x, unsigned int y, unsigned int z) {
#ifdef XXMPI
  DEBUG(14, "ParallelMesh: resize.");
  m_x = x; m_y = y; m_z = z; m_volume = x*y*z;
  m_x2 = m_x / 2;
  
  const unsigned int cache_volume = cache_size * m_y *m_z;
  
  if (mesh_left != NULL)
    FFTW3(free(mesh_left));
  
  if (cache_size != 0)
    mesh_left = (complex_number *) FFTW3(malloc(sizeof(FFTW3(complex)) * cache_volume));
  
  if (mesh_right != NULL)
    FFTW3(free(mesh_right));
  
  if (cache_size != 0)
    mesh_right = (complex_number *) FFTW3(malloc(sizeof(FFTW3(complex)) * cache_volume));
  
  if (mesh_tmp != NULL)
    FFTW3(free(mesh_tmp));
  
  if (cache_size != 0)
    mesh_tmp = (complex_number *) FFTW3(malloc(sizeof(FFTW3(complex)) * cache_volume));
  
  
  // the parallel FFTW library takes care of all the ranges and returns
  // them.
  // local_alloc: number of complex numbers in the slice
  // local_n0: the width of the slice (along x)
  // local_0_start: the index (x) at which the slice starts
  ptrdiff_t local_alloc, local_n0, local_0_start;
  local_alloc = FFTW3(mpi_local_size_3d(m_x, m_y, m_z, MPI::COMM_WORLD,
                                              &local_n0, &local_0_start));
  DEBUG(12,"local_n0: " << local_n0 << " local_0_start" << local_0_start);
  slice_width = local_n0;
  slice_start = local_0_start;
  slice_end = local_0_start + local_n0;
  slice_end_inclusive = slice_end - 1;
  DEBUG(12,"slice width: " << slice_width << " left: " << slice_start <<
          " right: " << slice_end);

  if (slice_width < cache_size) {
    io::messages.add("Cache bigger than a slice of the grid. Reduce number of CPUs.",
            "Lattice Sum", io::message::error);
    return;
  }
  
  if (m_mesh != NULL)
    FFTW3(free(m_mesh));
  
  m_mesh = (complex_number *) FFTW3(malloc(sizeof(FFTW3(complex)) * local_alloc));
  
  if (plan_forward != NULL)
    FFTW3(destroy_plan(plan_forward));
  plan_forward = FFTW3(mpi_plan_dft_3d(m_x, m_y, m_z,
              reinterpret_cast<FFTW3(complex)*> (m_mesh), reinterpret_cast<FFTW3(complex)*> (m_mesh),
              MPI::COMM_WORLD, FFTW_FORWARD, FFTW_ESTIMATE));
  
  if (plan_backward != NULL)
    FFTW3(destroy_plan(plan_backward));
  
  plan_backward = FFTW3(mpi_plan_dft_3d(m_x, m_y, m_z,
              reinterpret_cast<FFTW3(complex)*> (m_mesh), reinterpret_cast<FFTW3(complex)*> (m_mesh),
              MPI::COMM_WORLD, FFTW_BACKWARD, FFTW_ESTIMATE));
#endif
}

void ParallelMesh::zero() {
  // zero the mesh
  const unsigned int slice_volume = slice_width * m_y * m_z;
  for(unsigned int i = 0; i < slice_volume; ++i)
    m_mesh[i] = 0.0;
  // and the caches
  const unsigned int cache_volume = cache_size * m_y * m_z;
  for(unsigned int i = 0; i < cache_volume; ++i) {
    mesh_left[i] = 0.0;
    mesh_right[i] = 0.0;
  }
}

void ParallelMesh::get_neighbors() {
#ifdef XXMPI
  if (cache_size == 0)
    return;
  const unsigned int cache_volume = cache_size * m_y * m_z;
  
  const unsigned int cpu_left = (rank + num_threads - 1) % num_threads;
  const unsigned int cpu_right = (rank + 1) % num_threads;
  
  // create the cache for the left CPU
  for(unsigned int x = 0; x < cache_size; ++x) {
    for(unsigned int y = 0; y < m_y; ++y) {
      for(unsigned int z = 0; z < m_z; ++z) {
        const unsigned int index = z + m_z*(y + m_y*x);
        mesh_tmp[index] = m_mesh[index];
      }
    }
  }
  
  // send it to the left cpu (nonblocking)
  MPI::Request r = MPI::COMM_WORLD.Isend(mesh_tmp, cache_volume * 2, MPI::DOUBLE, cpu_left, 0);
  // receive the right cache (blocking)
  MPI::COMM_WORLD.Recv(mesh_right, cache_volume * 2, MPI::DOUBLE, cpu_right, 0);
  // wait before overwriting the mesh_tmp variable
  r.Wait();
 
  
  // create the cache for the right CPU
  for(unsigned int x = 0; x < cache_size; ++x) {
    const unsigned int mesh_x = slice_width - cache_size + x;
    for(unsigned int y = 0; y < m_y; ++y) {
      for(unsigned int z = 0; z < m_z; ++z) {
        mesh_tmp[z + m_z*(y + m_y*x)] = m_mesh[z + m_z*(y + m_y*mesh_x)];
      }
    }
  }
 
  // send it to the right cpu (nonblocking)
  r = MPI::COMM_WORLD.Isend(mesh_tmp, cache_volume * 2, MPI::DOUBLE, cpu_right, 1);
  // receive the left cache (blocking)
  MPI::COMM_WORLD.Recv(mesh_left, cache_volume * 2, MPI::DOUBLE, cpu_left, 1);
  // wait before overwriting the mesh_tmp variable
  r.Wait();  
#endif
}

void ParallelMesh::add_neighbors_caches() {
  if (cache_size == 0)
    return;
#ifdef XXMPI
  const unsigned int cache_volume = cache_size * m_y * m_z;
  
  const unsigned int cpu_left = (rank + num_threads - 1) % num_threads;
  const unsigned int cpu_right = (rank + 1) % num_threads;
  
  // send the left cache to the left cpu (nonblocking)
  MPI::Request r = MPI::COMM_WORLD.Isend(mesh_left, cache_volume * 2, MPI::DOUBLE, cpu_left, 2);
  // receive the left cache from the right cpu (blocking)
  MPI::COMM_WORLD.Recv(mesh_tmp, cache_volume * 2, MPI::DOUBLE, cpu_right, 2);
  // wait before overwriting the mesh_tmp variable
  r.Wait();

  // add the cache from the right cpu
  for(unsigned int x = 0; x < cache_size; ++x) {
    const unsigned int mesh_x = slice_width - cache_size + x;
    for(unsigned int y = 0; y < m_y; ++y) {
      for(unsigned int z = 0; z < m_z; ++z) {
        m_mesh[z + m_z*(y + m_y*mesh_x)] += mesh_tmp[z + m_z*(y + m_y*x)];
      }
    }
  }

  // send the right cache to the right cpu (nonblocking)
  r = MPI::COMM_WORLD.Isend(mesh_right, cache_volume * 2, MPI::DOUBLE, cpu_right, 3);
  // receive the right cache from the left cpu (blocking)
  MPI::COMM_WORLD.Recv(mesh_tmp, cache_volume * 2, MPI::DOUBLE, cpu_left, 3);
  // wait before overwriting the mesh_tmp variable
  r.Wait();

  // add the cache from the left cpu
  for(unsigned int x = 0; x < cache_size; ++x) {
    for(unsigned int y = 0; y < m_y; ++y) {
      for(unsigned int z = 0; z < m_z; ++z) {
        m_mesh[z + m_z*(y + m_y*x)] += mesh_tmp[z + m_z*(y + m_y*x)];
      }
    }
  }
#endif
}

} // namespace configuration


