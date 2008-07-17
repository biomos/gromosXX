/**
 * @file mesh.cc
 * Mesh class
 */

#include <stdheader.h>

#include <configuration/mesh.h>

#include <util/error.h>
#include <util/debug.h>

#include <fftw3.h>

// for memcpy
#include <cstdlib>

#undef MODULE
#undef SUBMODULE
#define MODULE configuration
#define SUBMODULE configuration

configuration::Mesh::Mesh() : m_x(0), m_y(0), m_z(0), m_volume(0), 
        m_mesh(NULL), plan_forward(NULL), plan_backward(NULL) {
}
configuration::Mesh::Mesh(unsigned int x, unsigned int y, unsigned int z) :
m_x(x), m_y(y), m_z(z), m_volume(x*y*z) {
  // allocate the grid arrays. Here we have to use a dynamic array
  // because static arrays would lie on the stack and are not aligned
  // in the SIMD way. 
  // Elements can be accessed by:
  // statarr[x][y][z] == dynarr[z + Nz * (y + Ny * x)]
  // for more details see 3.2 Multidimensional Array format in section
  
  m_mesh = (complex_type*) fftw_malloc(m_volume * sizeof(fftw_complex));
  plan_forward = fftw_plan_dft_3d(m_x, m_y, m_z,
              reinterpret_cast<fftw_complex*> (m_mesh), reinterpret_cast<fftw_complex*> (m_mesh),
              FFTW_FORWARD, FFTW_ESTIMATE);
  plan_backward = fftw_plan_dft_3d(m_x, m_y, m_z,
              reinterpret_cast<fftw_complex*> (m_mesh), reinterpret_cast<fftw_complex*> (m_mesh),
              FFTW_BACKWARD, FFTW_ESTIMATE);
}

configuration::Mesh::~Mesh() {
  // destory the plans
  if (plan_forward != NULL)
    fftw_destroy_plan(plan_forward);
  if (plan_backward != NULL)
    fftw_destroy_plan(plan_backward);
  
  if (m_mesh != NULL) {
    // make sure you free the mesh. It can be big.
    fftw_free((fftw_complex*) m_mesh);
  }
}

configuration::Mesh::Mesh(const Mesh & mesh) {
  m_x = mesh.x();
  m_y = mesh.y();
  m_z = mesh.z();
  m_volume = m_x * m_y * m_z;
  // reserve memory for the new mesh
  m_mesh = (complex_type*) fftw_malloc(m_volume * sizeof(fftw_complex));
  // just copy the memory
  std::memcpy(m_mesh, mesh.m_mesh, m_volume * sizeof(fftw_complex));
}

void configuration::Mesh::resize(unsigned int x, unsigned int y, unsigned int z){
  m_x = x;
  m_y = y;
  m_z = z;
  m_volume = m_x * m_y * m_z;
  
  if (plan_forward != NULL)
    fftw_destroy_plan(plan_forward);
  if (plan_backward != NULL)
    fftw_destroy_plan(plan_forward);  
  
  if (m_mesh != NULL)
    fftw_free(m_mesh);
  
  m_mesh = (complex_type*) fftw_malloc(m_volume * sizeof(fftw_complex));
  plan_forward = fftw_plan_dft_3d(m_x, m_y, m_z,
              reinterpret_cast<fftw_complex*> (m_mesh), reinterpret_cast<fftw_complex*> (m_mesh),
              FFTW_FORWARD, FFTW_ESTIMATE);
  plan_backward = fftw_plan_dft_3d(m_x, m_y, m_z,
              reinterpret_cast<fftw_complex*> (m_mesh), reinterpret_cast<fftw_complex*> (m_mesh),
              FFTW_BACKWARD, FFTW_ESTIMATE);
}

void configuration::Mesh::fft(configuration::Mesh::fft_type type) {
  switch (type) {
    case configuration::Mesh::fft_forward :
              fftw_execute(plan_forward);
      break;
    case configuration::Mesh::fft_backward :
              fftw_execute(plan_backward);
      break;
    default:
      io::messages.add("fft type invalid", "Lattice Sum", io::message::critical);
      return;
  }
}

void configuration::Mesh::zero() {
  for(unsigned int i = 0; i < m_volume; ++i)
    m_mesh[i] = 0;
}
