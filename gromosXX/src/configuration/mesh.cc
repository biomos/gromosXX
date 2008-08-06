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
  
  m_mesh = (complex_number*) fftw_malloc(m_volume * sizeof(fftw_complex));
  plan_forward = fftw_plan_dft_3d(m_x, m_y, m_z,
              reinterpret_cast<fftw_complex*> (m_mesh), reinterpret_cast<fftw_complex*> (m_mesh),
              FFTW_FORWARD, FFTW_ESTIMATE);
  plan_backward = fftw_plan_dft_3d(m_x, m_y, m_z,
              reinterpret_cast<fftw_complex*> (m_mesh), reinterpret_cast<fftw_complex*> (m_mesh),
              FFTW_BACKWARD, FFTW_ESTIMATE);
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
    fftw_destroy_plan(plan_forward);
  if (plan_backward != NULL)
    fftw_destroy_plan(plan_backward);
  
  if (m_mesh != NULL) {
    // make sure you free the mesh. It can be big.
    fftw_free((configuration::complex_number*) m_mesh);
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
  m_mesh = (configuration::complex_number*) fftw_malloc(m_volume * sizeof(configuration::complex_number));
  // just copy the memory
  std::memcpy(m_mesh, mesh.m_mesh, m_volume * sizeof(configuration::complex_number));
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
    fftw_destroy_plan(plan_forward);
  if (plan_backward != NULL)
    fftw_destroy_plan(plan_forward);  
  
  if (m_mesh != NULL)
    fftw_free(m_mesh);
  
  m_mesh = (complex_number*) fftw_malloc(m_volume * sizeof(fftw_complex));
  plan_forward = fftw_plan_dft_3d(m_x, m_y, m_z,
              reinterpret_cast<fftw_complex*> (m_mesh), reinterpret_cast<fftw_complex*> (m_mesh),
              FFTW_FORWARD, FFTW_ESTIMATE);
  plan_backward = fftw_plan_dft_3d(m_x, m_y, m_z,
              reinterpret_cast<fftw_complex*> (m_mesh), reinterpret_cast<fftw_complex*> (m_mesh),
              FFTW_BACKWARD, FFTW_ESTIMATE);
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
              fftw_execute(plan_forward);
      break;
    case GenericMesh<complex_number>::fft_backward :
              fftw_execute(plan_backward);
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
}

// make sure you instantiate the template with all types you need here.
// if you fail to do so, linking will fail.
template class configuration::GenericMesh<double>;
template class configuration::GenericMesh<configuration::complex_number>;
template class configuration::GenericMesh<math::Matrix>;

