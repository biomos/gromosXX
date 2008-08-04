/**
 * @file mesh.h
 * the mesh class
 */

#ifndef INCLUDED_MESH_H
#define	INCLUDED_MESH_H
#include <complex>
#include <fftw3.h>
namespace configuration{
  /**
   * @class GenericMesh 
   * a generic gridded quantity
   *
   * all data types have to be manually instantiated in @file mesh.cc
   */
  
  template<typename complex_type>
  class GenericMesh {
  public:
    /**
     * Constructor (empty mesh)
     */
    GenericMesh();
            
    /**
     * Constructor
     */
    GenericMesh(unsigned int size_x, unsigned int size_y, unsigned int size_z);
    /**
     * Destructor
     */
    ~GenericMesh();
    /**
     * copy constructor
     */
    GenericMesh(const GenericMesh<complex_type> & mesh);
    
    /**
     * resize the Mesh.
     * This will delete the current content create a new (non zero) Mesh
     */
    void resize(unsigned int x, unsigned int y, unsigned int z);

    /**
     * accessor to a grid element
     */
    inline complex_type & operator()(unsigned int x, unsigned int y, unsigned int z) {
      assert(x < m_x && y < m_y && z < m_z);
      return m_mesh[z + m_z*(y + m_y * x)];
    }
    /**
     * accessor to grid element
     */
    inline complex_type & element(unsigned int x, unsigned int y, unsigned int z) {
      return (*this)(x,y,z);
    }
    /**
     * const accessor to a grid element
     */
    inline const complex_type & operator()(unsigned int x, unsigned int y, unsigned int z) const {
      assert(x < m_x && y < m_y && z < m_z);
      return m_mesh[z + m_z*(y + m_y * x)];
    }
    /**
     * const accessor to grid element
     */
    inline const complex_type & element(unsigned int x, unsigned int y, unsigned int z) const {
      return (*this)(x,y,z);
    }
    /**
     * accessor to x dimension
     */
    inline unsigned int x() const {
      return m_x;
    }
    /**
     * accessor to y dimension
     */
    inline unsigned int y() const {
      return m_y;
    }
    /**
     * accessor to z dimension
     */
    inline unsigned int z() const {
      return m_z;
    }
    /**
     * accessor to grid volume
     */
    inline unsigned int volume() const {
      return m_volume;
    }
    /**
     * accessor to the mesh
     */
    inline complex_type * mesh() {
      return m_mesh;
    }
    /**
     * const accessor to the mesh
     */
    inline const complex_type * mesh() const {
      return m_mesh;
    }
    /**
     * @enum fft_type
     *  FFT forward or backward
     */
    enum fft_type { fft_forward = 1, fft_backward = -1};
    /**
     * fft the mesh
     * This will result in a runtime error for non complex meshs.
     */
    void fft(fft_type type);
    
    /**
     * accessor to a grid point via an unsigned int vector argument
     * periodicity is taken into account by a modulo operation
     */
    inline complex_type & operator()(const math::GenericVec<int> &p)  {
      // convert it to an unsigned int and do the modulo operation
      // to take the periodicity into account.
      return (*this)(
                (p(0) + m_x) % m_x,
                (p(1) + m_y) % m_y,
                (p(2) + m_z) % m_z);
    }
    
    /**
     * set the numbers to zero
     */
    void zero() {
      for(unsigned int i = 0; i < m_volume; ++i)
        m_mesh[i] = 0;
    }
    
    
  protected:
    /**
     * size of grid in x dimension
     */
    unsigned int m_x;
    /**
     * size of grid in y dimension
     */
    unsigned int m_y;
    /**
     * size of grid in z dimension
     */
    unsigned int m_z;
    /**
     * volume of the grid
     */
    unsigned int m_volume;
    /**
     * the mesh itself
     */
    complex_type * m_mesh;
    /**
     * plan to do a forward fft
     */
    fftw_plan plan_forward;
    /**
     * plan to do a backward fft
     */
    fftw_plan plan_backward;
  };
  
  typedef std::complex<double> complex_number;
  /**
   * @class Mesh
   * the complex mesh used for FFTW calculations.
   */
  typedef GenericMesh<complex_number> Mesh;
}
#endif	/* INCLUDED_MESH_H */


