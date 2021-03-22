/**
 * @file mesh.h
 * the mesh class
 */

#ifndef INCLUDED_MESH_H
#define	INCLUDED_MESH_H
#include "../math/fft.h"

#undef MODULE
#undef SUBMODULE
#define MODULE configuration
#define SUBMODULE configuration

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
    /**
     * get the boundaries. Only for parallelization
     */
    unsigned int left_boundary() const {
      return 0;
    }
    /**
     * get the boundaries. Only for parallelization
     */
    unsigned int right_boundary() const {
      return m_x;
    }
    /**
     * get the data from the other processors and fill the caches. Only
     * for compatibility.
     */
    void get_neighbors() { return; }
    /**
     * get the caches from the other processors and add them
     * to the real data. Only for compatibility.
     */
    void add_neighbors_caches() {
      return;
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
    FFTW3(plan) plan_forward;
    /**
     * plan to do a backward fft
     */
    FFTW3(plan) plan_backward;
  };
  
  typedef std::complex<double> complex_number;
  /**
   * @class Mesh
   * the complex mesh used for FFTW calculations.
   */
  typedef GenericMesh<complex_number> Mesh;
  
  /**
   * @class ParallelMesh
   * A mesh for MPI parallelization
   */
  class ParallelMesh : public GenericMesh<complex_number> {
  public:
    /**
     * constructor
     */
    ParallelMesh(unsigned int size, unsigned int arank, unsigned int acache_size);
    /**
     * another constructor
     */
    ParallelMesh(unsigned int size, unsigned int arank, unsigned int acache_size,
            unsigned int x, unsigned int y, unsigned int z);
    /**
     * resize the mesh
     */
    void resize(unsigned int x, unsigned int y, unsigned int z);
    /**
     * destruct the mesh
     */
    ~ParallelMesh();
    /**
     * zero the mesh and cache
     */
    void zero();
    /** 
     * access an element from the grid
     * we do no error checking in here. Make sure that you don't access elements
     * that are outside the caches.
     */
    inline complex_number & operator() (unsigned int x, unsigned int y, unsigned int z) {
      assert(x < m_x && y < m_y && z < m_z);
      
      const int l = slice_start;
      const int & r = slice_end_inclusive;
      
      int l2 = (l <= m_x2) ? l + m_x : l - m_x;
      int r2 = (r <= m_x2) ? r + m_x : r - m_x;
      
      int dist_l = std::min(abs(l-int(x)), abs(l2-int(x)));
      int dist_r = std::min(abs(int(x)-r), abs(int(x)-r2));

      if (dist_l + dist_r <= int(slice_width)) {
        DEBUG(30, "rank("<<rank<<"): middle");
        return m_mesh[z + m_z*(y + m_y * (x - slice_start))];
      } else if (dist_l < dist_r) {
        DEBUG(30, "rank(" << rank << "): left");
        // x is just to the left of the slice
        return mesh_left[z + m_z * (y + m_y * (cache_size - dist_l))];
      } else {
        DEBUG(30, "rank(" << rank << "): right");
        // x is just to the right of the slice
        return mesh_right[z + m_z * (y + m_y * (dist_r - 1))];
      }
    }
    /**
     * accessor to a grid point via an unsigned int vector argument
     * periodicity is taken into account by a modulo operation
     */
    inline complex_number & operator()(const math::GenericVec<int> &p)  {
      // convert it to an unsigned int and do the modulo operation
      // to take the periodicity into account.
      return (*this)(
                (p(0) + m_x) % m_x,
                (p(1) + m_y) % m_y,
                (p(2) + m_z) % m_z);
    }
    /**
     * get the data from the other processors and fill the caches
     */
    void get_neighbors();
    /**
     * get the caches from the other processors and add them
     * to the real data
     */
    void add_neighbors_caches();
    /**
     * accessor to the left slice boundary
     */
    unsigned int left_boundary() const {
      return slice_start;
    }
    /**
     * accessor to the left slice boundary
     */
    unsigned int right_boundary() const {
      return slice_end;
    }
  protected :
    /**
     * the mesh to the leftern side of the real mesh
     */
    complex_number *mesh_left;
    /**
     * the mesh to the righern side of the real mesh
     */
    complex_number *mesh_right;
    /**
     * a temporary mesh
     */
    complex_number *mesh_tmp;
    /**
     * number of threads/cpus
     */
    unsigned int num_threads;
    /**
     * rank of the thread
     */
    unsigned int rank;
    /**
     * size of the cach to the left and the right side of the mesh
     */
    unsigned int cache_size;
    /**
     * slice width
     */
    unsigned int slice_width;
    /**
     * integer where slice starts (inclusive)
     */
    int slice_start;
    /**
     * integer where slice ends (exclusive)
     */
    int slice_end;
    /**
     * integer where slice ends (inclusive)
     */
    int slice_end_inclusive;
    /**
     * half of the box length along x
     */
    int m_x2;
  };
}
#endif	/* INCLUDED_MESH_H */


