/**
 * @file container.h
 * @author Poliak (peter.poliak@boku.ac.at)
 * @date 17.06.2023
 * data container for GPU (pairlist, exclusions, ...)
 * holds data and pointers to the 2D container on GPU
 * this is a jagged array and is a light-weight equivalent of std::vector< std::vector<T> >
 * stored in the GPU memory
 */

#ifndef INCLUDED_CUCONTAINER_H
#define INCLUDED_CUCONTAINER_H

#include "utils.h"

namespace cukernel {
  /**
   * @class Container
   * holds data and pointers to the 2D container on GPU
   * this is a jagged array and is a light-weight equivalent of std::vector< std::vector<T> > for GPU
   */
  template <typename T>
  class Container {
      public:
          typedef T num_type;
          /**
           * default constructor
           */
          Container(bool persistent = true);

          /**
           * allocation constructor
           */
          Container(size_t height, size_t width, bool persistent = true);

          /**
           * destructor
           */
          ~Container() {
              // it is not allowed to call CUDA API in destructors
              // the CUDA context might get destroyed earlier throwing error on normal termination
              // if you really need to destruct the Container, call destruct() explicitly
              // or construct with persistent = false
              if (!m_persistent) {
                  this->deallocate();
              }
          };
          
          /**
           * destruct, free up resources manually
           */
          void destruct() {
              this->deallocate();
          }

          /**
           * resize the container
           */
          void resize(size_t height, size_t width);

          /**
           * get the pointer to the beginning
           */
          __host__ __device__ T* begin() const {
              return this->dev_ptr;
          };

          /**
           * get i-th row
           */
          __host__ __device__ T* operator[](const unsigned i) const {
              assert(i < this->m_height);
              return reinterpret_cast<T*>(reinterpret_cast<char*>(this->dev_ptr) + i * this->dev_pitch);
          };

          /**
           * get i-th row
           */
          __host__ __device__ T* operator()(const unsigned i) const {
              return this->operator[](i);
          };

          /**
           * get the element from i-th row and j-th column
           */
          __host__ __device__ T operator()(const unsigned i, const unsigned j) const;

          /**
           * modify the element from i-th row and j-th column
           */
          __device__ T& operator()(const unsigned i, const unsigned j);

          /**
           * get the height of the array
           */
          __device__ __host__ unsigned height() const {
              return this->m_height;
          }

          /**
           * get the width of the array
           */
          __device__ __host__ unsigned width() const {
              return this->m_width;
          }

          /**
           * add value to the array
           */
          //__device__ int push_back(unsigned row, T val);

          /**
           * add values to the array
           */
          __device__ int push_back(unsigned row, T* val, size_t size = 1);

          /**
           * reserve a strip in a row, return its offset, and increase corresponding size
           * if writing exclusively (no race condition possible), use <true>
           */
          template <bool exclusive = false>
          __device__ int reserve_strip(unsigned i, size_t size);

          /**
           * update the row from host
           */
          int update_row(unsigned row, T* val, unsigned size);

          /**
           * find item in a row
           * 
           * we can maybe distinguish between trivial exclusions (just subsequent indices) and complex (non-consecutive indices)
           */
          __device__ int find(unsigned row, T val) const;

          /**
           * get size of the row
           */
          __device__ __host__ unsigned size(unsigned row) const;

          /**
           * clear all data
           */
          __device__ __host__ void clear();

          /**
           * copy data to host and return as 2D vector
           */
          std::vector< std::vector<T> > copy_to_host() const;

          /**
           * Free up the memory (we have to do this manually as it is too late in destructor)
           * use this method to avoid memory leaks if the instance lifespan is shorter than the program
           */
          void free() {
              this->deallocate();
          }

          /**
           * set overflow value from device
           */
          __device__ void set_overflow() {
              *this->dev_overflow = true;
          }

          /**
           * access the overflow value
           */
          __device__ __host__ bool overflown() const;

          /**
           * print the container content 
           */
          template <typename U>
          friend std::ostream& operator<<(std::ostream &os, const Container<U> &c);

      private:
          /**
           * allocate the memory
           */
          void allocate();

          /**
           * deallocate the memory
           */
          void deallocate();

          /**
           * report (de)allocation only in the debug mode
           */
          template <typename U>
          void report(U* p, std::size_t n, bool alloc = true) const;

          /**
           * reset the overflow flag, only privately, if clearing or resizing
           */
          __device__ __host__ void reset_overflow();

          /**
           * device pointer to the the array
           */
          T* dev_ptr;

          /**
           * pitch on the device (total width of row in bytes)
           */
          size_t dev_pitch;

          /**
           * pointer to the array of sizes of the rows (number of valid elements in every row)
           */
          unsigned* dev_sizes;

          /**
           * overflow indicator
           */
          bool* dev_overflow;

          /**
           * height of the array (number of rows)
           */
          unsigned m_height;

          /**
           * width of the array (number of columns)
           */
          unsigned m_width;

          /**
           * pitch on the host (total width of row in bytes)
           */
          //size_t h_pitch;

          /**
           * persistent, resides in memory for the app lifetime
           */
          bool m_persistent;
  };

  /**
   * print the container content
   */
  template <typename T>
  std::ostream& operator<<(std::ostream &os, const Container<T> &c);
}

#include "container.cu"
#endif
