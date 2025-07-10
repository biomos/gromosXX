/**
 * @file container.cu
 * @author Poliak (peter.poliak@boku.ac.at)
 * @date 17.06.2023
 * container implementation
 */


#include "gpu.h"

#include "container.h"

template <typename T>
gpu::Container<T>::Container(bool persistent) :
                  dev_ptr(nullptr),
                  dev_pitch(0),
                  dev_sizes(nullptr),
                  dev_overflow(nullptr),
                  m_height(0),
                  m_width(0),
                  m_persistent(persistent) {};

template <typename T>
gpu::Container<T>::Container(size_t height, size_t width, bool persistent) :
                  dev_ptr(nullptr),
                  dev_pitch(0),
                  dev_sizes(nullptr),
                  dev_overflow(nullptr),
                  m_height(height),
                  m_width(width),
                  m_persistent(persistent) {this->allocate();};

/**
* resize the container
*/
template <typename T>
void gpu::Container<T>::resize(size_t height, size_t width) {
    // if increasing the size and is not overflown, copy
    if (height >= this->m_height && width >= this->m_width && !this->overflown()) {
        Container<T> newc(height, width);
        CHECK(cudaMemcpy(newc.dev_sizes, this->dev_sizes, this->m_height*sizeof(unsigned), cudaMemcpyDeviceToDevice));
        CHECK(cudaMemcpy2D(newc.dev_ptr, newc.dev_pitch, this->dev_ptr, this->dev_pitch, this->m_width*sizeof(T), this->m_height, cudaMemcpyDeviceToDevice));
        this->deallocate();
        *this = newc;
    } else {
        this->deallocate();
        this->m_height = height;
        this->m_width = width;
        this->allocate();
    }
};

/**
* get the element from i-th row and j-th column
*/
template <typename T>
__host__ __device__ T gpu::Container<T>::operator()(const unsigned i, const unsigned j) const {
    assert(i < this->m_height);
    assert(j < this->dev_sizes[i]);
    #ifdef __CUDA_ARCH__
        return *(this->operator[](i) + j);
    #else
        // inefficient for __host__ but maybe useful
        T v;
        cudaMemcpy(&v, this->operator[](i) + j, sizeof(T), cudaMemcpyDeviceToHost);
        return v;
#endif
};

/**
 * modify the element from i-th row and j-th column
 */
template <typename T>
__device__ T& gpu::Container<T>::operator()(const unsigned i, const unsigned j) {
    //assert(j < this->m_width);
    assert(j < this->dev_sizes[i]);
    return *(this->operator[](i) + j);
};

/**
* add value to the array
*/
/*__device__ int push_back(unsigned row, T val) {
    assert(row < this->m_height);
    unsigned offset = atomicAdd(&this->dev_sizes[row], 1);
    if (offset >= this->m_width) {
        this->set_overflow();
        return 1;
    }
    this->operator()(row,offset) = val;
    return 0;
}*/

/**
* add value(s) to the array
*/
template <typename T>
__device__ int gpu::Container<T>::push_back(unsigned row, T* val, size_t size) {
    assert(row < this->m_height);
    unsigned offset = atomicAdd(&this->dev_sizes[row], size);
    if (offset + size > this->m_width) {
        this->set_overflow();
        return 1;
    }
    for (unsigned j = 0; j < size; ++j) {
        this->operator()(row,offset+j) = val[j];
    }
    return 0;
}

/**
 * reserve a strip in a row, return its offset, and increase corresponding size
 * if writing exclusively (no race condition possible), use <true>
 */
template <typename T>
template <bool exclusive>
__device__ int gpu::Container<T>::reserve_strip(unsigned i, size_t size) {
    assert(i < this->m_height);
    if (size == 0) return this->dev_sizes[i];
    unsigned offset = 0;
    if (exclusive) {
        offset = this->dev_sizes[i];
        this->dev_sizes[i] += size;
    } else {
        offset = atomicAdd(&this->dev_sizes[i], size);
    }
    if (offset + size > this->m_width) {
        this->set_overflow();
        return -1;
    }
    return offset;
}

/**
 * update the row from host
 */
template <typename T>
int gpu::Container<T>::update_row(unsigned row, T* val, unsigned size) {
    if (size > this->m_width) {
        this->resize(this->m_height, size);
    }
    CHECK(cudaMemcpy(this->dev_sizes + row, &size, sizeof(unsigned), cudaMemcpyHostToDevice));
    CHECK(cudaMemcpy(this->operator[](row), val, size*sizeof(T), cudaMemcpyHostToDevice));
    return 0;
}

/**
 * find item in a row
 */
template <typename T>
__device__ int gpu::Container<T>::find(unsigned row, T val) const {
    const unsigned size = this->dev_sizes[row];
    if (not size) return -1;
    if (val > this->operator()(row, size - 1)) return -1;
    for (unsigned i = 0; i < this->dev_sizes[row]; ++i) {
        if (this->operator()(row, i) == val) {
            return static_cast<int>(i);
        } /*else if (this->operator()(row, i) > val) { // exit if we passed, works only for sorted data
            return -1;
        }*/
    }
    return -1;
}

/**
 * get size of the row
 */
template <typename T>
__device__ __host__ unsigned gpu::Container<T>::size(unsigned row) const {
    unsigned s;
    #ifdef __CUDA_ARCH__
        s = this->dev_sizes[row];
    #else
        CHECK(cudaMemcpy(&s, this->dev_sizes + row, sizeof(unsigned), cudaMemcpyDeviceToHost));
    #endif
    return s;
}

/**
 * get sizes of all the rows
 */
template <typename T>
std::vector<unsigned> gpu::Container<T>::sizes() const {
    std::vector<unsigned> s(this->m_height,0);
    CHECK(cudaMemcpy(s.data(), this->dev_sizes, this->m_height * sizeof(unsigned), cudaMemcpyDeviceToHost));
    return s;
}

/**
 * clear all data
 */
template <typename T>
__host__ __device__ void gpu::Container<T>::clear() {
#ifdef __CUDA_ARCH__
    const unsigned tid = threadIdx.x + blockDim.x * blockIdx.x;
    const unsigned grid_size = blockDim.x * gridDim.x;
    for (unsigned i = tid; i < this->m_height; i+=grid_size) {
        this->dev_sizes[i] = 0;
    }
#else
    #ifndef NDEBUG
        // this does not have to be zeroed
        CHECK(cudaMemset2D(this->dev_ptr, this->dev_pitch, 0, this->m_width * sizeof(T), this->m_height));
    #endif
    CHECK(cudaMemset(this->dev_sizes, 0, this->m_height * sizeof(unsigned)));
#endif
    this->reset_overflow();
};

/**
 * copy data to host and return as 2D vector
 */
template <typename T>
__host__ std::vector< std::vector<T> > gpu::Container<T>::copy_to_host() const {
    std::vector<unsigned> s = this->sizes();

    std::vector< std::vector<T> > arr;
    arr.resize(this->m_height);
    for (unsigned i = 0; i < this->m_height; ++i) {
        arr[i].reserve(this->m_width);
        arr[i].resize(s[i],0);
        CHECK(cudaMemcpy(arr[i].data(), this->operator[](i), s[i] * sizeof(T), cudaMemcpyDeviceToHost));
    }
    return arr;
}

/**
 * access the overflow value
 */
template <typename T>
__device__ __host__ bool gpu::Container<T>::overflown() const {
#ifdef __CUDA_ARCH__
    return *this->dev_overflow;
#else
    bool overflow;
    CHECK(cudaMemcpy(&overflow, this->dev_overflow, sizeof(bool), cudaMemcpyDeviceToHost));
    return overflow;
#endif
}

/**
 * allocate the memory
 */
template <typename T>
void gpu::Container<T>::allocate() {
  if (this->dev_overflow == nullptr) {
    CHECK(cudaMalloc(&this->dev_overflow, sizeof(bool)));
    report(this->dev_overflow, 1);
  }
  CHECK(cudaMallocPitch(&this->dev_ptr,
                        &this->dev_pitch,
                        this->m_width * sizeof(T),
                        this->m_height));
  this->m_width = this->dev_pitch / sizeof(T); // use also the padding area
  report(this->dev_ptr, this->dev_pitch * this->m_height / sizeof(T));
  CHECK(cudaMalloc(&this->dev_sizes, this->m_height * sizeof(unsigned)));
  CHECK(cudaMemset(this->dev_sizes, 0, this->m_height * sizeof(unsigned)));
  report(this->dev_sizes, this->m_height);
  reset_overflow();
}

/**
 * deallocate the memory
 */
template <typename T>
void gpu::Container<T>::deallocate() {
    CHECK(cudaFree(this->dev_ptr));
    //cudaFree(this->dev_ptr);
    report(this->dev_ptr, this->dev_pitch * this->m_height / sizeof(T), false);
    this->dev_ptr = nullptr;
    this->dev_pitch = 0;
    this->m_width = 0;
    CHECK(cudaFree(this->dev_sizes));
    //cudaFree(this->dev_sizes);
    report(this->dev_sizes, this->m_height, false);
    this->dev_sizes = nullptr;
    this->m_height = 0;
    CHECK(cudaFree(this->dev_overflow));
    //cudaFree(this->dev_overflow);
    report(this->dev_overflow, 1, false);
    this->dev_overflow = nullptr;
}


/**
 * report (de)allocation only in the debug mode
 * TODO: set DEBUG level
 */
template <typename T>
template <typename U>
void gpu::Container<T>::report(U* p, std::size_t n, bool alloc) const {
#ifndef NDEBUG
    std::cout << "Container::" << (alloc ? "alloc: " : "dealloc: ") << sizeof(U) * n
        << " bytes at " << std::hex << std::showbase
        << reinterpret_cast<void*>(p) << std::dec << '\n';
#endif
}

/**
 * reset the overflow flag
 */
template <typename T>
__device__ __host__ void gpu::Container<T>::reset_overflow() {
#ifdef __CUDA_ARCH__
    *this->dev_overflow = false;
#else
    CHECK(cudaMemset(this->dev_overflow, 0, sizeof(bool)));
#endif
}

template <typename T>
std::ostream& operator<<(std::ostream &os, const gpu::Container<T> &c) {
    std::vector<unsigned> s;
    s.resize(c.m_height,0);
    CHECK(cudaMemcpy(s.data(), c.dev_sizes, c.m_height * sizeof(unsigned), cudaMemcpyDeviceToHost));

    std::vector<T> arr;
    arr.reserve(c.m_width);
    for (unsigned i = 0; i < c.m_height; ++i) {
        os << i << " (" << s[i] << "): ";
        arr.resize(s[i],0);
        CHECK(cudaMemcpy(arr.data(), c[i], s[i] * sizeof(T), cudaMemcpyDeviceToHost));
        for (unsigned j = 0; j < s[i]; ++j) {
            os << arr[j] << " ";
        }
        os << std::endl;
    }
    return os;
};