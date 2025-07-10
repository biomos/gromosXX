/**
 * @file grid.h
 * @author Poliak (peter.poliak@boku.ac.at)
 * @brief Container for atom grid cell assignment
 * @version 0.1
 * @date 2023-06-22
 * 
 * @copyright Copyright (c) 2023
 * 
 */

#pragma once

namespace cuinteraction = cukernel::interaction;

template <typename T>
cuinteraction::GridT<T>::GridT(T num_atoms, dim3 grid_dims) :  m_grid_dim(grid_dims),
                                      dev_grid(nullptr),
                                      dev_ptrs(nullptr),
                                      dev_sizes(nullptr),
                                      dev_cells(nullptr),
                                      m_num_atoms(num_atoms),
                                      m_num_cells(grid_dims.x * grid_dims.y * grid_dims.z) {
    this->allocate();
#ifndef NDEBUG
    std::cout << "m_grid_dim: " << m_grid_dim.x << " " << m_grid_dim.y << " " << m_grid_dim.z << std::endl;
    std::cout << "dev_grid: " << dev_grid << std::endl;
    std::cout << "dev_ptrs: " << dev_ptrs << std::endl;
    std::cout << "dev_sizes: " << dev_sizes << std::endl;
    std::cout << "dev_cells: " << dev_cells << std::endl;
    std::cout << "m_num_atoms: " << m_num_atoms << std::endl;
    std::cout << "m_num_cells: " << m_num_cells << std::endl;
#endif
}

template <typename T>
cuinteraction::GridT<T>::GridT(T num_atoms, unsigned grid_dim) :
                                      GridT(num_atoms, make_uint3(grid_dim, grid_dim, grid_dim)) {};

template <typename T>
cuinteraction::GridT<T>::GridT(T num_atoms, unsigned num_cells_x, unsigned num_cells_y, unsigned num_cells_z) :
                                      GridT(num_atoms, make_uint3(num_cells_x, num_cells_y, num_cells_z)) {};
      
template <typename T>
__device__ __host__ void cuinteraction::GridT<T>::clear() {
#ifdef __CUDA_ARCH__
    const unsigned tid = threadIdx.x + blockDim.x * blockIdx.x;
    const unsigned grid_size = blockDim.x * gridDim.x;
    for (unsigned i = tid; i < this->num_cells(); i+=grid_size) {
        this->dev_ptrs[i] = nullptr;
        this->dev_sizes[i] = 0;
    }
#else
    // this does not have to be zeroed
    //CHECK(cudaMemset(this->dev_grid, 0, this->m_num_atoms * sizeof(T)));
    CHECK(cudaMemset(this->dev_sizes, 0, this->num_cells() * sizeof(T*)));
    CHECK(cudaMemset(this->dev_ptrs, nullptr, this->num_cells() * sizeof(unsigned)));
#endif
};


template <typename T>
typename std::vector<unsigned> cuinteraction::GridT<T>::sizes() const {
  std::vector<unsigned> sizes(this->m_num_cells);
  CHECK(cudaMemcpy(sizes.data(), this->dev_sizes, this->m_num_cells*sizeof(unsigned), cudaMemcpyDeviceToHost));
  return sizes;
}

template <typename T>
typename std::vector<T> cuinteraction::GridT<T>::cell_atoms(unsigned i) const {
  unsigned v_size;
  CHECK(cudaMemcpy(&v_size, this->dev_sizes + i, sizeof(unsigned), cudaMemcpyDeviceToHost));
  T* cell_ptr;
  CHECK(cudaMemcpy(&cell_ptr, this->dev_ptrs + i, sizeof(T*), cudaMemcpyDeviceToHost));
  std::vector<T> v(v_size);
  CHECK(cudaMemcpy(v.data(), cell_ptr, v_size * sizeof(T), cudaMemcpyDeviceToHost));
  return v;
}

template <typename T>
__device__ __host__ uint3 cuinteraction::GridT<T>::unpack_indices(unsigned cell_idx) const {
  const unsigned k = cell_idx / (this->m_grid_dim.y * this->m_grid_dim.x);
  const unsigned tmp = cell_idx % (this->m_grid_dim.y * this->m_grid_dim.x);
  const unsigned j = tmp / this->m_grid_dim.x;
  const unsigned i = tmp % this->m_grid_dim.x;
  return uint3{i,j,k};
}

template <typename T>
__device__ __host__ unsigned cuinteraction::GridT<T>::pack_indices(int i, int j, int k) const {
  const int dim_x = this->m_grid_dim.x;
  const int dim_y = this->m_grid_dim.y;
  const int dim_z = this->m_grid_dim.z;
  // wrap around for both positive and negative
  i = (i%dim_x + dim_x) % dim_x;
  j = (j%dim_y + dim_y) % dim_y;
  k = (k%dim_z + dim_z) % dim_z;
  return i + j*dim_x + k*dim_x*dim_y;
}

template <typename T>
__device__ __host__ unsigned cuinteraction::GridT<T>::pack_indices(int3 d) const {
  return this->pack_indices(d.x, d.y, d.z);
}

template <typename T>
__device__ __host__ unsigned cuinteraction::GridT<T>::pack_indices(uint3 d) const {
  assert(d.x < this->m_grid_dim.x);
  assert(d.y < this->m_grid_dim.y);
  assert(d.z < this->m_grid_dim.z);
  return this->pack_indices(d.x, d.y, d.z);
}

template <typename T>
__device__ __host__ unsigned cuinteraction::GridT<T>::num_cells() const {
    return this->m_num_cells;
}

template <typename T>
__device__ __host__ T cuinteraction::GridT<T>::num_atoms() const {
    return this->m_num_atoms;
}

template <typename T>
__device__ unsigned cuinteraction::GridT<T>::cell_size(unsigned cell_idx) const {
    //printf("cell (%u) size: %u\n", cell_idx, this->dev_sizes[cell_idx]);
    return this->dev_sizes[cell_idx];
}

template <typename T>
__device__ unsigned cuinteraction::GridT<T>::atom(T i) const {
    return this->dev_cells[i];
}

template <typename T>
__device__ T cuinteraction::GridT<T>::get_atom_in_cell(T i, unsigned cell_idx) const {
    assert(cell_idx < this->m_num_cells);
    assert(i < this->dev_sizes[cell_idx]);
    return *(this->dev_ptrs[cell_idx] + i);
}

template <typename T>
__device__ unsigned cuinteraction::GridT<T>::assign_cell_to_atom(T atom_idx, uint3& cell) const {
    unsigned cell_idx = this->pack_indices(cell);
    this->dev_cells[atom_idx] = cell_idx;
    return atomicAdd(this->dev_sizes + cell_idx, 1);
}

template <typename T>
__device__ void cuinteraction::GridT<T>::store_atom_in_cell(T atom_idx, unsigned cell_idx) {
  unsigned offset = atomicAdd(this->dev_sizes + cell_idx, 1);
  *(this->dev_ptrs[cell_idx] + offset) = atom_idx; 
}

template <typename T>
__device__ T* cuinteraction::GridT<T>::begin() const {
  return this->dev_grid;
}

template <typename T>
__device__ void cuinteraction::GridT<T>::set_cell_pointer(T i, T* ptr) const {
  this->dev_ptrs[i] = ptr;
}

template <typename T>
template <unsigned num_threads>
int cuinteraction::GridT<T>::assign_atoms_to_grid(float3* pos) {
    // iterate over atoms
    const unsigned a_num_blocks = this->m_num_atoms / num_threads + 1;
    // iterate over cells
    const unsigned c_num_blocks = this->m_num_cells / num_threads + 1;
    // const unsigned num_threads = 1;
    // const unsigned a_num_blocks = 1, c_num_blocks = 1;
    cudaMemset(this->dev_cells, 0, this->m_num_atoms * sizeof(unsigned));
    cudaMemset(this->dev_sizes, 0, this->m_num_cells * sizeof(unsigned));
    assign_cells<<<a_num_blocks, num_threads>>>(*this, pos);
    set_cell_pointers<<<c_num_blocks, num_threads>>>(*this);
    // We zero the sizes to track number of written elements in next kernel
    CHECK(cudaMemset(this->dev_sizes, 0, this->m_num_cells * sizeof(unsigned)));
    group_cells<<<a_num_blocks, num_threads>>>(*this);
    // Should we also sort the atoms in cells? 
    // TODO: Try
    cudaDeviceSynchronize();
    return 0;
}


template <typename T>
void cuinteraction::GridT<T>::allocate() {
    CHECK(cudaMalloc(&this->dev_grid, this->m_num_atoms*sizeof(T)));
    report(this->dev_grid, this->m_num_atoms);
    CHECK(cudaMalloc(&this->dev_ptrs, this->m_num_cells*sizeof(T*)));
    report(this->dev_ptrs, this->m_num_cells);
    CHECK(cudaMalloc(&this->dev_sizes, this->m_num_cells*sizeof(unsigned)));
    report(this->dev_sizes, this->m_num_cells);
    CHECK(cudaMalloc(&this->dev_cells, this->m_num_atoms*sizeof(T)));
    report(this->dev_cells, this->m_num_atoms);
}

template <typename T>
void cuinteraction::GridT<T>::deallocate() {
    CHECK(cudaFree(&this->dev_grid));
    report(this->dev_grid, this->m_num_atoms, false);
    this->dev_grid = nullptr;
    CHECK(cudaFree(&this->dev_ptrs));
    report(this->dev_ptrs, this->m_num_cells, false);
    this->dev_ptrs = nullptr;
    CHECK(cudaFree(&this->dev_sizes));
    report(this->dev_sizes, this->m_num_cells, false);
    this->dev_sizes = nullptr;
    CHECK(cudaFree(&this->dev_cells));
    report(this->dev_cells, this->m_num_atoms, false);
    this->dev_cells = nullptr;
    this->m_num_atoms = 0;
    this->m_num_cells = 0;
}

template <typename T>
template <typename U>
void cuinteraction::GridT<T>::report(U* p, std::size_t n, bool alloc) const {
#ifndef NDEBUG
    std::cout << "Grid::" << (alloc ? "alloc: " : "dealloc: ") << sizeof(U) * n
                << " bytes at " << std::hex << std::showbase
                << reinterpret_cast<void*>(p) << std::dec << '\n';
#endif
}

template <typename T>
__global__ void cuinteraction::::assign_cells(cukernel::GridT<T> grid, float3 *pos) {
    unsigned tid = threadIdx.x + blockIdx.x * blockDim.x;
    unsigned grid_size = blockDim.x * gridDim.x;
    for (T idx = tid; idx < grid.num_atoms(); idx += grid_size) {
    float3 p = pos[idx];
    // Move to first quadrant box
    while (p.x < 0.0) {
    p.x += dev_box_size.x;
    }
    while (p.x >= dev_box_size.x) {
    p.x -= dev_box_size.x;
    }

    while (p.y < 0.0) {
    p.y += dev_box_size.y;
    }
    while (p.y >= dev_box_size.y) {
    p.y -= dev_box_size.y;
    }

    while (p.z < 0.0) {
    p.z += dev_box_size.z;
    }
    while (p.z >= dev_box_size.z) {
    p.z -= dev_box_size.z;
    }

    uint3 cell;
    cell.x = p.x / dev_cell_size.x;
    cell.y = p.y / dev_cell_size.y;
    cell.z = p.z / dev_cell_size.z;
    grid.assign_cell_to_atom(idx, cell);
    }
};

template <typename T>
__global__ void cuinteraction::::set_cell_pointers(const cukernel::GridT<T> grid) {
    unsigned tid = threadIdx.x + blockIdx.x * blockDim.x;
    unsigned grid_size = blockDim.x * gridDim.x;
    for (T idx = tid; idx < grid.num_cells(); idx += grid_size) {
        T* cell_ptr = grid.begin();
        for(unsigned i = 0; i < idx; ++i) {
            cell_ptr += grid.cell_size(i);
        }
        //printf("cell_ptr: 0x%p\n", cell_ptr);
        grid.set_cell_pointer(idx, cell_ptr);
    }
};

template <typename T>
__global__ void cuinteraction::::group_cells(cukernel::GridT<T> grid) {
    unsigned tid = threadIdx.x + blockIdx.x * blockDim.x;
    unsigned grid_size = blockDim.x * gridDim.x;
    for (T idx = tid; idx < grid.num_atoms(); idx += grid_size) {
        unsigned cell_idx = grid.atom(idx);
        grid.store_atom_in_cell(idx, cell_idx);
    }
};
