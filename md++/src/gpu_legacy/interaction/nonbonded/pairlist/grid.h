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

namespace cukernel {
    namespace interaction {
    template <typename T, typename std::enable_if< // allow only integer types
                                          std::is_integral<T>::value, bool>::type = true>
    class GridT {
      public:
        /**
         * atoms index type
         */
        typedef T num_type;

        /**
         * @brief Construct a new Grid object - rectangular
         * 
         * @param num_atoms total number of atoms
         * @param grid_dims number of cells in the x,y,z axes
         */
        GridT(T num_atoms, dim3 grid_dims);

        /**
         * @brief Construct a new Grid object - quadratic
         * 
         * @param num_atoms total number of atoms
         * @param grid_dim number of cells in the x,y,z axes
         */
        GridT(T num_atoms, unsigned grid_dim);

        /**
         * @brief Construct a new grid object - rectangular
         * 
         * @param num_atoms total number of atoms
         * @param num_cells_x number of cells in the x axis
         * @param num_cells_y number of cells in the y axis
         * @param num_cells_z number of cells in the z axis
         */
        GridT(T num_atoms, unsigned num_cells_x, unsigned num_cells_y, unsigned num_cells_z);

        // disable default constructor
        GridT() = delete;
        
        /**
         * clear all data
         */
        __device__ __host__ void clear();

        /**
         * get sizes of the cells
         */
        typename std::vector<unsigned> sizes() const;

        /**
         * get all atoms in the specified cell
         */
        __host__ typename std::vector<T> cell_atoms(unsigned i) const;

        /**
         * convert single cell index to 3d indices
         */
        __device__ __host__ uint3 unpack_indices(unsigned cell_idx) const;

        /**
         * convert 3d indices to single cell index
         */
        __device__ __host__ unsigned pack_indices(int i, int j, int k) const;
        
        /**
         * convert 3d signed indices to single cell index
         */
        __device__ __host__ unsigned pack_indices(int3 d) const;
        
        /**
         * convert 3d unsigned indices to single cell index
         */
        __device__ __host__ unsigned pack_indices(uint3 d) const;

        /**
         * total number of cells
         */
        __device__ __host__ unsigned num_cells() const;

        /**
         * total number of atoms
         */
        __device__ __host__ T num_atoms() const;

        /**
         * number of atoms in the cell
         */
        __device__ unsigned cell_size(unsigned cell_idx) const;

        /**
         * get cell of the atom
         */
        __device__ unsigned atom(T i) const;

        /**
         * get i-th atom from the cell
         */
        __device__ T get_atom_in_cell(T i, unsigned cell_idx) const;

        /**
         * loop over atoms and assign their cell index
         */
        __device__ unsigned assign_cell_to_atom(T atom_idx, uint3& cell) const;

        /**
         * sort all atoms to their cells and store in the cell grid
         */
        __device__ void store_atom_in_cell(T atom_idx, unsigned cell_idx);

        /**
         * get the pointer to the beginning of the grid
         */
        __device__ T* begin() const;

        /**
         * set the cell pointer to the first atom of the cell
         */
        __device__ void set_cell_pointer(T i, T* ptr) const;

        /**
         * assigns atom to the grid
         * generates linear array of atoms sorted by cells, with pointers and sizes to each of them
         * num_threads has to be tuned, for smaller systems, 32 seems optimal
         */
        template <unsigned num_threads>
        __host__ int assign_atoms_to_grid(float3* pos);

      private:
        /**
         * allocate the GPU memory
         */
        void allocate();

        /**
         * deallocate the GPU memory
         */
        void deallocate();

        /**
         * report (de)allocation (only in the debug mode)
         */
        template <typename U>
        void report(U* p, std::size_t n, bool alloc = true) const;
        
        /**
         * dimensions of the grid
         */
        const dim3 m_grid_dim;
        
        /**
         * actual storage array containing atom indices ordered by correspondence to the cell
         */
        T *dev_grid;
        
        /**
         * pointers to beginning of every cell in dev_grid
         */
        T **dev_ptrs;

        /**
         * sizes of cells in dev_grid
         */
        unsigned *dev_sizes;
        
        /**
         * array assigning atoms to cells
         */
        unsigned *dev_cells;

        /**
         * number of atoms
         */
        const T m_num_atoms;

        /**
         * number of cells
         */
        const unsigned m_num_cells;
    };

    /**
     * assign atoms to cells
     */
    template <typename T>
    __global__ void assign_cells(GridT<T> grid, float3 *pos);
    
    /**
     * set grid cell pointers
     */
    template <typename T>
    __global__ void set_cell_pointers(const GridT<T> grid);

    /**
     * fill dev_grid with atoms grouped in cells
     */
    template <typename T>
    __global__ void group_cells(GridT<T> grid);

    typedef GridT<unsigned> Grid;
  }
}

#include "grid.cu"
