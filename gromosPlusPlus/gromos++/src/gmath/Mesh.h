/*
 * This file is part of GROMOS.
 * 
 * Copyright (c) 2011, 2012, 2016, 2018, 2021, 2023 Biomos b.v.
 * See <https://www.gromos.net> for details.
 * 
 * GROMOS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

/**
 * @file Mesh.h
 * a 3D mesh class
 */

#ifndef INCLUDED_GMATH_MESH_H
#define	INCLUDED_GMATH_MESH_H

#include <vector>

#include "Vec.h"
#include "../gcore/Box.h"

namespace gmath
{
  /**
   * A generic 3D mesh class which can be used for spatial distributions
   * and can be written in formats readable by VMD
   *
   * @class Mesh
   * @author N. Schmid
   * @ingroup gmath
   */
  template<typename T>
  class Mesh {
    std::vector<T> data;
    int m_size[3];
    std::string m_title;
    gcore::Box m_box;
    double K_inv2, L_inv2, M_inv2;
    T dummy;
  public:
    /**
     * default constructor
     */
    Mesh() : m_title("GROMOS Mesh"), K_inv2(1.0), L_inv2(1.0), M_inv2(1.0) { clear(); }

    /**
     * constructor
     */
    Mesh(int x, int y, int z) : m_title("GROMOS Mesh"), K_inv2(1.0), L_inv2(1.0), M_inv2(1.0) { resize(x,y,z); }

    /**
     * copy constructor
     * @param mesh the mesh to copy
     */
    Mesh(const Mesh<T> & mesh);
    /**
     * assignment copy
     * @param mesh the mesh to copy
     * @return reference to the mesh
     */
    Mesh<T> & operator=(const Mesh<T> & mesh) {
      if (this != &mesh) {
        new(this) Mesh(mesh);
      }
      return *this;
    }

    /**
     * assignment of a single value
     * @param val the value to assign
     * @return reference to the mesh
     */
    Mesh<T> & operator=(const T & val) {
      for (typename std::vector<T>::iterator it = data.begin(), to = data.end();
              it != to; ++it) {
        *it = val;
      }
      return *this;
    }



    /**
     * clear the mesh, free the data, set the size to zero
     */
    void clear() {
      data.clear();
      m_size[0] = m_size[1] = m_size[2] = 0;
    }

    /**
     * accessor to size
     * @return a three membered array containing the dimensions along x, y and z.
     */
    const int* size() const {
      return m_size;
    }

    /**
     * number of grid points
     * @return the number of grid points
     */
    int numPoints() const {
      return int(data.size());
    }

    /**
     * accessor to the title
     * @return the title of the mesh
     */
    const std::string & title() const {
      return m_title;
    }

    /**
     * setter of the title
     * @param title the new title
     */
    void setTitle(const std::string & title) {
      m_title = title;
    }

    /**
     * the the dimensionality of the grid
     * @param x dimension along x
     * @param y dimension along y
     * @param z dimension along z
     */
    void resize(int x, int y, int z) {
      data.resize(x * y * z);
      m_size[0] = x;
      m_size[1] = y;
      m_size[2] = z;
    }

    /**
     * set the box
     * @param box the box
     */
    void setBox(const gcore::Box & box) {
      m_box = box;
      K_inv2 = 1.0 / m_box.K().abs2();
      L_inv2 = 1.0 / m_box.L().abs2();
      M_inv2 = 1.0 / m_box.M().abs2();
    }

    /**
     * access a data element (const version)
     * @param x x index
     * @param y y index
     * @param z z index
     * @return the data element at position (x,y,z)
     */
    inline const T & at(int x, int y, int z) const {
      assert(x < m_size[0] && y < m_size[1] && z < m_size[2]);
      return data[z + m_size[2]*(y + m_size[1] * x)];
    }
    /**
     * access a data element (const version)
     * @param i  index
     */
    inline const T & at(unsigned int i) const {
      assert(i < data.size());
      return data[i];
    }
    /**
     * access a data element (const version)
     * @param x x index
     * @param y y index
     * @param z z index
     * @return the data element at position (x,y,z)
     */
    inline const T & operator()(int x, int y, int z) const {
      return at(x, y, z);
    }

    /**
     * access a data element
     * @param x x index
     * @param y y index
     * @param z z index
     * @return the data element at position (x,y,z)
     */
    inline T & operator()(int x, int y, int z) {
      return const_cast<T&>(at(x, y, z));
    }

    /**
     * access a data element (const version) at a position in space.
     * A box is needed for this function. Set that first.
     * @param v position
     * @return the data element at position v
     */
    inline const T & at(const gmath::Vec & v) const {
      int x = int(m_box.K().dot(v) * m_size[0] * K_inv2);
      int y = int(m_box.L().dot(v) * m_size[1] * L_inv2);
      int z = int(m_box.M().dot(v) * m_size[2] * M_inv2);

      if (x < 0 || x >= m_size[0] ||
              y < 0 || y >= m_size[1] ||
              z < 0 || z >= m_size[2])
        return dummy;

      return (*this)(x, y, z);
    }
    /**
     * access a data element (const version) at a position in space.
     * A box is needed for this function. Set that first.
     * @param v position
     * @return the data element at position v
     */
    inline const T & operator()(const gmath::Vec & v) const {
      return at(v);
    }

    /**
     * access a data element at a position in space.
     * A box is needed for this function. Set that first.
     * @param v position
     * @return the data element at position v
     */
    inline T & operator()(const gmath::Vec & v) {
      return const_cast<T&>(at(v));
    }

    /**
     * access a data element (const version)
     * @param i index
     * @return the data element at index i
     */
    inline const T & operator()(int i) const {
      return at(i);
    }

    /**
     * access a data element
     * @param i index
     * @return the data element at index i
     */
    inline T & operator()(int i) {
      return const_cast<T&>(at(i));
    }

    /**
     * Get the position of a grid point depending on the box.
     * @param x x coordinate of the grid point
     * @param y y coordinate of the grid point
     * @param z z coordinate of the grid point
     * @return position of the grid point
     */
    inline gmath::Vec pos(int x, int y, int z) {
      return m_box.K() * (double(x) / m_size[0]) +
              m_box.L() * (double(y) / m_size[1]) +
              m_box.M() * (double(z) / m_size[2]);
    }

    /**
     * write the mesh to a file in Gaussians CUBE format.
     * A box is required for this.
     * In addition one dummy hydrogen atom at the origin will be added.
     * @param os the stream to write the stuff to.
     */
    void write(std::ostream & os) const;

  };
}

#endif	/* INCLUDED_GMATH_MESH_H */

