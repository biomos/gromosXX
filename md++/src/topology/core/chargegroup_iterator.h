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
 * @file chargegroup_iterator.h
 * iterator over the chargegroups.
 */

#ifndef INCLUDED_CHARGEGROUP_ITERATOR_H
#define INCLUDED_CHARGEGROUP_ITERATOR_H

#undef MODULE
#undef SUBMODULE
#define MODULE topology
#define SUBMODULE topology

#include "gpu/cuda/cuhostdevice.h"


namespace topology
{
  /**
   * @class Chargegroup_Iterator
   * iterates over the chargegroups.
   */
  template <template <typename...> class VecType = std::vector>
  class Chargegroup_IteratorT
    : public Atomgroup_Iterator<typename VecType<int>::const_iterator>
  {
  public:
    using Base = Atomgroup_Iterator<typename VecType<int>::const_iterator>;

    HOSTDEVICE Chargegroup_IteratorT(typename VecType<int>::const_iterator cg_it)
      : Base(cg_it)
    {
    }

    HOSTDEVICE Chargegroup_IteratorT()
    {
    }

  //   void cog(const math::VArray &pos, math::Vec &v)const
  //   {
  //     DEBUG(10, "cog: " << *this->begin() << " - " << *this->end() << " of " << pos.size());
  //     v = 0.0;
  //     for(Atom_Iterator it=this->begin(), to=this->end(); it!=to; ++it){
	// assert(pos.size() > *it);
	// v += pos(int(*it));
  //     }
      
  //     v /= this->num_atoms();
  //   }

    template <typename VArrayType, typename VType>
    HOSTDEVICE void cog(const VArrayType &pos, VType &v)const
    {
      // Check that pos(int) is valid and yields a type we can add to v
      static_assert(
          std::is_arithmetic<decltype(pos(int{}))>::value ||
          std::is_same<decltype(v += pos(int{})), VType&>::value,
          "VArrayType must support pos(int) returning a type compatible with VType"
      );

      // Check that v supports /= with size_t (num_atoms())
      static_assert(
          std::is_arithmetic<VType>::value || 
          std::is_same<decltype(v /= size_t{}), VType&>::value,
          "VType must support /= with size_t"
      );

      DEBUG(10, "cog: " << *this->begin() << " - " << *this->end() << " of " << pos.size());
      v = 0.0;
      for(Atom_Iterator it=this->begin(), to=this->end(); it!=to; ++it){
        assert(pos.size() > *it);
        v += pos(int(*it));
      }
      
      v /= this->num_atoms();
    }
      
  private:

  };

  using Chargegroup_Iterator = Chargegroup_IteratorT<std::vector>;
}

#endif
  
