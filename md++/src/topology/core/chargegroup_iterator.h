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


namespace topology
{
  /**
   * @class Chargegroup_Iterator
   * iterates over the chargegroups.
   */
  class Chargegroup_Iterator
    : public Atomgroup_Iterator<std::vector<int>::const_iterator>
  {
  public:
    Chargegroup_Iterator(std::vector<int>::const_iterator cg_it)
      : Atomgroup_Iterator<std::vector<int>::const_iterator>(cg_it)
    {
    }

    Chargegroup_Iterator()
    {
    }

    void cog(const math::VArray &pos, math::Vec &v)const
    {
      DEBUG(10, "cog: " << *begin() << " - " << *end() << " of " << pos.size());
      v = 0.0;
      for(Atom_Iterator it=begin(), to=end(); it!=to; ++it){
	assert(pos.size() > *it);
	v += pos(int(*it));
      }
      
      v /= num_atoms();
    }
      
  private:

  };
  
}

#endif
  
