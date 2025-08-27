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
 * @file periodicity_gpu.cu
 * implementation of the periodic boundary condition functions.
 */

#include "stdheader.h"
#include "gpu/cuda/cuheader.h"

#include "math/gmath.h"
#include "math/boundary_implementation.h"

#include "algorithm/algorithm.h"

#include "gpu/cuda/memory/topology_struct.h"
#include "gpu/cuda/memory/configuration_struct.h"

#include "topology/topology.h"

#include "configuration/configuration.h"

#include "periodicity_gpu.h"


#include "gpu/cuda/kernels/hello_world.h"
#include "gpu/cuda/kernels/periodicity.h"

#undef MODULE
#undef SUBMODULE
#define MODULE math
#define SUBMODULE math

#define NUM_THREADS_PER_BLOCK 256

template<math::boundary_enum b>
gpu::PeriodicityGpu<b>::PeriodicityGpu(math::Box const & bb) 
  : math::Boundary_Implementation<b>(bb)
{
}

template<math::boundary_enum b>
template<typename VecType>
void gpu::PeriodicityGpu<b>::put_into_box(VecType &v)const
{
  VecType o(0, 0, 0);
  this->nearest_image(v, o, v);
}

template<math::boundary_enum b>
void gpu::PeriodicityGpu<b>::put_into_positive_box(math::Vec &v)const {
  //this is not good in the triclinic case...  
  //Vec o(abs(this->m_box(0)), abs(this->m_box(1)), abs(this->m_box(2)));
  //o /= 2;
  math::Vec o = this->m_box(0) + this->m_box(1) + this->m_box(2);
  o /= 2;
 
  this->nearest_image(v, o, v);
  v += o;
}

template <math::boundary_enum b>
void gpu::PeriodicityGpu<b>
::put_chargegroups_into_box(configuration::Configuration & conf, 
			    topology::Topology const & topo)const
{
    // iterate over chargegroups
    // gpu kernels cannot do iterators
    // A. single thread treats single chargegroup ?
    // B. a warp treats single chargegroup ?
    // C. a warp treats 4 chargegroups ?

    /**
     * do this:
     * 1. calculate centers of geometry
     *  - we pull pos of our atoms to shared memory
     *  - we can use a single transaction, since the atoms are adjacent
     *  - we calculate cog of every cg
     * 2. put them into box
     * 3. translate all atoms of chargegroup together
     *  - write from shared to global
     * 4. repeat for solvent (cog is the first atom)
     */


    dim3 dimGrid(1);
    dim3 dimBlock(NUM_THREADS_PER_BLOCK);

    conf.copy_to_gpu();

    gpu::hello_world<<<dimGrid, dimBlock>>>(topo.get_gpu_view(), conf.get_gpu_view());
    gpu::put_chargegroups_into_box<<<dimGrid, dimBlock>>>(topo.get_gpu_view(), conf.get_gpu_view());
}

template<math::boundary_enum b>
void gpu::PeriodicityGpu<b>
::put_chargegroups_into_box_saving_shifts(configuration::Configuration & conf, 
			    topology::Topology const & topo)const
{
  math::VArray &pos = conf.current().pos;
  math::VArray &shift = conf.special().lattice_shifts;
  math::Vec v, v_box, trans;
  
  const math::Box & my_box = this->box();
  math::Matrix L(my_box(0),my_box(1),my_box(2), true);
  const math::Matrix & cartesian_to_oblique = math::inverse(L);

  DEBUG(10, "num cg = " << topo.num_chargegroups());
  DEBUG(10, "num atoms = " << topo.num_atoms());
  DEBUG(10, "pos.size() = " << pos.size());
  
  topology::Chargegroup_Iterator cg_it = topo.chargegroup_begin(),
    cg_to = topo.chargegroup_end();

  // solute chargegroups...
  unsigned int i = 0;
  for( ; i < topo.num_solute_chargegroups(); ++cg_it, ++i){
    DEBUG(11, "cg cog " << i);
    cg_it.cog(pos, v);
    // gather on first atom...
    // v = pos(*cg_it.begin());
    v_box = v;
    put_into_box(v_box);
    trans = v_box - v;
    const math::Vec & trans_shift = math::product(cartesian_to_oblique, trans);
    
    // atoms in a chargegroup
    topology::Atom_Iterator at_it = cg_it.begin(),
      at_to = cg_it.end();
    for( ; at_it != at_to; ++at_it){
      assert(pos.size() > *at_it && shift.size() > *at_it);
      pos(*at_it) += trans;
      shift(*at_it) += trans_shift;
    } // loop over atoms
  } // loop over solute cg's

  // solvent chargegroups
  for( ; cg_it != cg_to; ++cg_it){
    // on first atom
    v = pos(**cg_it);
    v_box = v;
    put_into_box(v_box);
    trans = v_box - v;
    const math::Vec & trans_shift = math::product(cartesian_to_oblique, trans);
    
    // loop over the atoms
    topology::Atom_Iterator at_it = cg_it.begin(),
      at_to = cg_it.end();
    for( ; at_it != at_to; ++at_it){
      assert(pos.size() > *at_it && shift.size() > *at_it);
      pos(*at_it) += trans;
      shift(*at_it) += trans_shift;
    } // atoms
  } // solvent cg's

}

template<math::boundary_enum b>
void gpu::PeriodicityGpu<b>
::gather_chargegroups(configuration::Configuration & conf, 
		      topology::Topology const & topo)const
{
  math::VArray &pos = conf.current().pos;
  math::Vec v, v_box, trans;

  DEBUG(10, "num cg = " << topo.num_chargegroups());
  DEBUG(10, "num atoms = " << topo.num_atoms());
  DEBUG(10, "pos.size() = " << pos.size());
  
  topology::Chargegroup_Iterator cg_it = topo.chargegroup_begin(),
    cg_to = topo.chargegroup_end();

  for( ; cg_it != cg_to; ++cg_it){

    v_box = pos(**cg_it);
    put_into_box(v_box);

    // std::cout << "--- cg ---" << std::endl;
    
    // loop over the atoms
    topology::Atom_Iterator at_it = cg_it.begin(),
      at_to = cg_it.end();
    for( ; at_it != at_to; ++at_it){

      assert(pos.size() > *at_it);
      
      this->nearest_image(pos(*at_it), v_box, v);
      pos(*at_it) = v_box + v;
      
      // std::cout << "  " << math::v2s(v_box + v) << "\n";

    } // atoms
  } // solvent cg's
  
}

template<math::boundary_enum b>
void gpu::PeriodicityGpu<b>
::gather_molecules_into_box(configuration::Configuration & conf, 
			    topology::Topology const & topo)const
{
  math::Vec cog, o, trans;
  o = 0.0;
  
  for(size_t i=0; i<topo.molecules().size()-1; ++i){

    // first atom
    cog = conf.current().pos(topo.molecules()[i]);
      
    // put into box
    DEBUG(12, "mol " << i << " cog      = " << math::v2s(cog));
    this->nearest_image(conf.current().pos(topo.molecules()[i]), o, trans);
    conf.current().pos(topo.molecules()[i]) = trans;
    
    DEBUG(12, "mol " << i << " cog(box) = " << math::v2s(cog));
    
    // put the molecule into the box
    // using nearest image with respect to the previous atom!
    for(unsigned int a=topo.molecules()[i] + 1;
	a < topo.molecules()[i+1]; ++a){

      if (a > topo.num_atoms()){
	io::messages.add("Periodicity", "(SUB)MOLECULE information wrong", io::message::critical);
	exit(1);
      }

      this->nearest_image(conf.current().pos(a), conf.current().pos(a-1), trans);
      DEBUG(12, "atom " << a << " pos " << math::v2s(conf.current().pos(a)));
      DEBUG(12, "\tni = " << math::v2s(trans));
      
      conf.current().pos(a) = conf.current().pos(a-1) + trans;

    } // loop over atoms in molecule
  } // loop over molecules
}

template class gpu::PeriodicityGpu<math::vacuum>;
template class gpu::PeriodicityGpu<math::rectangular>;
template class gpu::PeriodicityGpu<math::triclinic>;
