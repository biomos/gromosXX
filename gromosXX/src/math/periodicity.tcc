/**
 * @file periodicity.tcc
 * implementation of the periodic boundary condition functions.
 */

#undef MODULE
#undef SUBMODULE
#define MODULE math

#include "../debug.h"

template<math::boundary_enum b>
math::Periodicity<b>::Periodicity(boundary_enum boundary)
  : Boundary_Implementation<b>(boundary)
{
}

template<math::boundary_enum b>
void math::Periodicity<b>::put_into_box(math::Vec &v)const
{
  Vec o(0, 0, 0);
  nearest_image(v, o, v);
}

template<math::boundary_enum b>
void math::Periodicity<b>::put_into_positive_box(math::Vec &v)const
{
  Vec o(m_box(0)(0), m_box(1)(1), m_box(2)(2));
  o /= 2;
  nearest_image(v, o, v);
  v += o;  
}

template<math::boundary_enum b> template< typename t_simulation>
void math::Periodicity<b>::put_chargegroups_into_box(t_simulation &sim)const
{
  math::VArray &pos = sim.system().pos();
  math::Vec v, v_box, trans;

  simulation::chargegroup_iterator cg_it = sim.topology().chargegroup_begin(),
    cg_to = sim.topology().chargegroup_end();

  // solute chargegroups...
  size_t i = 0;
  for( ; i < sim.topology().num_solute_chargegroups(); ++cg_it, ++i){
    cg_it.cog(pos, v);
    // gather on first atom...
    // v = pos(*cg_it.begin());
    v_box = v;
    sim.system().periodicity().put_into_box(v_box);
    trans = v_box - v;
    
    // atoms in a chargegroup
    simulation::Atom_Iterator at_it = cg_it.begin(),
      at_to = cg_it.end();
    for( ; at_it != at_to; ++at_it){
      pos(*at_it) += trans;
    } // loop over atoms
  } // loop over solute cg's

  // solvent chargegroups
  for( ; cg_it != cg_to; ++cg_it){
    // on first atom
    v = pos(**cg_it);
    v_box = v;
    sim.system().periodicity().put_into_box(v_box);
    trans = v_box - v;
    
    // loop over the atoms
    simulation::Atom_Iterator at_it = cg_it.begin(),
      at_to = cg_it.end();
    for( ; at_it != at_to; ++at_it){
      pos(*at_it) += trans;
    } // atoms
  } // solvent cg's

}
