/**
 * @file periodicity.cc
 * implementation of the periodic boundary condition functions.
 */

#undef MODULE
#undef SUBMODULE
#define MODULE math


template<math::boundary_enum b>
math::Periodicity<b>::Periodicity(math::Box const & bb) 
  : Boundary_Implementation<b>(bb)
{
}

template<math::boundary_enum b>
void math::Periodicity<b>::put_into_box(math::Vec &v)const
{
  Vec o(0, 0, 0);
  this->nearest_image(v, o, v);
}

template<math::boundary_enum b>
void math::Periodicity<b>::put_into_positive_box(math::Vec &v)const
{
  Vec o(this->m_box(0)(0), this->m_box(1)(1), this->m_box(2)(2));
  o /= 2;
  this->nearest_image(v, o, v);
  v += o;  
}

template<math::boundary_enum b>
void math::Periodicity<b>
::put_chargegroups_into_box(configuration::Configuration & conf, 
			    topology::Topology const & topo)const
{
  math::VArray &pos = conf.current().pos;
  math::Vec v, v_box, trans;

  topology::Chargegroup_Iterator cg_it = topo.chargegroup_begin(),
    cg_to = topo.chargegroup_end();

  // solute chargegroups...
  size_t i = 0;
  for( ; i < topo.num_solute_chargegroups(); ++cg_it, ++i){
    cg_it.cog(pos, v);
    // gather on first atom...
    // v = pos(*cg_it.begin());
    v_box = v;
    put_into_box(v_box);
    trans = v_box - v;
    
    // atoms in a chargegroup
    topology::Atom_Iterator at_it = cg_it.begin(),
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
    put_into_box(v_box);
    trans = v_box - v;
    
    // loop over the atoms
    topology::Atom_Iterator at_it = cg_it.begin(),
      at_to = cg_it.end();
    for( ; at_it != at_to; ++at_it){
      pos(*at_it) += trans;
    } // atoms
  } // solvent cg's

}

template<math::boundary_enum b>
void math::Periodicity<b>
::gather_molecules_into_box(configuration::Configuration & conf, 
			    topology::Topology & topo)const
{
  Vec cog, o, trans;
  o = 0.0;
  
  // std::cout.precision(12);
  // std::cout.setf(std::ios_base::fixed, std::ios_base::floatfield);

  for(size_t i=0; i<topo.molecules().size()-1; ++i){

    // first atom
    cog = conf.current().pos(topo.molecules()[i]);
      
    // put into box
    this->nearest_image(cog, o, trans);
    // cog = o + trans;
    cog = trans;
    
    // put the molecule into the box
    for(size_t a=topo.molecules()[i];
	a<topo.molecules()[i+1]; ++a){

      nearest_image(conf.current().pos(a), cog, trans);
      conf.current().pos(a) = cog + trans;

    } // loop over atoms in molecule
  
  
  } // loop over molecules
  
}

