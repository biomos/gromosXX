/**
 * @file periodicity.cc
 * implementation of the periodic boundary condition functions.
 */

#undef MODULE
#undef SUBMODULE
#define MODULE math
#define SUBMODULE math

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
void math::Periodicity<b>::put_into_positive_box(math::Vec &v)const {
  //this is not good in the triclinic case...  
  //Vec o(abs(this->m_box(0)), abs(this->m_box(1)), abs(this->m_box(2)));
  //o /= 2;
  Vec o = this->m_box(0) + this->m_box(1) + this->m_box(2);
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
    
    // atoms in a chargegroup
    topology::Atom_Iterator at_it = cg_it.begin(),
      at_to = cg_it.end();
    for( ; at_it != at_to; ++at_it){
      assert(pos.size() > *at_it);
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
      assert(pos.size() > *at_it);
      pos(*at_it) += trans;
    } // atoms
  } // solvent cg's

}

template<math::boundary_enum b>
void math::Periodicity<b>
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
void math::Periodicity<b>
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
void math::Periodicity<b>
::gather_molecules_into_box(configuration::Configuration & conf, 
			    topology::Topology const & topo)const
{
  Vec cog, o, trans;
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
