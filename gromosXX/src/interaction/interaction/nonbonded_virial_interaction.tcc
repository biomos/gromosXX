/**
 * @file nonbonded_virial_interaction.tcc
 * template methods of Nonbonded_Virial_Interaction.
 */

#undef MODULE
#undef SUBMODULE

#define MODULE interaction
#define SUBMODULE interaction

#include "../../debug.h"

/**
 * Constructor.
 */
template<typename t_simulation, typename t_pairlist>
inline interaction::Nonbonded_Virial_Interaction<t_simulation, t_pairlist>
::Nonbonded_Virial_Interaction(virial_enum virial)
  : Nonbonded_Interaction<t_simulation, t_pairlist>::Nonbonded_Interaction(),
    m_virial_type(virial)
{
}

/**
 * Destructor.
 */
template<typename t_simulation, typename t_pairlist>
inline interaction::Nonbonded_Virial_Interaction<t_simulation, t_pairlist>
::~Nonbonded_Virial_Interaction()
{
  DEBUG(4, "Nonbonded_Virial_Interaction::destructor");
}

/**
 * calculate nonbonded forces and energies.
 */
template<typename t_simulation, typename t_pairlist>
inline void interaction::Nonbonded_Virial_Interaction<t_simulation, t_pairlist>
::calculate_interactions(t_simulation &sim)
{
  DEBUG(7, "nonbonded virial interactions");
  
  // prepare for the virial calculation
  m_com_pos.resize(sim.system().pos().size());
  
  simulation::Molecule_Iterator m_it = sim.topology().molecule_begin(),
    m_to = sim.topology().molecule_end();
  
  math::Vec com_pos;
  math::Matrix com_ekin;
  
  sim.system().molecular_kinetic_energy() = 0.0;

  for( ; m_it != m_to; ++m_it){
    m_it.com(sim.system(), sim.topology().mass(), com_pos, com_ekin);

    for(int i=0; i<3; ++i)
      for(int j=0; j<3; ++j)
	sim.system().molecular_kinetic_energy()(i,j) += com_ekin(i,j);
    
    simulation::Atom_Iterator a_it = m_it.begin(),
      a_to = m_it.end();
    
    math::VArray &pos = sim.system().pos();

    for( ; a_it != a_to; ++a_it){
      assert(unsigned(m_com_pos.size()) > *a_it);
      // m_com_pos(*a_it) = pos(*a_it) - com_pos;
      sim.system().periodicity().nearest_image(pos(*a_it), com_pos, m_com_pos(*a_it));
    }

  }

  // reset the virial
  if(!(sim.steps() % sim.nonbonded().update())){
    m_longrange_virial = 0.0;
  }
  sim.system().virial() = 0.0;

  // do the work...
  Nonbonded_Interaction<t_simulation, t_pairlist>::calculate_interactions(sim);

  for(int i=0; i<3; ++i)
    for(int j=0; j<3; ++j)
      sim.system().virial()(i,j) = -0.5*(sim.system().virial()(i,j)+
			       m_longrange_virial(i,j));

}

/**
 * helper function to calculate forces and energies, 
 * stores them in the arrays pointed to by parameters
 * to make it usable for longrange calculations.
 */
template<typename t_simulation, typename t_pairlist>
inline void interaction::Nonbonded_Virial_Interaction<t_simulation, t_pairlist>
::do_interactions(t_simulation &sim, typename t_pairlist::iterator it, 
		  typename t_pairlist::iterator to,
		  typename Nonbonded_Interaction<t_simulation, t_pairlist>
		  ::nonbonded_type_enum range)
{
  math::Vec r, f, r_com;
  double e_lj, e_crf;
  
  math::VArray &pos = sim.system().pos();
  math::SArray &charge = sim.topology().charge();

  math::VArray *force;
  if (range == shortrange) force = &sim.system().force();
  else force = &m_longrange_force;
  
  math::Matrix *virial;
  if (range == shortrange) virial = &sim.system().virial();
  else virial = &m_longrange_virial;

  DEBUG(7, "\tcalculate interactions with virial");  

  for( ; it != to; ++it){
    
    DEBUG(10, "\tpair\t" << it.i() << "\t" << *it);

    assert(pos.size() > it.i() && pos.size() > *it);
    assert(m_com_pos.size() > it.i() && m_com_pos.size() > *it);
    
    sim.system().periodicity().nearest_image(pos(it.i()), pos(*it), r);
    // sim.system().periodicity().nearest_image(m_com_pos(it.i()), m_com_pos(*it), r_com);

    const lj_parameter_struct &lj = 
      lj_parameter(sim.topology().iac(it.i()),
		   sim.topology().iac(*it));

    DEBUG(11, "\tlj-parameter c6=" << lj.c6 << " c12=" << lj.c12);

    assert(charge.size() > it.i() && charge.size() > *it);
    
    lj_crf_interaction(r, lj.c6, lj.c12,
		       charge(it.i()) * 
		       charge(*it),
		       f, e_lj, e_crf);

    (*force)(it.i()) += f;
    (*force)(*it) -= f;

    for(int i=0; i<3; ++i)
      for(int j=0; j<3; ++j)
	(*virial)(i, j) += (r(i) - m_com_pos(it.i())(i) + m_com_pos(*it)(i)) 
	  * f(j);

    // energy
    sim.system().energies().lj_energy[sim.topology().atom_energy_group(it.i())]
      [sim.topology().atom_energy_group(*it)] += e_lj;

    sim.system().energies().crf_energy[sim.topology().atom_energy_group(it.i())]
      [sim.topology().atom_energy_group(*it)] += e_crf;

  }
    
}

/**
 * helper function to calculate 1,4 forces and energies.
 */
template<typename t_simulation, typename t_pairlist>
inline void interaction::Nonbonded_Virial_Interaction<t_simulation, t_pairlist>
::do_14_interactions(t_simulation &sim)
{
  
  Nonbonded_Interaction<t_simulation, t_pairlist>::do_14_interactions(sim);
  
}


