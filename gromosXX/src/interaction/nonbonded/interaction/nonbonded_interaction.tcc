/**
 * @file nonbonded_interaction.tcc
 * template methods of Nonbonded_Interaction.
 */

#undef MODULE
#undef SUBMODULE

#define MODULE interaction
#define SUBMODULE interaction

#include <util/debug.h>

/**
 * Constructor.
 */
template<typename t_interaction_spec>
inline
interaction::Nonbonded_Interaction<t_interaction_spec>
::Nonbonded_Interaction()
  : Interaction("NonBonded"),
    Nonbonded_Base(),
    Storage(),
    t_interaction_spec::nonbonded_innerloop_type(*dynamic_cast<Nonbonded_Base *>(this)),
    m_pairlist_algorithm()
{
}

/**
 * Destructor.
 */
template<typename t_interaction_spec>
inline 
interaction::Nonbonded_Interaction<t_interaction_spec>
::~Nonbonded_Interaction()
{
  DEBUG(4, "Nonbonded_Interaction::destructor");
}

/**
 * calculate nonbonded forces and energies.
 */
template<typename t_interaction_spec>
inline int 
interaction::Nonbonded_Interaction<t_interaction_spec>
::calculate_interactions(topology::Topology & topo,
			 configuration::Configuration & conf,
			 simulation::Simulation & sim)
{
  DEBUG(4, "Nonbonded_Interaction::calculate_interactions");

  // initialize the constants
  if (!sim.steps())
    initialize(topo, conf, sim);

  // need to update pairlist?
  if(!(sim.steps() % sim.param().pairlist.skip_step)){
    // create a pairlist
    // zero the longrange forces, energies, virial
    force = 0.0;
    energies.zero();
    DEBUG(15, "zero the longrange lambda energies");
    perturbed_energy_derivatives.zero();
    virial_tensor = 0.0;
    
    DEBUG(7, "\tupdate the parlist");
    m_pairlist_algorithm.update(topo, conf, sim, *this);
    DEBUG(7, "\tpairlist updated");
  }

  // calculate forces / energies
  DEBUG(7, "\tshort range");

  do_interactions(topo, conf, sim, 
		  m_pairlist.begin(),
		  m_pairlist.end() );
  
  // add long-range force
  DEBUG(7, "\tadd long range forces and energies");

  conf.current().force += force;
  
  // and long-range energies
  for(size_t i = 0; i < energies.lj_energy.size(); ++i){
    for(size_t j = 0; j < energies.lj_energy.size(); ++j){
      conf.current().energies.lj_energy[i][j] += 
	energies.lj_energy[i][j];
      conf.current().energies.crf_energy[i][j] += 
	energies.crf_energy[i][j];
    }
  }
  
  // add longrange virial
  if (t_interaction_spec::do_virial){
    DEBUG(7, "\tadd long range virial");
    for(size_t i=0; i<3; ++i)
      for(size_t j=0; j<3; ++j)
	conf.current().virial_tensor(i,j) =
	  conf.current().virial_tensor(i,j) + virial_tensor(i,j);
  }
  
  // add 1,4 - interactions
  DEBUG(7, "\t1,4 - interactions");
  do_14_interactions(topo, conf, sim);

  // possibly do the RF contributions due to excluded atoms
  if(sim.param().longrange.rf_excluded){
    DEBUG(7, "\tRF excluded interactions and self term");
    do_RF_excluded_interactions(topo, conf, sim);
  }

  return 0;
  
}

/**
 * add a shortrange interaction
 */
template<typename t_interaction_spec>
inline void
interaction::Nonbonded_Interaction<t_interaction_spec>
::add_shortrange_pair(topology::Topology & topo,
		      configuration::Configuration & conf,
		      simulation::Simulation & sim, 
		      size_t const i, size_t const j)
{
  assert(pairlist().size() > i);
  pairlist()[i].push_back(j);
}

/**
 * add a shortrange interaction
 */
template<typename t_interaction_spec>
inline void
interaction::Nonbonded_Interaction<t_interaction_spec>
::add_shortrange_pair(topology::Topology & topo,
		      configuration::Configuration & conf,
		      simulation::Simulation & sim,
		      size_t const i, size_t const j,
		      int pc)
{
  DEBUG(10, "\tadding shortrange pair i=" << i << " j=" << j << " pc=" << pc);

  assert(pairlist().size() > i);
  assert(pc >=0 && pc <= 26);
  
  pairlist()[i].push_back((pc << 26) + j);
}

/**
 * add a longrange interaction
 */
template<typename t_interaction_spec>
inline void
interaction::Nonbonded_Interaction<t_interaction_spec>
::add_longrange_pair(topology::Topology & topo,
		     configuration::Configuration & conf,
		     simulation::Simulation & sim, 
		     size_t const i, size_t const j,
		     math::Periodicity<t_interaction_spec::boundary_type>
		     const & periodicity)
{
  interaction_innerloop(topo, conf, i, j, *this, periodicity);
}

/**
 * add a longrange interaction
 */
template<typename t_interaction_spec>
inline void
interaction::Nonbonded_Interaction<t_interaction_spec>
::add_longrange_pair(topology::Topology & topo,
		     configuration::Configuration & conf,
		     simulation::Simulation & sim,
		     size_t const i, size_t const j,
		     math::Periodicity<t_interaction_spec::boundary_type>
		     const & periodicity, int pc)
{
  DEBUG(10, "\tadding long range pair i=" << i << " j=" << j << " pc=" << pc);

  assert(pc >=0 && pc <= 26);

  interaction_innerloop(topo, conf, i, j, *this, periodicity, pc);
}

//==================================================
// interaction loops
//==================================================

/**
 * helper function to calculate forces and energies, 
 * stores them in the arrays pointed to by parameters
 * to make it usable for longrange calculations.
 */
template<typename t_interaction_spec>
inline void interaction::Nonbonded_Interaction<t_interaction_spec>
::do_interactions(topology::Topology & topo,
		  configuration::Configuration & conf,
		  simulation::Simulation & sim, 
		  Pairlist::iterator it, 
		  Pairlist::iterator to)
{  
  DEBUG(7, "\tcalculate interactions");  

  math::Periodicity<t_interaction_spec::boundary_type> 
    periodicity(conf.current().box);

  if (t_interaction_spec::do_bekker){

    periodicity.recalc_shift_vectors();

    int pc;
    size_t j;
    // translate the atom j
    DEBUG(9, "nonbonded_interaction: grid based pairlist");

    for( ; it != to; ++it){
      
      pc = (*it >> 26);
      j = (*it & 67108863);
      
      DEBUG(10, "\tnonbonded_interaction: i " << it.i() << " j " << j
	    << " pc " << pc);
      
      interaction_innerloop(topo, conf, it.i(), j, 
			    conf.current(), periodicity, pc);
    }
  }
  else{ // no grid based pairlist
    DEBUG(9, "nonbonded_interaction: no grid based pairlist");
    
    for( ; it != to; ++it){    

      DEBUG(10, "\tnonbonded_interaction: i " << it.i() << " j " << *it);

      // shortrange, therefore store in simulation.system()
      interaction_innerloop(topo, conf, it.i(), *it, 
			    conf.current(), periodicity);
      
    }
    
  }
  
}

/**
 * helper function to calculate the forces and energies from the
 * 1,4 interactions.
 */
template<typename t_interaction_spec>
inline void interaction::Nonbonded_Interaction<t_interaction_spec>
::do_14_interactions(topology::Topology & topo,
		     configuration::Configuration & conf,
		     simulation::Simulation & sim)
{
  DEBUG(7, "\tcalculate 1,4-interactions");

  math::Periodicity<t_interaction_spec::boundary_type> 
    periodicity(conf.current().box);

  std::set<int>::const_iterator it, to;
  
  for(size_t i=0; i<topo.num_solute_atoms(); ++i){
    it = topo.one_four_pair(i).begin();
    to = topo.one_four_pair(i).end();
    
    for( ; it != to; ++it){

      one_four_interaction_innerloop(topo, conf, i, *it, periodicity);

    } // loop over 1,4 pairs
  } // loop over solute atoms
}  

/**
 * helper function to calculate the forces and energies from the
 * RF contribution of excluded atoms and self term
 */
template<typename t_interaction_spec>
inline void interaction::Nonbonded_Interaction<t_interaction_spec>
::do_RF_excluded_interactions(topology::Topology & topo,
			      configuration::Configuration & conf,
			      simulation::Simulation & sim)
{
  
  DEBUG(7, "\tcalculate RF excluded interactions");

  math::Periodicity<t_interaction_spec::boundary_type> 
    periodicity(conf.current().box);
  
  for(size_t i=0; i<topo.num_solute_atoms(); ++i){
    
    RF_excluded_interaction_innerloop(topo, conf, i, periodicity);
    
  } // loop over solute atoms

  // Solvent
  topology::Chargegroup_Iterator cg_it = topo.chargegroup_begin(),
    cg_to = topo.chargegroup_end();
  cg_it += topo.num_solute_chargegroups();
  
  for( ; cg_it != cg_to; ++cg_it){

    RF_solvent_interaction_innerloop(topo, conf, cg_it, periodicity);

  } // loop over solvent charge groups
}  

/**
 * initialize the arrays
 */
template<typename t_interaction_spec>
inline void interaction::Nonbonded_Interaction<t_interaction_spec>
::initialize(topology::Topology & topo,
	     configuration::Configuration & conf,
	     simulation::Simulation & sim)
{
  DEBUG(15, "nonbonded_interaction::initialize");
  
  Nonbonded_Base::initialize(sim);

  force.resize(conf.current().force.size());

  energies.resize(conf.current().energies.bond_energy.size(),
		    conf.current().energies.kinetic_energy.size());

  perturbed_energy_derivatives.resize
    (conf.current().perturbed_energy_derivatives.bond_energy.size(),
     conf.current().perturbed_energy_derivatives.kinetic_energy.size());
  
}
