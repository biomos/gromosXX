/**
 * @file range_filter.cc
 * methods of Range_Filter
 */

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE filter

#include <util/debug.h>

template<typename t_interaction_spec, typename t_perturbation_spec>
inline
interaction::Range_Filter<t_interaction_spec, t_perturbation_spec>
::Range_Filter()
  : interaction::Filter()
{
}

template<typename t_interaction_spec, typename t_perturbation_spec>
inline void
interaction::Range_Filter<t_interaction_spec, t_perturbation_spec>
::set_cutoff(double const cutoff_short, double const cutoff_long)
{
  m_cutoff_long = cutoff_long;
  m_cutoff_short = cutoff_short;
  m_cutoff_short_2 = cutoff_short * cutoff_short;
  m_cutoff_long_2  = cutoff_long * cutoff_long;
}

template<typename t_interaction_spec, typename t_perturbation_spec>
inline void
interaction::Range_Filter<t_interaction_spec, t_perturbation_spec>
::prepare_cog(topology::Topology & topo,
	      configuration::Configuration & conf,
	      simulation::Simulation & sim)
{
  DEBUG(11, "Range_Filter::prepare_cog");
  
  // first put the chargegroups into the box
  math::Periodicity<t_interaction_spec::boundary_type> 
    periodicity(conf.current().box);
  
  periodicity.put_chargegroups_into_box(conf, topo);

  // calculate cg cog's
  m_cg_cog.resize(topo.num_chargegroups());
  math::VArray const &pos = conf.current().pos;

  // calculate all center of geometries
  topology::Chargegroup_Iterator
    cg1 =   topo.chargegroup_begin(),
    cg_to = topo.chargegroup_end();

  unsigned int i,
		num_cg = topo.num_solute_chargegroups();

  // solute
  for(i=0; i < num_cg; ++cg1, ++i){
    cg1.cog(pos, m_cg_cog(i));
  }
  // solvent
  for( ; cg1 != cg_to; ++cg1, ++i){
    m_cg_cog(i) = pos(**cg1);
  }  
}

template<typename t_interaction_spec, typename t_perturbation_spec>
inline void
interaction::Range_Filter<t_interaction_spec, t_perturbation_spec>
::grid_cog(topology::Topology & topo,
	   configuration::Configuration & conf,
	   simulation::Simulation & sim,
	   Chargegroup_Grid_type & grid)
{
  math::Periodicity<t_interaction_spec::boundary_type>
    periodicity(conf.current().box);

  DEBUG(8, "box:\n\t" << periodicity.box()(0)(0)
	<< "\t" << periodicity.box()(0)(1)
	<< "\t" << periodicity.box()(0)(2)
	<< "\n\t" << periodicity.box()(1)(0)
	<< "\t" << periodicity.box()(1)(1)
	<< "\t" << periodicity.box()(1)(2)
	<< "\n\t" << periodicity.box()(2)(0)
	<< "\t" << periodicity.box()(2)(1)
	<< "\t" << periodicity.box()(2)(2));
	
  for(unsigned int i = 0, i_to = unsigned(m_cg_cog.size()); i != i_to; ++i)
    grid.add(m_cg_cog(i), i);
  
}

template<typename t_interaction_spec, typename t_perturbation_spec>
inline bool
interaction::Range_Filter<t_interaction_spec, t_perturbation_spec>
::range_chargegroup_pair(topology::Topology & topo,
			 configuration::Configuration & conf,
			 simulation::Simulation & sim,
			 Nonbonded_Set<t_interaction_spec, t_perturbation_spec> & nbs,
			 unsigned int i, unsigned int j,
			 topology::Chargegroup_Iterator const & it_i,
			 topology::Chargegroup_Iterator const & it_j,
			 Periodicity_type const & periodicity)
{
  DEBUG(10, "Range_Filter::range_chargegroup_pair " << i << " - " << j);
  
  math::Vec p;
  
  assert(unsigned(m_cg_cog.size()) > i &&
	 unsigned(m_cg_cog.size()) > j);
  
  periodicity.nearest_image(m_cg_cog(i), m_cg_cog(j), p);
 
  // the distance
  const double d = dot(p, p);

  DEBUG(11, "\tdistance: " << d);

  if (d > m_cutoff_long_2){        // OUTSIDE: filter
    DEBUG(11, "cg pair " << i << " - " << j << " outside range");
    return true;
  }
  
  if (d < m_cutoff_short_2){       // SHORTRANGE: no filter
    DEBUG(11, "cg pair " << i << " - " << j << " short range");
    return false;
  }

  DEBUG(11, "cg pair " << i << " - " << j << " long range");  
    
  // LONGRANGE: interactions and filter

  topology::Atom_Iterator a1 = it_i.begin(),
    a1_to = it_i.end();
  
  for( ; a1 != a1_to; ++a1){
    for(topology::Atom_Iterator
	  a2 = it_j.begin(),
	  a2_to = it_j.end();
	a2 != a2_to; ++a2){

      // the interactions
      nbs.add_longrange_pair
	(topo, conf, sim, *a1, *a2, periodicity);

    } // loop over atom 2 of cg1
  } // loop over atom 1 of cg1

  return true;
}

template<typename t_interaction_spec, typename t_perturbation_spec>
inline bool
interaction::Range_Filter<t_interaction_spec, t_perturbation_spec>
::range_chargegroup_pair(topology::Topology & topo,
			 configuration::Configuration & conf,
			 simulation::Simulation & sim,
			 Nonbonded_Set<t_interaction_spec, t_perturbation_spec> & nbs,
			 unsigned int i, unsigned int j,
			 topology::Chargegroup_Iterator const & it_i,
			 topology::Chargegroup_Iterator const & it_j,
			 int pc,
			 Periodicity_type const & periodicity)
{
  DEBUG(11, "Range_Filter::range_chargegroup_pair (shift) " << i << " - " << j);
  
  math::Vec p;
  
  assert(unsigned(m_cg_cog.size()) > i &&
	 unsigned(m_cg_cog.size()) > j);
  
  p = m_cg_cog(i) + periodicity.shift(pc).pos - m_cg_cog(j);
  
  // the distance
  const double d = dot(p, p);
  DEBUG(11, "\tdistance: " << d);

  if (d > m_cutoff_long_2){        // OUTSIDE: filter
    DEBUG(11, "cg pair " << i << " - " << j << " outside range");
    return true;
  }
  
  if (d <= m_cutoff_short_2){       // SHORTRANGE: no filter
    DEBUG(11, "cg pair " << i << " - " << j << " short range");
    return false;
  }

  DEBUG(11, "cg pair " << i << " - " << j << " long range");  
  // LONGRANGE: interactions and filter

  topology::Atom_Iterator a1 = it_i.begin(),
    a1_to = it_i.end();
  
  for( ; a1 != a1_to; ++a1){
    for(topology::Atom_Iterator
	  a2 = it_j.begin(),
	  a2_to = it_j.end();
	a2 != a2_to; ++a2){

      // the interactions
      nbs.add_longrange_pair
	(topo, conf, sim, *a1, *a2, periodicity, pc);

    } // loop over atom 2 of cg1
  } // loop over atom 1 of cg1

  return true;
}

template<typename t_interaction_spec, typename t_perturbation_spec>
inline bool
interaction::Range_Filter<t_interaction_spec, t_perturbation_spec>
::range_atom_pair(topology::Topology & topo,
		  configuration::Configuration & conf,
		  simulation::Simulation & sim,
		  Nonbonded_Set<t_interaction_spec, t_perturbation_spec> & nbs,
		  unsigned int i, unsigned int j,
		  Periodicity_type const & periodicity)
{

  DEBUG(11, "Range_Filter::range_atom_pair " << i << " - " << j);
  
  math::Vec p;
  
  periodicity.nearest_image(conf.current().pos(i),
			      conf.current().pos(j), p);
 
  // the distance
  const double d = dot(p, p);

  if (d > m_cutoff_long_2){        // OUTSIDE: filter
    DEBUG(11, "atom pair " << i << " - " << j << " outside range");
    return true;
  }
  
  if (d < m_cutoff_short_2){       // SHORTRANGE: no filter
    DEBUG(11, "atom pair " << i << " - " << j << " short range");
    return false;
  }

  DEBUG(11, "atom pair " << i << " - " << j << " long range");  
    
  // LONGRANGE: interactions and filter

  // the interactions
  nbs.add_longrange_pair(topo, conf, sim, i, j, periodicity);

  return true;
}

template<typename t_interaction_spec, typename t_perturbation_spec>
inline bool
interaction::Range_Filter<t_interaction_spec, t_perturbation_spec>
::range_atom_pair(topology::Topology & topo,
		  configuration::Configuration & conf,
		  simulation::Simulation & sim,
		  Nonbonded_Set<t_interaction_spec, t_perturbation_spec> & nbs,
		  unsigned int i, unsigned int j,
		  int pc,
		  Periodicity_type const & periodicity)
{

  DEBUG(11, "Range_Filter::range_atom_pair " << i << " - " << j);
  
  math::Vec p;
  
  p = conf.current().pos(i) + periodicity.shift(pc).pos - 
    conf.current().pos(j);
 
  // the distance
  const double d = dot(p, p);

  if (d > m_cutoff_long_2){        // OUTSIDE: filter
    DEBUG(11, "atom pair " << i << " - " << j << " outside range");
    return true;
  }
  
  if (d < m_cutoff_short_2){       // SHORTRANGE: no filter
    DEBUG(11, "atom pair " << i << " - " << j << " short range");
    return false;
  }

  DEBUG(11, "atom pair " << i << " - " << j << " long range");  
    
  // LONGRANGE: interactions and filter

  // the interactions
  nbs.add_longrange_pair(topo, conf, sim, i, j,
			 periodicity, pc);

  return true;
}