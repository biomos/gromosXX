/**
 * @file twinrange_filter.tcc
 * methods of Twinrange_Filter
 */

template<typename t_simulation, typename t_base, typename t_innerloop>
inline
interaction::Twinrange_Filter<t_simulation, t_base, t_innerloop>
::Twinrange_Filter(t_base &base)
  : Basic_Filter<t_simulation, t_base>(base),
    t_innerloop(base, *this)
{
}

template<typename t_simulation, typename t_base, typename t_innerloop>
inline void
interaction::Twinrange_Filter<t_simulation, t_base, t_innerloop>
::prepare(t_simulation &sim)
{
  m_cutoff_short_2 = sim.nonbonded().cutoff_short() * 
    sim.nonbonded().cutoff_short();
  m_cutoff_long_2 = sim.nonbonded().cutoff_long() *
    sim.nonbonded().cutoff_long();
  force().resize(sim.system().force().size());
  
  force() = 0.0;
  energies().resize(sim.system().energies().bond_energy.size());
  virial() = 0.0;
  
}

template<typename t_simulation, typename t_base, typename t_innerloop>
inline bool
interaction::Twinrange_Filter<t_simulation, t_base, t_innerloop>
::range_pair(t_simulation const &sim, size_t const i, size_t const j)
{
  math::Vec p;
  sim.system().periodicity().nearest_image(sim.system().pos()(i), 
					   sim.system().pos()(j), p);
  const double d = dot(p, p);
  
  if (d > m_cutoff_long_2) return true;
  
  if (d <= m_cutoff_short_2) return false;
  
  // ok, we are in the middle range
  // calculate the interaction...
  interaction_inner_loop(sim, i, j);
  // ...and filter
  return true;
}
