/**
 * perturbation_filter.tcc
 * filter for perturbed atoms.
 */

template<typename t_simulation, typename t_nonbonded_spec>
inline
interaction::Perturbation_Filter<t_simulation, t_nonbonded_spec>
::Perturbation_Filter()
  : interaction::Filter<t_simulation, t_nonbonded_spec>()
{
}

template<typename t_simulation, typename t_nonbonded_spec>
inline bool
interaction::Perturbation_Filter<t_simulation, t_nonbonded_spec>
::perturbed_atom(t_simulation const &sim,
		 size_t const i)
{
  if (t_nonbonded_spec::do_perturbation){
    assert(sim.topology().perturbed_atom().size() > i);
    if (sim.topology().perturbed_atom()[i])
      return true;
  }
  return false;
}
