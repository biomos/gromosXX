/**
 * perturbation_filter.tcc
 * filter for perturbed atoms.
 */

template<typename t_nonbonded_spec>
inline
interaction::Perturbation_Filter<t_nonbonded_spec>
::Perturbation_Filter()
  : interaction::Filter<t_nonbonded_spec>()
{
}

template<typename t_nonbonded_spec>
inline bool
interaction::Perturbation_Filter<t_nonbonded_spec>
::perturbed_atom(topology::Topology & topo,
		 configuration::Configuration & conf,
		 simulation::Simulation & sim,
		 size_t const i)
{
  if (t_nonbonded_spec::do_perturbation){
    if (topo.is_perturbed(i))
      return true;
  }
  return false;
}
