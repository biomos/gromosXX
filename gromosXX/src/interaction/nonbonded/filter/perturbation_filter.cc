/**
 * perturbation_filter.cc
 * filter for perturbed atoms.
 */

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE filter

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
		 unsigned int i)
{
  if (t_nonbonded_spec::do_perturbation){
    if (i < topo.num_solute_atoms() &&
	topo.is_perturbed(i))
      return true;
  }
  return false;
}
