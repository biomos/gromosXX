/**
 * @file perturbation_filter.tcc
 * template methods of the Perturbation_Filter
 */

template<typename t_simulation, typename t_base, bool do_perturbation>
inline interaction::Perturbation_Filter<t_simulation, t_base, do_perturbation>
::Perturbation_Filter(t_base &base)
  : Basic_Filter<t_simulation, t_base>(base)
{
}

template<typename t_simulation, typename t_base, bool do_perturbation>
inline void 
interaction::Perturbation_Filter<t_simulation, t_base, do_perturbation>
::prepare(t_simulation &sim)
{
  Basic_Filter<t_simulation, t_base>::prepare(sim);
}

template<typename t_simulation, typename t_base, bool do_perturbation>
inline bool
interaction::Perturbation_Filter<t_simulation, t_base, do_perturbation>
::perturbed_pair(t_simulation const &sim, size_t const i, size_t const j)
{
  if (do_perturbation){
    if (sim.topology().perturbed_atom()[i] ||
	sim.topology().perturbed_atom()[j])
      return true;
    return false;
  }
  return false;
}

