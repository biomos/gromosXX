/**
 * @file perturbation_md.tcc
 * perturbed MD implementation.
 */
#undef MODULE
#undef SUBMODULE
#define MODULE algorithm
#define SUBMODULE constraint

#include "../../debug.h"

template<typename t_simulation,
	 typename t_temperature,
	 typename t_pressure,
	 typename t_distance_constraint,
	 typename t_integration>
algorithm::Perturbation_MD<t_simulation, t_temperature, t_pressure,
			   t_distance_constraint, t_integration>
::Perturbation_MD(t_simulation &sim)
  : MD<t_simulation, t_temperature, t_pressure,
       t_distance_constraint, t_integration>(sim)
{
}

template<typename t_simulation,
	 typename t_temperature,
	 typename t_pressure,
	 typename t_distance_constraint,
	 typename t_integration>
void algorithm::Perturbation_MD<t_simulation, t_temperature, t_pressure, 
				t_distance_constraint, t_integration>
::init_input(io::Argument &args, io::InTopology &topo,
	     io::InTrajectory &sys, io::InInput &input)
{
  // nothing to do...
  MD<t_simulation, t_temperature, t_pressure,
    t_distance_constraint, t_integration>::init_input(args, topo, sys, input);
}

/**
 * read the input and setup a standard simulation.
 */
template<typename t_simulation,
	 typename t_temperature,
	 typename t_pressure,
	 typename t_distance_constraint,
	 typename t_integration>
void algorithm::Perturbation_MD<t_simulation, t_temperature,
			       t_pressure, t_distance_constraint,
			       t_integration>
::read_input(io::Argument &args, io::InTopology &topo,
	     io::InTrajectory &sys, io::InInput &input)
{
  // std::cerr << "read input" << std::endl;
  MD<t_simulation, t_temperature, t_pressure,
    t_distance_constraint, t_integration>::read_input(args, topo, sys, input);
  
  int ntg;
  double rlam, dlamt;
  input.read_PERTURB(ntg, rlam, dlamt);
  
  if (ntg == 1){
    io::messages.add("Perturbation enabled", "Perturbation_MD",
		     io::message::notice);

    m_simulation.topology().lambda(rlam);
    // initialize
    init_perturbation(args, input);
  }
  else{
    io::messages.add("using Perturbation_MD, but perturbation disabled",
		     "Perturbation_MD", io::message::notice);
  }
  
}

template<typename t_simulation,
	 typename t_temperature,
	 typename t_pressure,
	 typename t_distance_constraint,
	 typename t_integration>
int algorithm::Perturbation_MD<t_simulation, t_temperature,
			       t_pressure, t_distance_constraint,
			       t_integration>
::init_perturbation(io::Argument &args, io::InInput &input)
{
  if (args.count("pert") != 1){
    // migh also be an error...
    io::messages.add("init_perturbation called but no perturbation topology",
		     "algorithm::md",
		     io::message::warning);
    return 0;
  }

  std::ifstream pert_file(args["pert"].c_str());
  if (!pert_file.good()){
    io::messages.add("unable to open perturbation topology: "
		     + args["pert"],
		     "algorithm::md",
		     io::message::error);
    return 1;
  }
  else
    io::messages.add("parsing perturbation topology file: "
		     + args["pert"], "algorithm::md",
		     io::message::notice);

  io::InPerturbationTopology pert_topo(pert_file);

  // it'd better be a perturbation topology!
  pert_topo >> m_simulation.topology();

  // resize the energy array
  m_simulation.system().lambda_energies().
    resize(m_simulation.system().energies().bond_energy.size());
  m_simulation.system().lambda_energies().
    kinetic_energy.resize(m_simulation.system().energies().kinetic_energy.size());
    
  return 0;

}

template<typename t_simulation,
	 typename t_temperature,
	 typename t_pressure,
	 typename t_distance_constraint,
	 typename t_integration>
void algorithm::Perturbation_MD<t_simulation, t_temperature, t_pressure, 
				t_distance_constraint, t_integration>
::G96Forcefield(io::InTopology &topo,
		io::InInput &input)
{
  // call parent
  MD<t_simulation, t_temperature, t_pressure,
    t_distance_constraint, t_integration>::G96Forcefield(topo, input);
  
  // add the perturbed interactions
  // add the perturbed interactions to the forcefield
  assert(m_qbond_interaction != NULL);

  interaction::Perturbed_Quartic_Bond_Interaction<t_simulation>
    *the_perturbed_qbond_interaction = 
    new interaction::Perturbed_Quartic_Bond_Interaction<t_simulation>
    (*m_qbond_interaction);

  m_forcefield.push_back(the_perturbed_qbond_interaction);
}


template<typename t_simulation,
	 typename t_temperature,
	 typename t_pressure,
	 typename t_distance_constraint,
	 typename t_integration>
void algorithm::Perturbation_MD<t_simulation, t_temperature, t_pressure, 
		   t_distance_constraint, t_integration>
::do_energies()
{
  // do the normal energies
  MD<t_simulation, t_temperature, t_pressure,
    t_distance_constraint, t_integration>::do_energies();

  // and sum up the energy lambda derivative arrays
  m_simulation.system().lambda_energies().calculate_totals();
  
  if (m_print_energy && m_simulation.steps() % m_print_energy == 0){
    io::print_ENERGY(std::cout, m_simulation.system().lambda_energies(),
		     m_simulation.topology().energy_groups(), "dE/dLAMBDA");
  }
}
