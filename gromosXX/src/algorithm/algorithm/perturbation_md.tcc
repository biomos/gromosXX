/**
 * @file perturbation_md.tcc
 * perturbed MD implementation.
 */

#undef MODULE
#undef SUBMODULE
#define MODULE algorithm
#define SUBMODULE algorithm

#include "../../debug.h"

template<typename t_md_spec, typename t_interaction_spec>
algorithm::Perturbation_MD<t_md_spec, t_interaction_spec>
::Perturbation_MD()
  : MD<t_md_spec, t_interaction_spec>(),
    m_d_lambda(0)
{
}

template<typename t_md_spec, typename t_interaction_spec>
inline int
algorithm::Perturbation_MD<t_md_spec, t_interaction_spec>
::pre_md(io::InInput &input)
{
  int r = MD<t_md_spec, t_interaction_spec>::pre_md(input);

  // resize the energy derivative array
  m_simulation.system().lambda_energies().
    resize(m_simulation.system().energies().bond_energy.size(),
	   m_simulation.system().energies().kinetic_energy.size());
  
  // initialize the lambda derivative fluctuations
  m_simulation.system().lambda_derivative_averages().
    resize(m_simulation.system().energies().bond_energy.size(),
	   m_simulation.system().energies().kinetic_energy.size());

  return r;
}



template<typename t_md_spec, typename t_interaction_spec>
void algorithm::Perturbation_MD<t_md_spec, t_interaction_spec>
::read_input(io::Argument &args, io::InTopology &topo,
	     io::InTrajectory &sys, io::InInput &input)
{
 
  MD<t_md_spec, t_interaction_spec>::read_input(args, topo, sys, input);

  int ntg, nlam;
  double rlam, dlamt;
  
  input.read_PERTURB(ntg, rlam, dlamt, nlam);

  if (ntg == 1){
    io::messages.add("Perturbation enabled", "Perturbation MD",
		     io::message::notice);

    m_do_perturbation = true;

    m_simulation.topology().lambda(rlam);
    m_simulation.topology().nlam(nlam);
    m_d_lambda = dlamt;
    
    if (m_d_lambda){
      io::messages.add("Perturbation: changing lambda during simulation",
		       "Perturbation MD", io::message::notice);
    }
    
  }
  
}

template<typename t_md_spec, typename t_interaction_spec>
void algorithm::Perturbation_MD<t_md_spec, t_interaction_spec>
::init_input(io::Argument &args, io::InTopology &topo,
	     io::InTrajectory &sys, io::InInput &input)
{
  // call parent
  MD<t_md_spec, t_interaction_spec>::init_input(args, topo, sys, input);

  if (m_do_perturbation){
    // initialize
    init_perturbation(args, input);
  }
  else{
    io::messages.add("using Perturbation_MD, but perturbation disabled",
		     "Perturbation_MD", io::message::notice);
  }
  
}

template<typename t_md_spec,
	 typename t_interaction_spec>
int algorithm::Perturbation_MD<t_md_spec, t_interaction_spec>
::init_perturbation(io::Argument &args, io::InInput &input)
{
  if (args.count("pert") != 1){
    // migh also be an error...
    io::messages.add("init_perturbation called but no perturbation topology",
		     "algorithm::md",
		     io::message::warning);
 
    if (m_do_perturbation) 
      throw std::string("No perturbation topology given");

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

  std::cout << "PERTURBATION\n\n"
	    << std::setw(15) << "lambda" 
	    << std::setw(20) << m_simulation.topology().lambda() << "\n"
	    << std::setw(15) << "nlam"
	    << std::setw(20) << m_simulation.topology().nlam() << "\n"
	    << std::setw(15) << "dlamt"
	    << std::setw(20) << m_d_lambda
	    << "\n\n";

  io::InPerturbationTopology pert_topo(pert_file);

  // it'd better be a perturbation topology!
  pert_topo >> m_simulation.topology();
  
  // initialize topology for lambda = ??
  m_simulation.topology().update_for_lambda();

  std::cout << "\nEND\n";
  
  return 0;

}

template<typename t_md_spec, typename t_interaction_spec>
void algorithm::Perturbation_MD<t_md_spec, t_interaction_spec>
::G96Forcefield(io::InTopology &topo,
		io::InInput &input,
		io::Argument &args)
{
  DEBUG(7, "md: create perturbed forcefield");

  Perturbed_G96_Forcefield(m_forcefield, m_simulation, topo, input, args);

  m_simulation.pressure_calculation(m_calculate_pressure);

  // decide on SHAKE
  m_distance_constraint.init(m_simulation, args, topo, input);
  
  DEBUG(7, "forcefield created");
}

template<typename t_md_spec, typename t_interaction_spec>
void algorithm::Perturbation_MD<t_md_spec, t_interaction_spec>
::print_pairlists()
{
  
  typename std::vector<typename interaction::Interaction<
    typename t_md_spec::simulation_type, t_interaction_spec> *>
    ::const_iterator it = m_forcefield.begin(),
    to = m_forcefield.end();
	
  for( ; it != to; ++it){
	  
    if ((*it)->name == "NonBonded"){
      
      (*m_print_file) << "shortrange\n" 
		      << dynamic_cast<
	typename t_interaction_spec::perturbed_nonbonded_interaction_type *>
	(*it)->pairlist()
		      << "perturbed shortrange pairlist\n"
		      << dynamic_cast<
	typename t_interaction_spec::perturbed_nonbonded_interaction_type *>
	(*it)->perturbed_pairlist()	  
		      << std::endl;
    }
    
  }
  
}

template<typename t_md_spec, typename t_interaction_spec>
void algorithm::Perturbation_MD<t_md_spec, t_interaction_spec>
::do_perturbed_energies()
{
  // and calculate the kinetic energy lamba derivative
  m_temperature.
    calculate_kinetic_energy_lambda_derivative(m_simulation);
  
  // and sum up the energy lambda derivative arrays
  m_simulation.system().lambda_energies().calculate_totals();
  
  // calculate averages
  m_simulation.system().lambda_derivative_averages().
    update(m_simulation.system().lambda_energies(), m_dt);

  if (m_print_energy && 
      (m_simulation.steps() ) % m_print_energy == 0){
    std::cout << "LAMBDA\t" << m_simulation.topology().lambda() << std::endl;
    io::print_ENERGY(std::cout, 
		     m_simulation.system().lambda_energies(),
		     m_simulation.topology().energy_groups(),
		     "dE/dLAMBDA");
  }
}



template<typename t_md_spec, typename t_interaction_spec>
void algorithm::Perturbation_MD<t_md_spec, t_interaction_spec>
::post_step()
{

  MD<t_md_spec, t_interaction_spec>::post_step();

  DEBUG(8, "md: calculate and print the perturbed energies");
  do_perturbed_energies();

  // allow for a change of lambda
  if (m_d_lambda){
    
    m_simulation.topology().lambda(m_simulation.topology().lambda() + 
				   m_d_lambda * m_dt);
    
    m_simulation.topology().update_for_lambda();
    
  }


}

template<typename t_md_spec, typename t_interaction_spec>
void algorithm::Perturbation_MD<t_md_spec, t_interaction_spec>
::post_md()
{
  MD<t_md_spec, t_interaction_spec>::post_md();

  // lambda derivatives
  if (m_do_perturbation){

    simulation::Energy energy, fluctuation;

    MD<t_md_spec, t_interaction_spec>::simulation().system().
      lambda_derivative_averages().
      average(energy, fluctuation);
  
    io::print_ENERGY(std::cout, energy,
		     MD<t_md_spec, t_interaction_spec>
		     ::simulation().topology().energy_groups(),
		     "AVERAGE ENERGY LAMBDA DERIVATIVES");

    io::print_MULTIBATH(std::cout, MD<t_md_spec, t_interaction_spec>
			::simulation().multibath(),
			energy);

    io::print_ENERGY(std::cout, fluctuation,
		     MD<t_md_spec, t_interaction_spec>
		     ::simulation().topology().energy_groups(),
		     "ENERGY LAMBDA DERIVATIVE FLUCTUATIONS");
    
    io::print_MULTIBATH(std::cout, MD<t_md_spec, t_interaction_spec>
			::simulation().multibath(),
			fluctuation);

  }
}
