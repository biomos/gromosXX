/**
 * @file perturbation_md.tcc
 * perturbed MD implementation.
 */
#undef MODULE
#undef SUBMODULE
#define MODULE algorithm
#define SUBMODULE algorithm

#include "../../debug.h"

template<typename t_spec>
algorithm::Perturbation_MD<t_spec>
::Perturbation_MD()
  : MD<t_spec>()
{
}

template<typename t_spec>
void algorithm::Perturbation_MD<t_spec>
::init_input(io::Argument &args, io::InTopology &topo,
	     io::InTrajectory &sys, io::InInput &input)
{
  // nothing to do...
  MD<t_spec>::init_input(args, topo, sys, input);
}

/**
 * read the input and setup a standard simulation.
 */
template<typename t_spec>
void algorithm::Perturbation_MD<t_spec>
::read_input(io::Argument &args, io::InTopology &topo,
	     io::InTrajectory &sys, io::InInput &input)
{
  // std::cerr << "read input" << std::endl;
  MD<t_spec>::read_input(args, topo, sys, input);
  
  int ntg, nlam;
  double rlam, dlamt;
  
  input.read_PERTURB(ntg, rlam, dlamt, nlam);

  if (ntg == 1){
    io::messages.add("Perturbation enabled", "Perturbation_MD",
		     io::message::notice);

    m_do_perturbation = true;

    m_simulation.topology().lambda(rlam);
    m_simulation.topology().nlam(nlam);
    
    // initialize
    init_perturbation(args, input);
  }
  else{
    io::messages.add("using Perturbation_MD, but perturbation disabled",
		     "Perturbation_MD", io::message::notice);
    // resize the energy derivative array
    m_simulation.system().lambda_energies().
      resize(m_simulation.system().energies().bond_energy.size(),
	   m_simulation.system().energies().kinetic_energy.size());
    
    // initialize the lambda derivative fluctuations
    m_simulation.system().lambda_derivative_averages().
      resize(m_simulation.system().energies().bond_energy.size(),
	     m_simulation.system().energies().kinetic_energy.size());
  }
  
}

template<typename t_spec>
int algorithm::Perturbation_MD<t_spec>
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

  std::cout << "PERTURBATION\n\n"
	    << std::setw(15) << "lambda" 
	    << std::setw(20) << m_simulation.topology().lambda() << "\n"
	    << std::setw(15) << "nlam"
	    << std::setw(20) << m_simulation.topology().nlam()
	    << "\n\n";

  io::InPerturbationTopology pert_topo(pert_file);

  // it'd better be a perturbation topology!
  pert_topo >> m_simulation.topology();

  // resize the energy derivative array
  m_simulation.system().lambda_energies().
    resize(m_simulation.system().energies().bond_energy.size(),
	   m_simulation.system().energies().kinetic_energy.size());

    
  // initialize the lambda derivative fluctuations
  m_simulation.system().lambda_derivative_averages().
    resize(m_simulation.system().energies().bond_energy.size(),
	   m_simulation.system().energies().kinetic_energy.size());
  
  // initialize topology for lambda = ??
  m_simulation.topology().update_for_lambda();

  std::cout << "\nEND\n";
  
  return 0;

}

template<typename t_spec>
void algorithm::Perturbation_MD<t_spec>
::G96Forcefield(io::InTopology &topo,
		io::InInput &input,
		io::Argument &args)
{
  DEBUG(7, "md: create perturbed forcefield");
  // check which interactions to add
  int do_bond, do_angle, do_dihedral, do_improper, do_nonbonded;
  input.read_FORCE(do_bond, do_angle, do_improper,
		   do_dihedral, do_nonbonded);

  int bond_term, angle_term;
  bool have_DIRK =
    input.read_FORCEFIELD(bond_term, angle_term);
  
  const std::vector<interaction::bond_type_struct> * bond_param = NULL;
  
  if (do_bond){

    if ((have_DIRK && (bond_term == 0))
	|| ((!have_DIRK) && (do_bond == 1))){
      
      if (do_bond == 2){
	io::messages.add("FORCEFIELD and FORCE block contradictory "
			 "for bond term", "md", io::message::error);
      }
      else
	io::messages.add("using Gromos96 quartic bond term", 
			 "md", io::message::notice);
      // bonds: quartic
      typename t_spec::quartic_bond_interaction_type * the_qbond_interaction =
	new typename t_spec::quartic_bond_interaction_type;
    
      topo >> *the_qbond_interaction;
      //bond_param = &m_qbond_interaction->parameter();

      m_forcefield.push_back(the_qbond_interaction); 
      
      typename t_spec::perturbed_quartic_bond_interaction_type *
	the_perturbed_qbond_interaction = 
	new typename t_spec::perturbed_quartic_bond_interaction_type
	(*the_qbond_interaction);

      m_forcefield.push_back(the_perturbed_qbond_interaction);
    
    }
    
    else if ((have_DIRK && (bond_term == 2))
	     || (do_bond == 2)){

      if (do_bond == 1){
	io::messages.add("FORCEFIELD and FORCE block contradictory "
			 "for bond term", "md", io::message::error);
      }
      else
	io::messages.add("using Gromos87 harmonic bond term", 
			 "md", io::message::notice);
      // bonds: harmonic
      typename t_spec::harmonic_bond_interaction_type *
	the_hbond_interaction =
	new typename t_spec::harmonic_bond_interaction_type;

      topo >> *the_hbond_interaction;
      // bond_param = &the_hbond_interaction->parameter();

      m_forcefield.push_back(the_hbond_interaction);
      typename t_spec::perturbed_harmonic_bond_interaction_type *
	the_perturbed_hbond_interaction = 
	new typename t_spec::perturbed_harmonic_bond_interaction_type
	(*the_hbond_interaction);
      
      m_forcefield.push_back(the_perturbed_hbond_interaction);
    }
    else{
      io::messages.add("FORCE or FORCEFIELD block wrong for bond term",
		       "md", io::message::error);
    }
  }
  
  if (do_angle){
    // angles

    if (have_DIRK && (angle_term != 0)){
      io::messages.add("FORCEFIELD harmonic angle term not supported",
		       "md", io::message::error);
    }

    typename t_spec::angle_interaction_type *the_angle_interaction = 
      new typename t_spec::angle_interaction_type;
        
    topo >> *the_angle_interaction;
 
    m_forcefield.push_back(the_angle_interaction);

    typename t_spec::perturbed_angle_interaction_type *
    the_perturbed_angle_interaction = 
    new typename t_spec::perturbed_angle_interaction_type
    (*the_angle_interaction);
  
    m_forcefield.push_back(the_perturbed_angle_interaction);
  }
  
  
  if (do_improper){
    // improper dihedrals
    DEBUG(7, "\timpropers");
    
    typename t_spec::improper_interaction_type *
      the_improper_interaction = 
      new typename t_spec::improper_interaction_type;

    topo >> *the_improper_interaction;
    
    m_forcefield.push_back(the_improper_interaction);
    DEBUG(7, "\tperturbed impropers");
    
    typename t_spec::perturbed_improper_interaction_type *
      the_perturbed_improper_interaction = 
      new typename t_spec::perturbed_improper_interaction_type
      (*the_improper_interaction);
    
    m_forcefield.push_back(the_perturbed_improper_interaction);
    DEBUG(7, "\tdone impropers");
    
  }
  
  if (do_dihedral){
    // dihedrals
    DEBUG(7, "\tdihedrals");
    
    typename t_spec::dihedral_interaction_type *
      the_dihedral_interaction =
      new typename t_spec::dihedral_interaction_type;
    
    topo >> *the_dihedral_interaction;
    
    m_forcefield.push_back(the_dihedral_interaction);
    DEBUG(7, "\tperturbed dihedrals");
    
    typename t_spec::perturbed_dihedral_interaction_type *
      the_perturbed_dihedral_interaction = 
      new typename t_spec::perturbed_dihedral_interaction_type
      (*the_dihedral_interaction);
    
    m_forcefield.push_back(the_perturbed_dihedral_interaction);
  }
  
  m_simulation.pressure_calculation(m_calculate_pressure);

  if (do_nonbonded){

    if (m_calculate_pressure){
      
      // nonbonded (with virial)
      DEBUG(8, "md (create_forcefield): nonbonded with pressure");

      typename t_spec::nonbonded_virial_interaction_type * 
	the_nonbonded_virial_interaction =
	new typename t_spec::nonbonded_virial_interaction_type(m_simulation);

      topo >> *the_nonbonded_virial_interaction;
      
      DEBUG(10, "md (create forcefield): nonbonded with pressure read in");

      m_forcefield.push_back(the_nonbonded_virial_interaction);
 
      typename t_spec::perturbed_nonbonded_virial_interaction_type *
	the_perturbed_nonbonded_interaction = 
	new typename t_spec::perturbed_nonbonded_virial_interaction_type 
	(m_simulation, *the_nonbonded_virial_interaction);
      
      m_forcefield.push_back(the_perturbed_nonbonded_interaction);
    }
    else{
      // nonbonded
      DEBUG(8, "md (create_forcefield): nonbonded without pressure");
      typename t_spec::nonbonded_interaction_type *
	the_nonbonded_interaction =
	new typename t_spec::nonbonded_interaction_type(m_simulation);
      
      topo >> *the_nonbonded_interaction;
      
      m_forcefield.push_back(the_nonbonded_interaction);

      typename t_spec::perturbed_nonbonded_interaction_type *
	the_perturbed_nonbonded_interaction = 
	new typename t_spec::perturbed_nonbonded_interaction_type
	(m_simulation, *the_nonbonded_interaction);
      
      m_forcefield.push_back(the_perturbed_nonbonded_interaction);

    }
    
  }

  // decide on SHAKE
  m_distance_constraint.init(m_simulation, args, topo, input);
  
  DEBUG(7, "forcefield created");
}

template<typename t_spec>
void algorithm::Perturbation_MD<t_spec>
::print_pairlists()
{
  
  typename std::vector<typename interaction::Interaction<
    typename t_spec::simulation_type> *>
    ::const_iterator it = m_forcefield.begin(),
    to = m_forcefield.end();
	
  for( ; it != to; ++it){
	  
    if ((*it)->name == "NonBonded"){
      
      if (m_calculate_pressure){
	std::cerr << "printing virial pairlist" << std::endl;

	(*m_print_file) << "shortrange\n" 
			<< dynamic_cast<
	  typename t_spec::nonbonded_virial_interaction_type *>
	  (*it)->pairlist()
			<< std::endl;
      }
      else {      
	std::cerr << "printing pairlist" << std::endl;
	
	(*m_print_file) << "shortrange\n" 
			<< dynamic_cast<typename t_spec::nonbonded_interaction_type *>
	  (*it)->pairlist()
			<< std::endl;
	
      }
      
    }
    
  }
  
}

template<typename t_spec>
void algorithm::Perturbation_MD<t_spec>
::do_energies()
{
  // do the normal energies
  MD<t_spec>::do_energies();

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
    io::print_ENERGY(std::cout, 
		     m_simulation.system().lambda_energies(),
		     m_simulation.topology().energy_groups(),
		     "dE/dLAMBDA");
  }
}
