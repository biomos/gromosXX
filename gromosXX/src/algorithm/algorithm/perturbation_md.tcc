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
  
  int ntg, nlam;
  double rlam, dlamt, alphlj, alphc;
  input.read_PERTURB(ntg, rlam, dlamt, alphlj, alphc, nlam);

  if (ntg == 1){
    io::messages.add("Perturbation enabled", "Perturbation_MD",
		     io::message::notice);

    m_simulation.topology().lambda(rlam);
    m_simulation.topology().alpha_lj(alphlj);
    m_simulation.topology().alpha_crf(alphc);
    m_simulation.topology().nlam(nlam);
    
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

  // resize the energy derivative array
  m_simulation.system().lambda_energies().
    resize(m_simulation.system().energies().bond_energy.size());
  DEBUG(7, "lambda baths: " << m_simulation.system().energies().kinetic_energy.size());
  m_simulation.system().lambda_energies().
    kinetic_energy.resize(m_simulation.system().energies().kinetic_energy.size());
    
  // initialize topology for lambda = ??
  m_simulation.topology().update_for_lambda();
  
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
  DEBUG(7, "md: create forcefield");
  // check which interactions to add
  int do_bond, do_angle, do_dihedral, do_improper, do_nonbonded;
  input.read_FORCE(do_bond, do_angle, do_improper,
		   do_dihedral, do_nonbonded);

  const std::vector<interaction::bond_type_struct> * bond_param = NULL;
  
  if (do_bond == 1){
    io::messages.add("using Gromos96 quartic bond term", "vtest", io::message::notice);
    // bonds: quartic
    m_qbond_interaction =
      new interaction::Quartic_bond_interaction<t_simulation>;
    
    topo >> *m_qbond_interaction;
    bond_param = &m_qbond_interaction->parameter();

    m_forcefield.push_back(m_qbond_interaction); 

    interaction::Perturbed_Quartic_Bond_Interaction<t_simulation>
    *the_perturbed_qbond_interaction = 
    new interaction::Perturbed_Quartic_Bond_Interaction<t_simulation>
    (*m_qbond_interaction);

    m_forcefield.push_back(the_perturbed_qbond_interaction);
    
  }
  
  if (do_bond == 2){
    io::messages.add("using Gromos87 harmonic bond term", "vtest", io::message::notice);
    // bonds: harmonic
    interaction::harmonic_bond_interaction<t_simulation>
      *the_hbond_interaction =
      new interaction::harmonic_bond_interaction<t_simulation>;

    topo >> *the_hbond_interaction;
    bond_param = &the_hbond_interaction->parameter();

    m_forcefield.push_back(the_hbond_interaction);
  }

  if (do_angle){
    // angles
    interaction::angle_interaction<t_simulation> *the_angle_interaction = 
      new interaction::angle_interaction<t_simulation>;
    
    topo >> *the_angle_interaction;
 
    m_forcefield.push_back(the_angle_interaction);

    interaction::Perturbed_Angle_Interaction<t_simulation>
    *the_perturbed_angle_interaction = 
    new interaction::Perturbed_Angle_Interaction<t_simulation>
    (*the_angle_interaction);
  
    m_forcefield.push_back(the_perturbed_angle_interaction);
  }
  
  
  if (do_improper){
    // improper dihedrals
    interaction::Improper_dihedral_interaction<t_simulation>
      *the_improper_interaction = 
      new interaction::Improper_dihedral_interaction<t_simulation>;

    topo >> *the_improper_interaction;
    
    m_forcefield.push_back(the_improper_interaction);
  }
  
  if (do_dihedral){
    // dihedrals
    interaction::Dihedral_interaction<t_simulation>
      *the_dihedral_interaction =
      new interaction::Dihedral_interaction<t_simulation>;
    
    topo >> *the_dihedral_interaction;
    
    m_forcefield.push_back(the_dihedral_interaction);
  }
  
  m_simulation.pressure_calculation(m_calculate_pressure);

  if (do_nonbonded){

    if (m_calculate_pressure){
      
      // nonbonded (with virial)
      DEBUG(8, "md (create_forcefield): nonbonded with pressure");

      interaction::Nonbonded_Interaction<t_simulation,
	perturbed_pairlist_virial_type,
	perturbed_innerloop_virial_type> 
	*the_nonbonded_virial_interaction =
	new interaction::Nonbonded_Interaction<t_simulation,
	perturbed_pairlist_virial_type,
	perturbed_innerloop_virial_type>(m_simulation);

      
      topo >> *the_nonbonded_virial_interaction;
      
      DEBUG(10, "md (create forcefield): nonbonded with pressure read in");

      m_forcefield.push_back(the_nonbonded_virial_interaction);
      interaction::Perturbed_Nonbonded_Interaction<t_simulation, 
	perturbed_pairlist_virial_type, 
	perturbed_innerloop_virial_type> 
	* the_perturbed_nonbonded_interaction = 
	new interaction::Perturbed_Nonbonded_Interaction<t_simulation, 
	perturbed_pairlist_virial_type, 
	perturbed_innerloop_virial_type>
	(m_simulation, *the_nonbonded_virial_interaction);
      
      m_forcefield.push_back(the_perturbed_nonbonded_interaction);
    }
    else{
      // nonbonded
      DEBUG(8, "md (create_forcefield): nonbonded without pressure");
      interaction::Nonbonded_Interaction<t_simulation,
	perturbed_pairlist_type, 
	perturbed_innerloop_type> 
	*the_nonbonded_interaction =
	new interaction::Nonbonded_Interaction<t_simulation,
	perturbed_pairlist_type,
	perturbed_innerloop_type>(m_simulation);
      
      topo >> *the_nonbonded_interaction;
      
      m_forcefield.push_back(the_nonbonded_interaction);

      interaction::Perturbed_Nonbonded_Interaction<t_simulation, 
	perturbed_pairlist_type, 
	perturbed_innerloop_type> 
	* the_perturbed_nonbonded_interaction = 
	new interaction::Perturbed_Nonbonded_Interaction<t_simulation, 
	perturbed_pairlist_type, 
	perturbed_innerloop_type>
	(m_simulation, *the_nonbonded_interaction);
      
      m_forcefield.push_back(the_perturbed_nonbonded_interaction);

    }
    
  }

  DEBUG(7, "md (create forcefield): decide about SHAKE");

  // decide on SHAKE
  int ntc;
  double tolerance;
  input.read_SHAKE(ntc, tolerance);

  // bonds: harmonic
  interaction::harmonic_bond_interaction<t_simulation>
    *shake_param_interaction = NULL;

  if (ntc > 1 && bond_param == NULL){
    shake_param_interaction =
      new interaction::harmonic_bond_interaction<t_simulation>;

    topo >> *shake_param_interaction;    
    bond_param = &shake_param_interaction->parameter();
  }
  
  switch(ntc){
    case 1:
      break;
    case 2: 
      std::cout << "SHAKE bonds containing hydrogen atoms" << std::endl;
      m_simulation.topology().solute().
	add_bond_length_constraints(1.0,
				    m_simulation.topology().mass(),
				    *bond_param);
      break;
    case 3: 
      std::cout << "SHAKE all bonds" << std::endl;
      m_simulation.topology().solute().
	add_bond_length_constraints(*bond_param);
      break;
    default:
      std::cout << "wrong ntc" << std::endl;
  }

  if (shake_param_interaction) delete shake_param_interaction;

  DEBUG(7, "forcefield created");
}

template<typename t_simulation,
	 typename t_temperature,
	 typename t_pressure,
	 typename t_distance_constraint,
	 typename t_integration>
void algorithm::Perturbation_MD<t_simulation, t_temperature, t_pressure, 
				t_distance_constraint, t_integration>
::print_pairlists()
{
  
  typename std::vector<typename interaction::Interaction<
    typename parent_type::simulation_type> *>
    ::const_iterator it = m_forcefield.begin(),
    to = m_forcefield.end();
	
  for( ; it != to; ++it){
	  
    if ((*it)->name == "NonBonded"){
      
      if (m_calculate_pressure){
	std::cerr << "printing virial pairlist" << std::endl;

	(*m_print_file) << "shortrange\n" 
			<< dynamic_cast<interaction::Nonbonded_Interaction
	  <typename parent_type::simulation_type,
	  perturbed_pairlist_virial_type,
	  typename parent_type::innerloop_virial_type> *>
	  (*it)->pairlist()
			<< std::endl;
      }
      else {      
	std::cerr << "printing pairlist" << std::endl;
	
	(*m_print_file) << "shortrange\n" 
			<< dynamic_cast<interaction::Nonbonded_Interaction
	  <typename parent_type::simulation_type,
	  perturbed_pairlist_type,
	  typename parent_type::innerloop_type> *>
	  (*it)->pairlist()
			<< std::endl;
	
      }
      
    }
    
  }
  
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

  // and calculate the kinetic energy lamba derivative
  m_temperature.calculate_kinetic_energy_lambda_derivative(m_simulation);
  
  // and sum up the energy lambda derivative arrays
  m_simulation.system().lambda_energies().calculate_totals();
  
  if (m_print_energy && m_simulation.steps() % m_print_energy == 0){
    io::print_ENERGY(std::cout, m_simulation.system().lambda_energies(),
		     m_simulation.topology().energy_groups(), "dE/dLAMBDA");
  }
}
