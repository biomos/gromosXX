/**
 * @file md_base.tcc
 * MD implementation
 */


#undef MODULE
#undef SUBMODULE
#define MODULE algorithm
#define SUBMODULE algorithm

#include "../../debug.h"

template<typename t_md_spec, typename t_interaction_spec>
algorithm::MD_Base<t_md_spec, t_interaction_spec>
::MD_Base()
  : title("\tgromosXX molecular dynamics program"),
    m_trajectory(),
    m_print_file(&std::cout),
    m_dt(0),
    m_time(0),
    m_print_energy(1),
    m_print_pairlist(0),
    m_print_force(0),
    m_remove_com(0),
    m_print_com(0),
    m_calculate_pressure(false),
    m_do_perturbation(false),
    m_simulation(),
    m_forcefield(),
    m_temperature(),
    m_pressure(),
    m_distance_constraint()
{
}

template<typename t_md_spec, typename t_interaction_spec>
algorithm::MD_Base<t_md_spec, t_interaction_spec>
::~MD_Base()
{
  m_print_file->flush();
  if (m_print_file != &std::cout){
    delete m_print_file;
  }
  if (m_trajectory)
    delete m_trajectory;
}

template<typename t_md_spec, typename t_interaction_spec>
typename t_md_spec::simulation_type &  algorithm::MD_Base<t_md_spec, t_interaction_spec>
::simulation()
{
  return m_simulation;
}

template<typename t_md_spec, typename t_interaction_spec>
interaction::Forcefield<typename t_md_spec::simulation_type, t_interaction_spec> & 
algorithm::MD_Base<t_md_spec, t_interaction_spec>
::forcefield()
{
  return m_forcefield;
}

template<typename t_md_spec, typename t_interaction_spec>
typename t_md_spec::temperature_type &  
algorithm::MD_Base<t_md_spec, t_interaction_spec>
::temperature_algorithm()
{
  return m_temperature;
}

template<typename t_md_spec, typename t_interaction_spec>
typename t_md_spec::pressure_type &  algorithm::MD_Base<t_md_spec, t_interaction_spec>
::pressure_algorithm()
{
  return m_pressure;
}

template<typename t_md_spec, typename t_interaction_spec>
typename t_md_spec::distance_constraint_type &  algorithm::MD_Base<t_md_spec, t_interaction_spec>
::distance_constraint_algorithm()
{
  return m_distance_constraint;
}

template<typename t_md_spec, typename t_interaction_spec>
typename t_md_spec::integration_type &  algorithm::MD_Base<t_md_spec, t_interaction_spec>
::integration_algorithm()
{
  return m_integration;
}

template<typename t_md_spec, typename t_interaction_spec>
void algorithm::MD_Base<t_md_spec, t_interaction_spec>
::parse_print_argument(io::Argument &args)
{

  m_print_pairlist = 0;
  m_print_force = 1;
  { // print
    io::Argument::const_iterator it = args.lower_bound("print"),
      to = args.upper_bound("print");
    if (it != to)
      std::cout << "printing\n";
    for( ; it != to; ++it){
      std::string s;
      int num;
      std::string::size_type sep = it->second.find(':');
      if (sep == std::string::npos){
	s = it->second;
	num = 1;
      }
      else{
	s = it->second.substr(0, sep);
	num = atoi(it->second.substr(sep+1, std::string::npos).c_str());
      }
      std::cout << "\t" << std::setw(15) << s << std::setw(6) << num << "\n";
      
      if (s == "pairlist") m_print_pairlist = num;
      else if (s == "force") m_print_force = num;
      else throw std::string("unknown @print argument");
    }
  }
}

template<typename t_md_spec, typename t_interaction_spec>
void algorithm::MD_Base<t_md_spec, t_interaction_spec>
::init_pos_vel(int init)
{

  if (init == 1){

    std::cout << "enforce position constraints for intial positions\n"
	      << "and remove initial velocity along the constraints\n";

    // shake positions and velocities!
    simulation().system().exchange_pos();
    simulation().system().exchange_vel();
    simulation().system().pos() = simulation().system().old_pos();
    simulation().system().vel() = simulation().system().old_vel();
    
    DEBUG(7, "shake initial coordinates -- solute");
    // shake the current coordinates with respect to themselves
    distance_constraint_algorithm().solute(simulation().topology(),
					   simulation().system(),
					   m_dt);
    DEBUG(7, "shake initial coordinates -- solvent");
    distance_constraint_algorithm().solvent(simulation().topology(),
					    simulation().system(),
					    m_dt);
  
    // restore the velocities
    simulation().system().vel() = simulation().system().old_vel();

    // take a step back (positions)
    DEBUG(10, "take a step back");
    simulation().system().exchange_pos();
    simulation().system().pos() = simulation().system().old_pos()
      - m_dt * simulation().system().vel();
    
    // shake again
    DEBUG(7, "shake initial velocities -- solute");
    distance_constraint_algorithm().solute(simulation().topology(),
					   simulation().system(),
					   m_dt);

    DEBUG(7, "shake initial velocities -- solvent");
    distance_constraint_algorithm().solvent(simulation().topology(),
					    simulation().system(),
					    m_dt);

    // restore the positions
    simulation().system().exchange_pos();

    // the velocities are negative (wrong time direction)
    simulation().system().vel() = -1.0 * simulation().system().vel();
    
    // and copy to old array
    simulation().system().exchange_vel();
    simulation().system().vel() = simulation().system().old_vel();
    
  }
  else if (init == 2){
    io::messages.add("init = 2 in START block no longer supported",
		     "md", io::message::error);
  }  
  
}

template<typename t_md_spec, typename t_interaction_spec>
io::OutG96Trajectory<typename t_md_spec::simulation_type> &  
algorithm::MD_Base<t_md_spec, t_interaction_spec>
::trajectory()
{
  return *m_trajectory;
}

