/**
 * @file md.tcc
 * MD implementation
 */


#undef MODULE
#undef SUBMODULE
#define MODULE algorithm
#define SUBMODULE algorithm

#include "../../debug.h"

template<typename t_spec>
algorithm::MD<t_spec>
::MD()
  : title("\tgromosXX molecular dynamics program"),
    m_simulation(),
    m_forcefield(),
    m_temperature(),
    m_pressure(),
    m_distance_constraint(),
    m_trajectory(),
    m_print_file(&std::cout),
    m_dt(0),
    m_time(0),
    m_print_energy(1),
    m_print_pairlist(0),
    m_print_force(0),
    m_remove_com(0),
    m_print_com(0),
    m_calculate_pressure(0),
    m_do_perturbation(false)
{
}

template<typename t_spec>
algorithm::MD<t_spec>
::~MD()
{
  if (m_print_file != &std::cout){
    m_print_file->flush();
    delete m_print_file;
  }
  if (m_trajectory)
    delete m_trajectory;
}

template<typename t_spec>
int algorithm::MD<t_spec>::initialize(io::Argument &args)
{

  std::cout << "==============\n"
	    << "INITIALIZATION\n"
	    << "==============\n";

  //----------------------------------------------------------------------------

  DEBUG(7, "constructing topo, sys & input");
  io::InTopology topo;
  io::InTrajectory sys;
  io::InInput input;

  DEBUG(7, "open_files");
  open_files(args, topo, sys, input);

  // read input
  DEBUG(7, "read_input");
  read_input(args, topo, sys, input);
  
  DEBUG(7, "init_input");
  init_input(args, topo, sys, input);
  DEBUG(7, "init_output");
  init_output(args, input);
  
  //----------------------------------------------------------------------------
  // prepare for the run

  // and create the forcefield
  DEBUG(7, "md: create forcefield");
  G96Forcefield(topo, input, args);

  pre_md(input);
  
  //----------------------------------------------------------------------------
  
  // see whether everything is all right  
  m_simulation.check_state();

  // messages?
  std::cout << "MESSAGES (startup)\n";
  if (io::messages.display(std::cout) > io::message::warning)
    return 1;
  std::cout << "END\n";
  io::messages.clear();

  DEBUG(7, "pre md done");
  
  return 0;

}


template<typename t_spec>
inline int algorithm::MD<t_spec>
::pre_md(io::InInput &input)
{
  // last minute initializations...

  // prepare temperature calculation
  DEBUG(7, "md: degrees of freedom");
  m_simulation.calculate_degrees_of_freedom();

  // initialize the energy fluctuations
  m_simulation.system().energy_averages().
    resize(m_simulation.system().energies().bond_energy.size(),
	   m_simulation.system().energies().kinetic_energy.size());
  
  // initialize the lambda derivative averages
  // to avoid a segfault because kinetic vector is 0 size
  m_simulation.system().lambda_derivative_averages().
    resize(m_simulation.system().energies().bond_energy.size(),
	   m_simulation.system().energies().kinetic_energy.size());

  int ntx, init;
  unsigned int ig;
  double tempi;
  input.read_START(ntx, init, tempi, ig);  

  if (tempi != 0){
    DEBUG(7, "generating initial velocities with T=" << tempi);
    std::cout << "generating initial velocities with T=" << tempi << "\n";

    m_simulation.system().generate_velocities(tempi, 
					      m_simulation.topology().mass(),
					      ig);
  }
  
  else if (ntx == 1){
    DEBUG(7, "setting initial velocities to 0");
    m_simulation.system().vel() = 0.0;
    m_simulation.system().exchange_vel();
    m_simulation.system().vel() = 0.0;
  }

  // shake intial positions / velocities if required.
  init_pos_vel(init);
  
  // centre of mass removal
  int ndfmin, ntcm;
  input.read_CENTREOFMASS(ndfmin, ntcm, m_remove_com);

  double e_kin_trans, e_kin_rot;
  if(init<4 && ntcm) {
    
    std::cout << "stopping initial center of mass motion\n";
  
    m_simulation.remove_com_motion(m_dt, true, true, 
				   e_kin_trans, e_kin_rot);
    io::print_CENTREOFMASS(std::cout , e_kin_trans, e_kin_rot);
  }
  else{
    m_simulation.remove_com_motion(m_dt, false, false, 
				   e_kin_trans, e_kin_rot);
    io::print_CENTREOFMASS(std::cout , e_kin_trans, e_kin_rot);
    
  }

  temperature_algorithm().calculate_kinetic_energy(m_simulation);
  std::cout << "INITIAL TEMPERATURES AND TEMPERATURE COUPLING\n";
  io::print_DEGREESOFFREEDOM(std::cout, m_simulation.multibath());
  io::print_MULTIBATH_COUPLING(std::cout, m_simulation.multibath());
  io::print_MULTIBATH(std::cout, m_simulation.multibath(),
		      m_simulation.system().energies());


  m_trajectory->print_title(title);

  DEBUG(8, "md: put chargegroups into box");
  simulation().system().periodicity().
    put_chargegroups_into_box(simulation());
    
  return 0;
}


template<typename t_spec>
typename t_spec::simulation_type &  algorithm::MD<t_spec>
::simulation()
{
  return m_simulation;
}

template<typename t_spec>
interaction::Forcefield<typename t_spec::simulation_type> & 
algorithm::MD<t_spec>
::forcefield()
{
  return m_forcefield;
}

template<typename t_spec>
typename t_spec::temperature_type &  algorithm::MD<t_spec>
::temperature_algorithm()
{
  return m_temperature;
}

template<typename t_spec>
typename t_spec::pressure_type &  algorithm::MD<t_spec>
::pressure_algorithm()
{
  return m_pressure;
}

template<typename t_spec>
typename t_spec::distance_constraint_type &  algorithm::MD<t_spec>
::distance_constraint_algorithm()
{
  return m_distance_constraint;
}

template<typename t_spec>
typename t_spec::integration_type &  algorithm::MD<t_spec>
::integration_algorithm()
{
  return m_integration;
}

template<typename t_spec>
io::OutG96Trajectory<typename t_spec::simulation_type> &  
algorithm::MD<t_spec>
::trajectory()
{
  return *m_trajectory;
}

template<typename t_spec>
void algorithm::MD<t_spec>
::G96Forcefield(io::InTopology &topo,
		io::InInput &input, io::Argument &args)
{

  DEBUG(7, "md: create forcefield");
  // check which interactions to add
  int do_bond, do_angle, do_dihedral, do_improper, do_nonbonded;
  input.read_FORCE(do_bond, do_angle, do_improper,
		   do_dihedral, do_nonbonded);

  int bond_term, angle_term;
  bool have_DIRK =
    input.read_FORCEFIELD(bond_term, angle_term);

  // const std::vector<interaction::bond_type_struct> * bond_param = NULL;
  
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
      typename t_spec::quartic_bond_interaction_type *
	qbond_interaction =
	new typename t_spec::quartic_bond_interaction_type;      
    
      topo >> *qbond_interaction;
      // bond_param = &m_qbond_interaction->parameter();
      
      m_forcefield.push_back(qbond_interaction);
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
      typename t_spec::harmonic_bond_interaction_type
	*the_hbond_interaction =
	new typename t_spec::harmonic_bond_interaction_type;

      topo >> *the_hbond_interaction;
      m_forcefield.push_back(the_hbond_interaction);

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

    typename t_spec::angle_interaction_type
      *the_angle_interaction = 
      new     typename t_spec::angle_interaction_type;
    
    topo >> *the_angle_interaction;
 
    m_forcefield.push_back(the_angle_interaction);
  }
  
  if (do_improper){
    // improper dihedrals
    typename t_spec::improper_interaction_type
      *the_improper_interaction = 
      new typename t_spec::improper_interaction_type;

    topo >> *the_improper_interaction;
    
    m_forcefield.push_back(the_improper_interaction);
  }
  
  if (do_dihedral){
    // dihedrals
    typename t_spec::dihedral_interaction_type
      *the_dihedral_interaction =
      new typename t_spec::dihedral_interaction_type;
    
    topo >> *the_dihedral_interaction;
    
    m_forcefield.push_back(the_dihedral_interaction);
  }
  
  m_simulation.pressure_calculation(m_calculate_pressure);

  if (do_nonbonded){

    if (m_calculate_pressure){
      
      // nonbonded (with virial)
      DEBUG(8, "md (create_forcefield): nonbonded with pressure");

      typename t_spec::nonbonded_virial_interaction_type
	*the_nonbonded_virial_interaction =
	new typename t_spec::nonbonded_virial_interaction_type(m_simulation);

      topo >> *the_nonbonded_virial_interaction;
      
      DEBUG(10, "md (create forcefield): nonbonded with pressure read in");

      m_forcefield.push_back(the_nonbonded_virial_interaction);
    }
    else{
      // nonbonded
      DEBUG(8, "md (create_forcefield): nonbonded without pressure");
     typename t_spec::nonbonded_interaction_type
	* the_nonbonded_interaction =
	new typename t_spec::nonbonded_interaction_type(m_simulation);
      
      topo >> *the_nonbonded_interaction;
      
      m_forcefield.push_back(the_nonbonded_interaction);
    }
    
  }

  // decide on SHAKE
  m_distance_constraint.init(m_simulation, args, topo, input);

  DEBUG(7, "forcefield created");

}

/**
 * perform an MD simulation.
 */
template<typename t_spec>
int algorithm::MD<t_spec>::do_md(io::Argument &args)
{
  if(initialize(args)){
    return 1;
  }
  
  run();
  
  post_md();

  return 0;
}

/**
 * run an MD simulation.
 *
 * this is the main md loop.
 */
template<typename t_spec>
void algorithm::MD<t_spec>
::run(double time)
{

  std::cout << "========================\n";
  std::cout << "PERFORMING MD SIMULATION\n";
  std::cout << "========================\n";

  if (time == -1) time = m_time;
  double end_time = m_simulation.time() + time;
  
  while(m_simulation.time() < end_time){

    pre_step();
    
    do_step();

    post_step();
  
    DEBUG(8, "md: increase time");
    m_simulation.increase_time(m_dt);
     
  }    
}

template<typename t_spec>
void algorithm::MD<t_spec>
::pre_step()
{
  if (m_print_energy && m_simulation.steps() % m_print_energy == 0){
    io::print_TIMESTEP(std::cout, 
		       m_simulation.steps(), m_simulation.time());
  }
  
  DEBUG(8, "md: print trajectory");
  (*m_trajectory) << m_simulation;
}

template<typename t_spec>
void algorithm::MD<t_spec>
::do_step()
{
    // integrate
    DEBUG(8, "md: integrate");
    m_integration.step(m_simulation, m_forcefield, m_dt);
    
    if (m_print_pairlist && m_simulation.steps() % m_print_pairlist == 0){
      print_pairlists();
    }
      
    DEBUG(8, "md: shake");
    try{
      // std::cout << "shake solute:  " << 
      m_distance_constraint.solute(m_simulation.topology(), 
				   m_simulation.system(), m_dt);
      
            
      // std::cout << "shake solvent: " <<
      m_distance_constraint.solvent(m_simulation.topology(),
				    m_simulation.system(), m_dt);  
    }
    catch(std::runtime_error e){
      // go back to positions before SHAKE
      m_simulation.system().exchange_pos();
      (*m_trajectory) << m_simulation;
      throw;
    }
}

template<typename t_spec>
void algorithm::MD<t_spec>
::post_step()
{
  DEBUG(8, "md: calculate pressure");
  if (m_calculate_pressure){
    m_pressure.apply(m_simulation, m_dt);
  }
  if(m_remove_com && (m_simulation.steps()+1) % m_remove_com == 0){
    double ekin_trans, ekin_rot;
    
    DEBUG(8, "md: remove centre of mass");
    m_simulation.remove_com_motion(m_dt,true, true, ekin_trans, ekin_rot);
    if(m_print_com && (m_simulation.steps()+1 ) % m_print_com ==0){
      io::print_CENTREOFMASS(std::cout, ekin_trans, ekin_rot);
    }
  }
  else if(m_print_com &&( m_simulation.steps()+1 ) % m_print_com ==0){ 
    double ekin_trans, ekin_rot;
    DEBUG(8, "md: print centre of mass");
    
    m_simulation.remove_com_motion(m_dt,false, false, ekin_trans, ekin_rot);
    io::print_CENTREOFMASS(std::cout, ekin_trans, ekin_rot);
  }
  
  DEBUG(8, "md: calculate and print the energies");
  do_energies();
    
}

template<typename t_spec>
void algorithm::MD<t_spec>
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

template<typename t_spec>
void algorithm::MD<t_spec>
::open_files(io::Argument &args, io::InTopology &topo,
	     io::InTrajectory &sys, io::InInput &input)
{
  // read in the files - those are necessary
  DEBUG(7, "opening topology " << args["topo"]);
  std::ifstream *topo_file = new std::ifstream(args["topo"].c_str());
  if (!topo_file->good())
    io::messages.add("unable to open topology file: " + args["topo"], "md.tcc",
		     io::message::error);
  else 
    io::messages.add("parsing topology file: " + args["topo"], "md.tcc",
		     io::message::notice);

  topo.stream(*topo_file);
  DEBUG(7, "reading topology");
  topo.readStream();
  topo.auto_delete(true);
  
  DEBUG(7, "opening system");
  std::ifstream *sys_file = new std::ifstream(args["struct"].c_str());
  if (!sys_file->good())
    io::messages.add("unable to open initial structure file: " + args["struct"], 
                     "md.tcc",
		     io::message::error);
  else
    io::messages.add("parsing initial structure file: " + args["struct"], "md.tcc",
		     io::message::notice);
  sys.stream(*sys_file);
  sys.auto_delete(true);

  DEBUG(7, "opening input");
  std::ifstream *input_file = new std::ifstream(args["input"].c_str());
  if (!input_file->good())
    io::messages.add("unable to open input file: " + args["input"], 
                     "md.tcc",
		     io::message::error);
  else
    io::messages.add("parsing input file: " + args["input"], "md.tcc",
		     io::message::notice);
  input.stream(*input_file);
  DEBUG(7, "reading input");
  input.readStream();
  input.auto_delete(true);

  DEBUG(7, "read topology");
  topo.read_TOPOLOGY(m_simulation.topology());


  // decide whether we need velocities or not
  int ntx, init;
  unsigned int ig;
  double tempi;
  input.read_START(ntx, init, tempi, ig);  

  if (tempi != 0 || (tempi == 0 && ntx == 1)){
    sys.read_velocity = false;
    DEBUG(7, "not reading initial velocities from file");
  }
 
  int ntb, nrdbox;
  input.read_BOUNDARY(ntb, nrdbox);

  if (!nrdbox){
    switch(ntb){
      case 0:
	m_simulation.system().periodicity().boundary_condition(math::vacuum);
	io::messages.add("boundary conditions set to VACUUM", "InTrajectory", 
			 io::message::notice);
	break;
      case 1:
      case 2:
	m_simulation.system().periodicity().boundary_condition(math::triclinic);
	io::messages.add("boundary conditions set to TRICLINIC", "InTrajectory",
			 io::message::notice);
	break;
      default:
	throw std::runtime_error("bad boundary conditions in BOUNDARY block");
    }

    sys.read_box = false;
  }
  
  DEBUG(7, "read system");  
  sys >> m_simulation.system();

  // input, system and topology open and read.

}


template<typename t_spec>
void algorithm::MD<t_spec>
::init_input(io::Argument &args, io::InTopology &topo,
	     io::InTrajectory &sys, io::InInput &input)
{
  // add solvent
  DEBUG(7, "md: add solvent");
  int nsm;
  input.read_SYSTEM(nsm);
  if (nsm) m_simulation.solvate(0, nsm);

  // pressure coupling
  int ntp;
  double pres0, comp, tau;
  DEBUG(8, "md: read PCOUPLE");
  input.read_PCOUPLE(ntp, pres0, comp, tau);
  DEBUG(8, "md: PCOUPLE read");
  m_pressure.initialize(ntp, pres0, comp, tau);

}

template<typename t_spec>
void algorithm::MD<t_spec>
::read_input(io::Argument &args, io::InTopology &topo,
	     io::InTrajectory &sys, io::InInput &input)
{
  DEBUG(7, "parse print argument");
  parse_print_argument(args);


  DEBUG(7, "md: read input");
  input >> m_simulation;

  // pressure calculation
  DEBUG(7, "md: init pressure");
  int ntb, nrdbox;
  input.read_BOUNDARY(ntb, nrdbox);
  DEBUG(8, "md: boundary read");
  
  if (nrdbox != 1 && ntb != 0){
    io::messages.add("nrdbox!=1 only for vacuum runs supported","md.tcc",
		     io::message::error);
  }

  // do we need a virial calculation?
  if (abs(ntb) == 2)
    m_calculate_pressure = 1;
  else
    m_calculate_pressure = 0;

  // time to simulate
  int num_steps;
  double t0;
  input.read_STEP(num_steps, t0, m_dt);

  m_time = num_steps * m_dt;
  m_simulation.time(t0);
  
}

template<typename t_spec>
void algorithm::MD<t_spec>
::init_output(io::Argument &args, io::InInput &input)
{
  int print_trajectory, print_velocity_traj, print_energy_traj, 
    print_free_energy_traj;
  int conf_sel;
  int dihedral_monitoring;
  
  input.read_PRINT(m_print_energy, m_print_com, dihedral_monitoring);

  input.read_WRITE(print_trajectory, conf_sel, print_velocity_traj,
		   print_energy_traj, print_free_energy_traj);
  
  // some unhandled cases
  if (conf_sel != 0)
    io::messages.add("WRITE block: NTWSE != 0 not implemented",
		     "MD", io::message::error);

  std::ofstream *traj_file = new std::ofstream(args["trj"].c_str());
  std::ofstream *fin_file = new std::ofstream(args["fin"].c_str());

  // use a G96 trajectory
  m_trajectory = 
    new io::OutG96Trajectory<typename t_spec::simulation_type>
    (*traj_file, *fin_file, 
     print_trajectory, true);

  // optional files
  //================

  // velocity trajectory
  if (args.count("trv") == 1){
    std::ofstream *vel_file = new std::ofstream(args["trv"].c_str());
    m_trajectory->velocity_trajectory(*vel_file, print_velocity_traj);
  }

  // force trajectory
  if (args.count("trf") == 1){
    std::ofstream *force_file = new std::ofstream(args["trf"].c_str());
    m_trajectory->force_trajectory(*force_file, m_print_force);
  }

  // energy trajectory
  if (args.count("tre") == 1){
    std::ofstream *energy_file = new std::ofstream(args["tre"].c_str());
    m_trajectory->energy_trajectory(*energy_file, print_energy_traj);
  }

  // free energy trajectory
  if (args.count("trg") == 1){
    std::ofstream *free_energy_file = 
      new std::ofstream(args["trg"].c_str());
    m_trajectory->free_energy_trajectory(*free_energy_file, 
				    print_free_energy_traj);
  }
  
  // print file
  if (args.count("trp") == 1){
    m_print_file = new std::ofstream(args["trp"].c_str());
  }

}

template<typename t_spec>
void algorithm::MD<t_spec>
::print_pairlists()
{
  
  typename std::vector<typename interaction::Interaction<
    typename t_spec::simulation_type> *>
    ::const_iterator it = m_forcefield.begin(),
    to = m_forcefield.end();
	
  for( ; it != to; ++it){
	  
    if ((*it)->name == "NonBonded"){
      
      if (m_calculate_pressure){
	(*m_print_file) << "shortrange\n" 
			<< dynamic_cast<
	  typename t_spec::nonbonded_virial_interaction_type *>
	  (*it)->pairlist()
			<< std::endl;
      }
      else {      
	(*m_print_file) << "shortrange\n" 
			<< dynamic_cast<typename t_spec::nonbonded_interaction_type *>
	  (*it)->pairlist()
			<< std::endl;

      }
      
    }
    
  }
  
}

template<typename t_spec>
void algorithm::MD<t_spec>
::do_energies()
{
  // calculate the kinetic energy now (velocities adjusted for constraints)
  temperature_algorithm().calculate_kinetic_energy(m_simulation);
  // and sum up the energy arrays
  m_simulation.system().energies().calculate_totals();

  m_simulation.system().energy_averages().
    update(m_simulation.system().energies(), m_dt);
  
  if (m_print_energy && (m_simulation.steps()) % m_print_energy == 0){
    io::print_MULTIBATH(std::cout, m_simulation.multibath(),
			m_simulation.system().energies());
    io::print_ENERGY(std::cout, m_simulation.system().energies(),
		     m_simulation.topology().energy_groups());
    if (m_calculate_pressure)
      io::print_PRESSURE(std::cout, m_simulation.system());
  }
}

template<typename t_spec>
void algorithm::MD<t_spec>
::post_md()
{
  std::cout << "\nwriting final structure" << std::endl;
  trajectory() << io::final << simulation();
  
  simulation::Energy energy, fluctuation;

  simulation().system().energy_averages().
    average(energy, fluctuation);
  
  io::print_ENERGY(std::cout, energy,
		   simulation().topology().energy_groups(),
		   "AVERAGE ENERGIES");

  io::print_MULTIBATH(std::cout, simulation().multibath(),
		      energy);
  
  io::print_ENERGY(std::cout, fluctuation,
		   simulation().topology().energy_groups(),
		   "ENERGY FLUCTUATIONS");
  
  io::print_MULTIBATH(std::cout, simulation().multibath(),
		      fluctuation);

}

template<typename t_spec>
void algorithm::MD<t_spec>
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

