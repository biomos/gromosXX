/**
 * @file md.tcc
 * MD implementation
 */

template<typename t_simulation,
	 typename t_temperature,
	 typename t_pressure,
	 typename t_distance_constraint,
	 typename t_integration>
algorithm::MD<t_simulation, t_temperature, t_pressure,
	      t_distance_constraint, t_integration>
::MD(t_simulation &sim)
  : m_simulation(sim),
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
    m_calculate_pressure(0),
    m_qbond_interaction(NULL)
{
}

template<typename t_simulation,
	 typename t_temperature,
	 typename t_pressure,
	 typename t_distance_constraint,
	 typename t_integration>
algorithm::MD<t_simulation, t_temperature, t_pressure, t_distance_constraint, 
	      t_integration>
::~MD()
{
  if (m_print_file != &std::cout){
    m_print_file->flush();
    delete m_print_file;
  }
  if (m_trajectory)
    delete m_trajectory;
}

template<typename t_simulation,
	 typename t_temperature,
	 typename t_pressure,
	 typename t_distance_constraint,
	 typename t_integration>
int algorithm::MD<t_simulation, t_temperature, t_pressure, t_distance_constraint, t_integration>::initialize(io::Argument &args)
{

  //----------------------------------------------------------------------------
  // read input
  io::InTopology topo;
  io::InTrajectory sys;
  io::InInput input;
  
  init_input(args, topo, sys, input);
  
  //----------------------------------------------------------------------------
  // prepare for the output
  init_output(args, input);

  //----------------------------------------------------------------------------
  // prepare for the run

  // read in the input
  read_input(args, topo, sys, input);

  // and create the forcefield
  DEBUG(7, "md: create forcefield");
  G96Forcefield(topo, input);
  
  //----------------------------------------------------------------------------
  // see whether everything is all right  
  m_simulation.check_state();

  // messages?
  std::cout << "Messages (startup)\n";
  if (io::messages.display(std::cout) > io::message::warning)
    return 1;
  std::cout << "\n";
  io::messages.clear();

  DEBUG(7, "md initialized");
  return 0;

}


template<typename t_simulation,
	 typename t_temperature,
	 typename t_pressure,
	 typename t_distance_constraint,
	 typename t_integration>
t_simulation &  algorithm::MD<t_simulation, t_temperature, t_pressure,
			      t_distance_constraint, t_integration>
::simulation()
{
  return m_simulation;
}

template<typename t_simulation,
	 typename t_temperature,
	 typename t_pressure,
	 typename t_distance_constraint,
	 typename t_integration>
interaction::Forcefield<t_simulation> & 
algorithm::MD<t_simulation, t_temperature, t_pressure, 
	      t_distance_constraint, t_integration>
::forcefield()
{
  return m_forcefield;
}

template<typename t_simulation,
	 typename t_temperature,
	 typename t_pressure,
	 typename t_distance_constraint,
	 typename t_integration>
t_temperature &  algorithm::MD<t_simulation, t_temperature, t_pressure,
			       t_distance_constraint, t_integration>
::temperature_algorithm()
{
  return m_temperature;
}

template<typename t_simulation,
	 typename t_temperature,
	 typename t_pressure,
	 typename t_distance_constraint,
	 typename t_integration>
t_pressure &  algorithm::MD<t_simulation, t_temperature, t_pressure,
			    t_distance_constraint, t_integration>
::pressure_algorithm()
{
  return m_pressure;
}

template<typename t_simulation,
	 typename t_temperature,
	 typename t_pressure,
	 typename t_distance_constraint,
	 typename t_integration>
t_distance_constraint &  algorithm::MD<t_simulation, t_temperature, t_pressure,
				       t_distance_constraint, t_integration>
::distance_constraint_algorithm()
{
  return m_distance_constraint;
}

template<typename t_simulation,
	 typename t_temperature,
	 typename t_pressure,
	 typename t_distance_constraint,
	 typename t_integration>
t_integration &  algorithm::MD<t_simulation, t_temperature, t_pressure,
			       t_distance_constraint, t_integration>
::integration_algorithm()
{
  return m_integration;
}

template<typename t_simulation,
	 typename t_temperature,
	 typename t_pressure,
	 typename t_distance_constraint,
	 typename t_integration>
io::OutG96Trajectory<t_simulation> &  
algorithm::MD<t_simulation, t_temperature, t_pressure,
	      t_distance_constraint, t_integration>
::trajectory()
{
  return *m_trajectory;
}

template<typename t_simulation,
	 typename t_temperature,
	 typename t_pressure,
	 typename t_distance_constraint,
	 typename t_integration>
void algorithm::MD<t_simulation, t_temperature, t_pressure, 
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
    interaction::angle_interaction<t_simulation>
      *the_angle_interaction = 
      new interaction::angle_interaction<t_simulation>;
    
    topo >> *the_angle_interaction;
  
    m_forcefield.push_back(the_angle_interaction);
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
  
  if (do_nonbonded){

    if (m_calculate_pressure){
      // nonbonded (with virial)
      DEBUG(8, "md (create_forcefield): nonbonded with pressure");
      interaction::Nonbonded_Virial_Interaction<t_simulation, 
	interaction::twin_range_pairlist_cg<t_simulation> >
	*the_nonbonded_interaction =
	new interaction::Nonbonded_Virial_Interaction<t_simulation,
	interaction::twin_range_pairlist_cg<t_simulation> >;
    
      topo >> *the_nonbonded_interaction;
      
      DEBUG(10, "md (create forcefield): nonbonded with pressure read in");

      m_forcefield.push_back(the_nonbonded_interaction);
    }
    else{
      // nonbonded
      DEBUG(8, "md (create_forcefield): nonbonded without pressure");
      interaction::Nonbonded_Interaction<t_simulation, 
	interaction::twin_range_pairlist_cg<t_simulation> >
	*the_nonbonded_interaction =
	new interaction::Nonbonded_Interaction<t_simulation,
	interaction::twin_range_pairlist_cg<t_simulation> >;
      
      topo >> *the_nonbonded_interaction;
      
      m_forcefield.push_back(the_nonbonded_interaction);
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


/**
 * run an MD simulation.
 *
 * this is the main md loop.
 */
template<typename t_simulation,
	 typename t_temperature,
	 typename t_pressure,
	 typename t_distance_constraint,
	 typename t_integration>
void algorithm::MD<t_simulation, t_temperature, t_pressure, 
		   t_distance_constraint, t_integration>
::run(double time)
{
  if (time == -1) time = m_time;
  
  double end_time = m_simulation.time() + time;
  
  while(m_simulation.time() < end_time){

    io::print_TIMESTEP(std::cout, m_simulation.steps(), m_simulation.time());
    
    DEBUG(8, "md: put chargegroups into box");
    simulation().system().periodicity().
      put_chargegroups_into_box(simulation());
    

    DEBUG(8, "md: print trajectory");
    (*m_trajectory) << m_simulation;

    // integrate
    DEBUG(8, "md: integrate");
    m_integration.step(m_simulation, m_forcefield, m_dt);
    
    if (m_print_pairlist && m_simulation.steps() % m_print_pairlist == 0){
      print_pairlists();
    }
      
    DEBUG(8, "md: shake");
    try{
      std::cout << "shake solute:  " 
		<< m_distance_constraint.solute(m_simulation.topology(), 
						m_simulation.system(), m_dt)
		<< "\n";
      
      std::cout << "shake solvent: " 
		<< m_distance_constraint.solvent(m_simulation.topology(),
						 m_simulation.system(), m_dt)
		<< "\n";
    }
    catch(std::runtime_error e){
      // go back to positions before SHAKE
      m_simulation.system().exchange_pos();
      (*m_trajectory) << m_simulation;
      throw;
    }

    DEBUG(8, "md: calculate pressure");
    if (m_calculate_pressure){
      m_pressure.apply(m_simulation, m_dt);
    }
    
    DEBUG(8, "md: calculate and print the energies");
    do_energies();
    
    DEBUG(8, "md: increase time");
    m_simulation.increase_time(m_dt);
   
  }    
}

template<typename t_simulation,
	 typename t_temperature,
	 typename t_pressure,
	 typename t_distance_constraint,
	 typename t_integration>
void algorithm::MD<t_simulation, t_temperature, t_pressure, 
		   t_distance_constraint, t_integration>
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

template<typename t_simulation,
	 typename t_temperature,
	 typename t_pressure,
	 typename t_distance_constraint,
	 typename t_integration>
void algorithm::MD<t_simulation, t_temperature, t_pressure, 
		   t_distance_constraint, t_integration>
::open_files(io::Argument &args, io::InTopology &topo,
	     io::InTrajectory &sys, io::InInput &input)
{
  // std::cerr << "read in files" << std::endl;

  // read in the files - those are necessary
  std::ifstream *topo_file = new std::ifstream(args["topo"].c_str());
  if (!topo_file->good())
    io::messages.add("unable to open topology file: " + args["topo"], "md.tcc",
		     io::message::error);
  else 
    io::messages.add("parsing topology file: " + args["topo"], "md.tcc",
		     io::message::notice);
  topo.stream(*topo_file);
  topo.readStream();
  topo.auto_delete(true);
  
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

  std::ifstream *input_file = new std::ifstream(args["input"].c_str());
  if (!input_file->good())
    io::messages.add("unable to open input file: " + args["input"], 
                     "md.tcc",
		     io::message::error);
  else
    io::messages.add("parsing input file: " + args["input"], "md.tcc",
		     io::message::notice);
  input.stream(*input_file);
  input.readStream();
  input.auto_delete(true);
  
}

  
template<typename t_simulation,
	 typename t_temperature,
	 typename t_pressure,
	 typename t_distance_constraint,
	 typename t_integration>
void algorithm::MD<t_simulation, t_temperature, t_pressure, 
		   t_distance_constraint, t_integration>
::init_input(io::Argument &args, io::InTopology &topo,
	     io::InTrajectory &sys, io::InInput &input)
{
  // std::cerr << "init input" << std::endl;
  
  parse_print_argument(args);

  open_files(args, topo, sys, input);

  topo >> m_simulation.topology();
  sys >> m_simulation.system();
  
}

template<typename t_simulation,
	 typename t_temperature,
	 typename t_pressure,
	 typename t_distance_constraint,
	 typename t_integration>
void algorithm::MD<t_simulation, t_temperature, t_pressure, 
		   t_distance_constraint, t_integration>
::read_input(io::Argument &args, io::InTopology &topo,
	     io::InTrajectory &sys, io::InInput &input)
{
  DEBUG(7, "md: read input");
  input >> m_simulation;

  // add solvent
  DEBUG(7, "md: add solvent");
  int nsm;
  input.read_SYSTEM(nsm);
  if (nsm) m_simulation.solvate(0, nsm);

  // initialize SHAKE / ??
  DEBUG(7, "md init shake");
  int ntc;
  double tolerance;
  input.read_SHAKE(ntc, tolerance);
  m_distance_constraint.tolerance(tolerance);

  // pressure calculation
  DEBUG(7, "md: init pressure");
  int ntb, nrdbox;
  input.read_BOUNDARY(ntb, nrdbox);
  DEBUG(8, "md: boundary read");
  
  if (nrdbox != 1 && ntb != 0){
    io::messages.add("nrdbox!=1 only for vacuum runs supported","md.tcc",
		     io::message::error);
  }
  if (abs(ntb) == 2)
    m_calculate_pressure = 1;
  else
    m_calculate_pressure = 0;

  // prepare temperature calculation
  DEBUG(7, "md: degrees of freedom");
  m_simulation.calculate_degrees_of_freedom();
  temperature_algorithm().calculate_kinetic_energy(m_simulation);
  std::cout << "initial temperature:\n";
  io::print_MULTIBATH(std::cout, m_simulation.multibath());

  // pressure coupling (has to be done before 
  // constructing the forcefield!)
  int ntp;
  double pres0, comp, tau;
  DEBUG(8, "md: read PCOUPLE");
  input.read_PCOUPLE(ntp, pres0, comp, tau);
  DEBUG(8, "md: PCOUPLE read");
  
  m_pressure.initialize(ntp, pres0, comp, tau);

  // initialize the energy fluctuations
  m_simulation.system().energy_averages().
    resize(m_simulation.system().energies().bond_energy.size(),
	   m_simulation.system().energies().kinetic_energy.size());

  // time to simulate
  int num_steps;
  double t0;
  input.read_STEP(num_steps, t0, m_dt);

  m_time = num_steps * m_dt;
  m_simulation.time(t0);

}

template<typename t_simulation,
	 typename t_temperature,
	 typename t_pressure,
	 typename t_distance_constraint,
	 typename t_integration>
void algorithm::MD<t_simulation, t_temperature, t_pressure, 
		   t_distance_constraint, t_integration>
::init_output(io::Argument &args, io::InInput &input)
{
  // std::cerr << "init_output" << std::endl;
  
  int print_trajectory, print_velocity, print_energy_traj;
  input.read_PRINT(print_trajectory, print_velocity, m_print_energy);

  // std::cerr << "PRINT read" << std::endl;
  
  std::ofstream *traj_file = new std::ofstream(args["trj"].c_str());
  std::ofstream *fin_file = new std::ofstream(args["fin"].c_str());

  // std::cerr << "traj + fin file open" << std::endl;

  // use a G96 trajectory
  m_trajectory = 
    new io::OutG96Trajectory<simulation_type>(*traj_file, *fin_file, 
					      print_trajectory, true);

  // std::cerr << "OutG96Trajectory generated" << std::endl;  

  // optional files
  // velocity trajectory
  if (args.count("trv") == 1){
    // m_velocity_file.open(args["trv"].c_str());
    std::ofstream *vel_file = new std::ofstream(args["trv"].c_str());
    m_trajectory->velocity_trajectory(*vel_file, print_velocity);
    // std::cerr << "trv added" << std::endl;
  }

  // force trajectory
  if (args.count("trf") == 1){
    // m_force_file.open(args["trf"].c_str());
    std::ofstream *force_file = new std::ofstream(args["trf"].c_str());
    m_trajectory->force_trajectory(*force_file, m_print_force);
    // std::cerr << "trf added" << std::endl;
  }

  if (args.count("tre") == 1){
    // m_energy_file.open(args["tre"].c_str());
    std::ofstream *energy_file = new std::ofstream(args["tre"].c_str());
    m_trajectory->energy_trajectory(*energy_file, m_print_energy);
    // std::cerr << "tre added" << std::endl;
  }
  
  if (args.count("trp") == 1){
    m_print_file = new std::ofstream(args["trp"].c_str());
    // std::cerr << "trp added" << std::endl;
  }
  // std::cerr << "init output done" << std::endl;
  
}

template<typename t_simulation,
	 typename t_temperature,
	 typename t_pressure,
	 typename t_distance_constraint,
	 typename t_integration>
void algorithm::MD<t_simulation, t_temperature, t_pressure, 
		   t_distance_constraint, t_integration>
::print_pairlists()
{
  
  typename std::vector<typename interaction::Interaction<simulation_type> *>
    ::const_iterator it = m_forcefield.begin(),
    to = m_forcefield.end();
	
  for( ; it != to; ++it){
	  
    if ((*it)->name == "NonBonded"){

      (*m_print_file) << "shortrange\n" 
		      << dynamic_cast<interaction::Nonbonded_Interaction
	<simulation_type, pairlist_type> *>
	(*it)->pairlist().short_range()
		      << std::endl;

      (*m_print_file) << "longrange\n" 
		      << dynamic_cast<interaction::Nonbonded_Interaction
	<simulation_type, pairlist_type> *>
	(*it)->pairlist().long_range()
		      << std::endl;
	
    }
	  
  }
  
}

template<typename t_simulation,
	 typename t_temperature,
	 typename t_pressure,
	 typename t_distance_constraint,
	 typename t_integration>
void algorithm::MD<t_simulation, t_temperature, t_pressure, 
		   t_distance_constraint, t_integration>
::do_energies()
{
  // calculate the kinetic energy now (velocities adjusted for constraints)
  temperature_algorithm().calculate_kinetic_energy(m_simulation);
  // and sum up the energy arrays
  m_simulation.system().energies().calculate_totals();
  m_simulation.system().energy_averages().update(m_simulation.system().energies(), m_dt);
  
  if (m_print_energy && m_simulation.steps() % m_print_energy == 0){
    io::print_MULTIBATH(std::cout, m_simulation.multibath());
    io::print_ENERGY(std::cout, m_simulation.system().energies(),
		     m_simulation.topology().energy_groups());
    if (m_calculate_pressure)
      io::print_PRESSURE(std::cout, m_simulation.system());
  }
}
