/**
 * @file simulation.t.cc
 * check simple things in simulation.
 */

/**
 * check a pairlist.
 */
template<typename t_simulation>
int check_pairlist(std::ifstream &pl_file, 
		   interaction::simple_pairlist<t_simulation> &sr,
		   interaction::simple_pairlist<t_simulation> &lr)
{
  int i=-1, j=-1;
  string s;

  int errors = 0;

  typename interaction::simple_pairlist<t_simulation>::iterator *pl_it, *pl_to;
  typename interaction::simple_pairlist<t_simulation>::iterator 
    pl_sr_it = sr.begin(), pl_sr_to = sr.end();
  typename interaction::simple_pairlist<t_simulation>::iterator 
    pl_lr_it = lr.begin(), pl_lr_to = lr.end();

  while(pl_file >> s){

    // std::cout << s << std::endl;
    
    if (s == "") continue;
    
    if (s == "shortrange"){
      // std::cout << "shortrange" << std::endl;

      pl_it = &pl_sr_it;
      pl_to = &pl_sr_to;
      continue;
    }
    if (s == "longrange"){
      // std::cout << "longrange" << std::endl;

      pl_it = &pl_lr_it;
      pl_to = &pl_lr_to;
      continue;
    }
    
    if (s.find(':') == std::string::npos){
      j = atoi(s.c_str());
    }
    else{
      i = atoi(s.substr(0, s.find(':')).c_str());
      continue;
    }
    
    if (pl_it->i() != i || **pl_it != j){
      std::cout << "\nwrong: i=" << i << " j=" << j 
		<< " pl i=" << pl_it->i() << " pl j= " << **pl_it << std::endl;

      if (++errors == 3) return errors;
      
    }
    
    ++(*pl_it);
  }
  
  
  return errors;
}

/**
 * check forces
 */
int check_coordinates(std::ifstream &file,
		      std::string name,
		      math::VArray &coord,
		      const double epsilon)
{
  double f;
  int errors=0;
  std::string s = "";
  while (s != name)
    file >> s;

  std::cout.precision(18);

  for(int i=0; i<coord.size(); ++i){

    for(int j=0; j<3; ++j){
      
      file >> s;
      if (s[0] == '#'){
	getline(file,s);
	file >> s;
      }

      if (file.eof() || file.fail()){
	std::cout << "\n\treading of coord file failed." << std::endl;
	return ++errors;
      }
      
      f = atof(s.c_str());
      
      if (f - epsilon > coord(i)(j) || f + epsilon < coord(i)(j)){
	std::cout << "\n\tcoord " << i << "[" << j << "] = " << coord(i)(j)
		  << "\n\tshould be " << f << " +- " << epsilon << std::endl;

	if (++errors >= 1) return errors;
      }

    }

  }

  return errors;
  
}

int simulation_check()
{
  const double epsilon = 0.000000002;
  const int correct = 0;

  int result = 0;
  int last_result;
  int res;
  
  typedef simulation::System<math::any> system_type;
  typedef simulation::Topology topology_type;
  typedef simulation::Simulation<topology_type, system_type> simulation_type;
  typedef interaction::Forcefield<simulation_type> forcefield_type;
  typedef interaction::twin_range_pairlist_cg<simulation_type> pairlist_type;

  try{
    
    system_type the_system;
    topology_type the_topology;
    simulation_type the_simulation(the_topology, the_system);
    
    // read in a topology and a system and the input
    std::ifstream topo_file;
    GETFILE(topo_file, "rasn.topo");
    std::ifstream coord_file;
    GETFILE(coord_file, "rasn.coord");
    std::ifstream input_file;
    GETFILE(input_file, "rasn.in");
    
    io::InTopology topo(topo_file);
    io::InTrajectory coord(coord_file);
    io::InInput input(input_file);

    topo >> the_topology;
    coord >> the_system;

    CHECKING("InTopology", last_result);
    
    CHECK_EQUAL(the_topology.num_solute_atoms() + 1244 * 3,
		unsigned(the_system.pos().size()), last_result);
    
    CHECK_EQUAL(the_topology.num_solvents(), 0, last_result);
    
    CHECK_EQUAL(the_topology.num_chargegroups(), 3, last_result);
    
    CHECK_EQUAL(the_topology.num_solute_chargegroups(), 3, last_result);
  
    RESULT(last_result, result);

    CHECKING("InTrajectory", last_result);
    
    CHECK_EQUAL(the_system.periodicity().boundary_condition(), math::triclinic,
		last_result);
    
    RESULT(last_result, result);
    
    forcefield_type the_forcefield;
    
    // bonds: quartic
    interaction::Quartic_bond_interaction<simulation_type> 
      *the_qbond_interaction =
      new interaction::Quartic_bond_interaction<simulation_type>;

    // bonds: harmonic
    interaction::harmonic_bond_interaction<simulation_type>
      *the_hbond_interaction =
      new interaction::harmonic_bond_interaction<simulation_type>;

    // angles
    interaction::angle_interaction<simulation_type>
      *the_angle_interaction = 
      new interaction::angle_interaction<simulation_type>;
  
    // improper dihedrals
    interaction::Improper_dihedral_interaction<simulation_type>
      *the_improper_interaction = 
      new interaction::Improper_dihedral_interaction<simulation_type>;
    
    // dihedrals
    interaction::Dihedral_interaction<simulation_type>
      *the_dihedral_interaction =
      new interaction::Dihedral_interaction<simulation_type>;
    
    // nonbonded
    interaction::Nonbonded_Interaction<simulation_type, pairlist_type>
      *the_nonbonded_interaction =
      new interaction::Nonbonded_Interaction<simulation_type, pairlist_type>;
    
    interaction::Nonbonded_Virial_Interaction<simulation_type, pairlist_type>
      *the_nonbonded_virial_interaction =
      new interaction::Nonbonded_Virial_Interaction<simulation_type, pairlist_type>;
    
    // read parameter
    topo >> *the_qbond_interaction;
    topo >> *the_hbond_interaction;
    topo >> *the_angle_interaction;
    topo >> *the_improper_interaction;
    topo >> *the_dihedral_interaction;
    topo >> *the_nonbonded_interaction;
    topo >> *the_nonbonded_virial_interaction;

    input >> the_simulation;
    // solvate has to be called after reading the input (at least for now...)
    // because the SUBMOLECULES block has to be ordered...
    the_simulation.solvate(0, 1244);

    the_simulation.check_state();

    // messages?
    std::cout << "Messages (startup)\n";
    if (io::messages.display(std::cout) > io::message::warning)
      return 1;
    std::cout << "\n";
    io::messages.clear();

    // create a chargegroup based pairlist and compare
    pairlist_type the_pairlist;
    the_pairlist.update(the_simulation);
    
    std::ifstream corr_pl;
    GETFILE(corr_pl, "rasn.pl");
    
    CHECKING("pairlist construction", last_result);

    res = check_pairlist(corr_pl, 
			 the_pairlist.short_range(), 
			 the_pairlist.long_range());

    CHECK_EQUAL(res, correct, last_result);

    RESULT(last_result, result);
    
    std::ifstream check_file;
    GETFILE(check_file, "rasn.qbond");
    
    CHECKING("quartic bond forces", last_result);

    the_forcefield.push_back(the_qbond_interaction);
    the_forcefield.calculate_interactions(the_simulation);
    
    res = check_coordinates(check_file, "FORCERED", the_system.force(), epsilon);
    
    CHECK_EQUAL(res, correct, last_result);
    RESULT(last_result, result);

    the_forcefield.clear();
    check_file.close();

    CHECKING("harmonic bond forces", last_result);

    the_forcefield.push_back(the_hbond_interaction);
    the_forcefield.calculate_interactions(the_simulation);

    GETFILE(check_file, "rasn.hbond");
    
    res = check_coordinates(check_file, "FORCERED", the_system.force(), epsilon);
    
    CHECK_EQUAL(res, correct, last_result);
    RESULT(last_result, result);
    
    the_forcefield.clear();
    check_file.close();


    CHECKING("bond angle forces", last_result);

    the_forcefield.push_back(the_angle_interaction);
    the_forcefield.calculate_interactions(the_simulation);
    
    GETFILE(check_file, "rasn.angle");
    
    res = check_coordinates(check_file, "FORCERED", the_system.force(), epsilon);
    
    CHECK_EQUAL(res, correct, last_result);
    RESULT(last_result, result);
    
    the_forcefield.clear();
    check_file.close();


    CHECKING("improper forces", last_result);

    the_forcefield.push_back(the_improper_interaction);
    the_forcefield.calculate_interactions(the_simulation);
    
    GETFILE(check_file, "rasn.improper");
    
    res = check_coordinates(check_file, "FORCERED", the_system.force(), epsilon);
    
    CHECK_EQUAL(res, correct, last_result);
    RESULT(last_result, result);
    
    the_forcefield.clear();
    check_file.close();


    CHECKING("dihedral forces", last_result);

    the_forcefield.push_back(the_dihedral_interaction);
    the_forcefield.calculate_interactions(the_simulation);

    GETFILE(check_file, "rasn.dihedral");
    
    res = check_coordinates(check_file, "FORCERED", the_system.force(), epsilon);
    
    CHECK_EQUAL(res, correct, last_result);
    RESULT(last_result, result);
    
    the_forcefield.clear();
    check_file.close();


    CHECKING("nonbonded forces", last_result);

    the_forcefield.push_back(the_nonbonded_interaction);
    the_forcefield.calculate_interactions(the_simulation);
    
    GETFILE(check_file, "rasn.nonbonded");
    
    res = check_coordinates(check_file, "FORCERED", the_system.force(), epsilon);
    
    CHECK_EQUAL(res, correct, last_result);
    RESULT(last_result, result);
    
    the_forcefield.clear();
    check_file.close();


    CHECKING("molecular translational ekin", last_result);
    
    the_forcefield.push_back(the_nonbonded_virial_interaction);
    the_forcefield.calculate_interactions(the_simulation);
    
    // check whether we get the same forces
    GETFILE(check_file, "rasn.nonbonded");
    res = check_coordinates(check_file, "FORCERED", the_system.force(), epsilon);
    CHECK_EQUAL(res, correct, last_result);

    // check the virial
    /*
    std::cout << "pressure: " << the_simulation.system().pressure() << std::endl
	      << "virial: " << the_simulation.system().virial() << std::endl
	      << "molecular ekin: " << the_simulation.system().molecular_kinetic_energy() << std::endl;
    */

    math::Matrix &mol_ekin = the_simulation.system().molecular_kinetic_energy();

    CHECK_APPROX_EQUAL(mol_ekin(0,0), 1606.60325, 0.00001,  last_result);
    CHECK_APPROX_EQUAL(mol_ekin(0,1),  -64.84878, 0.00001,  last_result);
    CHECK_APPROX_EQUAL(mol_ekin(0,2),  -16.06713, 0.00001,  last_result);
    CHECK_APPROX_EQUAL(mol_ekin(1,0),  -64.84878, 0.00001,  last_result);
    CHECK_APPROX_EQUAL(mol_ekin(1,1), 1552.48036, 0.00001,  last_result);
    CHECK_APPROX_EQUAL(mol_ekin(1,2),  -22.77930, 0.00001,  last_result);
    CHECK_APPROX_EQUAL(mol_ekin(2,0),  -16.06713, 0.00001,  last_result);
    CHECK_APPROX_EQUAL(mol_ekin(2,1),  -22.77930, 0.00001,  last_result);
    CHECK_APPROX_EQUAL(mol_ekin(2,2), 1589.80563, 0.00001,  last_result);

    RESULT(last_result, result);

    CHECKING("virial calculation", last_result);

    math::Matrix &virial = the_simulation.system().virial();
    
    CHECK_APPROX_EQUAL(virial(0,0), -4416.96163, 0.00001,  last_result);
    CHECK_APPROX_EQUAL(virial(0,1),  -612.72666, 0.00001,  last_result);
    CHECK_APPROX_EQUAL(virial(0,2),  -461.87316, 0.00001,  last_result);
    CHECK_APPROX_EQUAL(virial(1,0),   -84.73985, 0.00001,  last_result);
    CHECK_APPROX_EQUAL(virial(1,1), -5210.89030, 0.00001,  last_result);
    CHECK_APPROX_EQUAL(virial(1,2), -1006.66333, 0.00001,  last_result);
    CHECK_APPROX_EQUAL(virial(2,0),  -462.23679, 0.00001,  last_result);
    CHECK_APPROX_EQUAL(virial(2,1),  -369.95442, 0.00001,  last_result);
    CHECK_APPROX_EQUAL(virial(2,2), -8084.57680, 0.00001,  last_result);

    // and the total pressure... -> after Berendsen_Barostat is called...
    // math::Matrix &pressure = the_simulation.system().pressure();
    // double press = pressure(0,0) + pressure(1,1) + pressure(2,2);
    // CHECK_APPROX_EQUAL(press, 402.73872, 0.00001, last_result);

    RESULT(last_result, result);
    
    the_forcefield.clear();

  }
  catch(std::runtime_error e){
    std::cout << "runtime error:\n";
    std::cout << e.what() << std::endl;
    throw;
  }
  catch(std::string s){
    std::cout << "argument error:\n";
    std::cout << s << std::endl;
    throw;
  }
  
  return result;
}
