/**
 * @file create_forcefield.h
 * creates a (G96) forcefield according to
 * input parameters and topology.
 */

template<typename t_simulation, typename t_pairlist >
interaction::Forcefield<t_simulation> & 
create_forcefield(io::InTopology &topo,
		  io::InInput &input,
		  typename t_simulation::topology_type &the_topology)
{

  // get a forcefield
  interaction::Forcefield<t_simulation> *the_forcefield
    = new interaction::Forcefield<t_simulation>;
  
  // check which interactions to add
  int do_bond, do_angle, do_dihedral, do_improper, do_nonbonded;
  input.read_FORCE(do_bond, do_angle, do_improper,
		   do_dihedral, do_nonbonded);

  const std::vector<interaction::bond_type_struct> * bond_param = NULL;
  
  if (do_bond == 1){
    io::messages.add("using Gromos96 quartic bond term", "vtest", io::message::notice);
    // bonds: quartic
    interaction::Quartic_bond_interaction<t_simulation> 
      *the_qbond_interaction =
      new interaction::Quartic_bond_interaction<t_simulation>;
    
    topo >> *the_qbond_interaction;
    bond_param = &the_qbond_interaction->parameter();

    the_forcefield->push_back(the_qbond_interaction);
  }

  if (do_bond == 2){
    io::messages.add("using Gromos87 harmonic bond term", "vtest", io::message::notice);
    // bonds: harmonic
    interaction::harmonic_bond_interaction<t_simulation>
      *the_hbond_interaction =
      new interaction::harmonic_bond_interaction<t_simulation>;

    topo >> *the_hbond_interaction;
    bond_param = &the_hbond_interaction->parameter();

    the_forcefield->push_back(the_hbond_interaction);
  }

  if (do_angle){
    // angles
    interaction::angle_interaction<t_simulation>
      *the_angle_interaction = 
      new interaction::angle_interaction<t_simulation>;
    
    topo >> *the_angle_interaction;
  
    the_forcefield->push_back(the_angle_interaction);
  }
  
  if (do_improper){
    // improper dihedrals
    interaction::Improper_dihedral_interaction<t_simulation>
      *the_improper_interaction = 
      new interaction::Improper_dihedral_interaction<t_simulation>;

    topo >> *the_improper_interaction;
    
    the_forcefield->push_back(the_improper_interaction);
  }
  
  if (do_dihedral){
    // dihedrals
    interaction::Dihedral_interaction<t_simulation>
      *the_dihedral_interaction =
      new interaction::Dihedral_interaction<t_simulation>;
    
    topo >> *the_dihedral_interaction;
    
    the_forcefield->push_back(the_dihedral_interaction);
  }
  
  if (do_nonbonded){
    // nonbonded (with virial)
    interaction::Nonbonded_Virial_Interaction<t_simulation, t_pairlist>
      *the_nonbonded_interaction =
      new interaction::Nonbonded_Virial_Interaction<t_simulation, t_pairlist>;
    
    topo >> *the_nonbonded_interaction;
    
    the_forcefield->push_back(the_nonbonded_interaction);
  }

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
      the_topology.solute().
	add_bond_length_constraints(1.0,
				    the_topology.mass(),
				    *bond_param);
      break;
    case 3: 
      std::cout << "SHAKE all bonds" << std::endl;
      the_topology.solute().
	add_bond_length_constraints(*bond_param);
      break;
    default:
      std::cout << "wrong ntc" << std::endl;
  }

  if (shake_param_interaction) delete shake_param_interaction;
  
  return *the_forcefield;
}

