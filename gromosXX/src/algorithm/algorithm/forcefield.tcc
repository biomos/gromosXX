/**
 * @file forcefield.tcc
 * global functions to create a gromos96 forcefield.
 */

#undef MODULE
#undef SUBMODULE
#define MODULE algorithm
#define SUBMODULE algorithm

#include "../../debug.h"

/**
 * create a Gromos96 like forcefield.
 */
template<typename t_simulation, typename t_interaction_spec>
void algorithm::G96_Forcefield(interaction::Forcefield<
			       t_simulation, t_interaction_spec> &ff,
			       t_simulation & sim,
			       io::InTopology &topo,
			       io::InInput &input,
			       io::Argument &args)
{
  int do_bond, do_angle, do_dihedral, do_improper, do_nonbonded;
  input.read_FORCE(do_bond, do_angle, do_improper,
		   do_dihedral, do_nonbonded);

  int bond_term, angle_term;
  bool have_DIRK =
    input.read_FORCEFIELD(bond_term, angle_term);

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
      typename t_interaction_spec::quartic_bond_interaction_type *
	qbond_interaction =
	new typename t_interaction_spec::quartic_bond_interaction_type;      
    
      topo >> *qbond_interaction;
      ff.push_back(qbond_interaction);
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
      typename t_interaction_spec::harmonic_bond_interaction_type
	*the_hbond_interaction =
	new typename t_interaction_spec::harmonic_bond_interaction_type;

      topo >> *the_hbond_interaction;
      ff.push_back(the_hbond_interaction);

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

    typename t_interaction_spec::angle_interaction_type
      *the_angle_interaction = 
      new     typename t_interaction_spec::angle_interaction_type;
    
    topo >> *the_angle_interaction;
 
    ff.push_back(the_angle_interaction);
  }
  
  if (do_improper){
    // improper dihedrals
    typename t_interaction_spec::improper_interaction_type
      *the_improper_interaction = 
      new typename t_interaction_spec::improper_interaction_type;

    topo >> *the_improper_interaction;
    
    ff.push_back(the_improper_interaction);
  }
  
  if (do_dihedral){
    // dihedrals
    typename t_interaction_spec::dihedral_interaction_type
      *the_dihedral_interaction =
      new typename t_interaction_spec::dihedral_interaction_type;
    
    topo >> *the_dihedral_interaction;
    
    ff.push_back(the_dihedral_interaction);
  }
  


  if (do_nonbonded){

    typename t_interaction_spec::nonbonded_interaction_type
      *the_nonbonded_interaction =
      new typename t_interaction_spec::nonbonded_interaction_type(sim);

    topo >> *the_nonbonded_interaction;
      
    DEBUG(10, "md (create forcefield): nonbonded with pressure read in");

    ff.push_back(the_nonbonded_interaction);
  }

}

/**
 * create a perturbed Gromos96 like forcefield.
 */
template<typename t_simulation, typename t_interaction_spec>
void algorithm::Perturbed_G96_Forcefield(interaction::Forcefield<
					 t_simulation, t_interaction_spec> &ff,
					 t_simulation & sim,
					 io::InTopology &topo,
					 io::InInput &input,
					 io::Argument &args)
{
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
      typename t_interaction_spec::quartic_bond_interaction_type * 
	the_qbond_interaction =
	new typename t_interaction_spec::quartic_bond_interaction_type;
    
      topo >> *the_qbond_interaction;
      ff.push_back(the_qbond_interaction); 
      
      typename t_interaction_spec::perturbed_quartic_bond_interaction_type *
	the_perturbed_qbond_interaction = 
	new typename t_interaction_spec::perturbed_quartic_bond_interaction_type
	(*the_qbond_interaction);

      ff.push_back(the_perturbed_qbond_interaction);
    
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
      typename t_interaction_spec::harmonic_bond_interaction_type *
	the_hbond_interaction =
	new typename t_interaction_spec::harmonic_bond_interaction_type;

      topo >> *the_hbond_interaction;
      ff.push_back(the_hbond_interaction);

      typename t_interaction_spec::perturbed_harmonic_bond_interaction_type *
	the_perturbed_hbond_interaction = 
	new typename t_interaction_spec::perturbed_harmonic_bond_interaction_type
	(*the_hbond_interaction);
      
      ff.push_back(the_perturbed_hbond_interaction);
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

    typename t_interaction_spec::angle_interaction_type *the_angle_interaction = 
      new typename t_interaction_spec::angle_interaction_type;
        
    topo >> *the_angle_interaction;
 
    ff.push_back(the_angle_interaction);

    typename t_interaction_spec::perturbed_angle_interaction_type *
    the_perturbed_angle_interaction = 
    new typename t_interaction_spec::perturbed_angle_interaction_type
    (*the_angle_interaction);
  
    ff.push_back(the_perturbed_angle_interaction);
  }
  
  
  if (do_improper){
    // improper dihedrals
    DEBUG(7, "\timpropers");
    
    typename t_interaction_spec::improper_interaction_type *
      the_improper_interaction = 
      new typename t_interaction_spec::improper_interaction_type;

    topo >> *the_improper_interaction;
    
    ff.push_back(the_improper_interaction);
    DEBUG(7, "\tperturbed impropers");
    
    typename t_interaction_spec::perturbed_improper_interaction_type *
      the_perturbed_improper_interaction = 
      new typename t_interaction_spec::perturbed_improper_interaction_type
      (*the_improper_interaction);
    
    ff.push_back(the_perturbed_improper_interaction);
    DEBUG(7, "\tdone impropers");
    
  }
  
  if (do_dihedral){
    // dihedrals
    DEBUG(7, "\tdihedrals");
    
    typename t_interaction_spec::dihedral_interaction_type *
      the_dihedral_interaction =
      new typename t_interaction_spec::dihedral_interaction_type;
    
    topo >> *the_dihedral_interaction;
    
    ff.push_back(the_dihedral_interaction);
    DEBUG(7, "\tperturbed dihedrals");
    
    typename t_interaction_spec::perturbed_dihedral_interaction_type *
      the_perturbed_dihedral_interaction = 
      new typename t_interaction_spec::perturbed_dihedral_interaction_type
      (*the_dihedral_interaction);
    
    ff.push_back(the_perturbed_dihedral_interaction);
  }
  
  if (do_nonbonded){

    typename t_interaction_spec::perturbed_nonbonded_interaction_type * 
      the_nonbonded_interaction =
      new typename t_interaction_spec
      ::perturbed_nonbonded_interaction_type(sim);
    
    topo >> *the_nonbonded_interaction;
    
    ff.push_back(the_nonbonded_interaction);
    
  }

}

