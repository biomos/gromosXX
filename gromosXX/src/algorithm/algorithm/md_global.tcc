/**
 * @file md_global.tcc
 * global functions to get an md simulation started.
 */

#undef MODULE
#undef SUBMODULE
#define MODULE algorithm
#define SUBMODULE algorithm

#include "../../debug.h"

/**
 * perform an MD simulation.
 */
int algorithm::do_md(io::Argument &args)
{
  // decide which code options to use
  // we need the input file!

  io::InInput input;
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

  

  // GRID PAIRLIST?
  bool do_grid;
  int nsnb;
  double rcutp, rcutl, size;
  input.read_PLIST(do_grid, nsnb, rcutp, rcutl, size);

  //==================================================
  // create the algorithm
  //==================================================

  if (do_grid){
    return do_md_grid<true>(args, input);
  }
  else{
    return do_md_grid<false>(args, input);
  }
  
  return 10;

}

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


template<bool do_grid>
int algorithm::do_md_grid(io::Argument &args, io::InInput &input)
{
  // INTERACTION SCALING
  int ntg, nlam;
  double rlam, dlamt;
  bool do_scaled;
  
  input.read_PERTURB(ntg, rlam, dlamt, nlam, do_scaled);

  if (do_scaled){
    return do_md_scaled<do_grid, true>(args, input);
  }
  else{
    return do_md_scaled<do_grid, false>(args, input);
  }
  
}

template<bool do_grid, bool do_scaled>
int algorithm::do_md_scaled(io::Argument &args, io::InInput &input)
{

  // VIRIAL?
  bool calc;
  int ntp;
  double comp, tau;
  math::Matrix pres0;
  interaction::virial_enum do_vir;
  
  input.read_PCOUPLE(calc, ntp, pres0, comp, tau, do_vir);

  switch (do_vir){
    case interaction::no_virial:
      {
	return do_md_virial<do_grid, do_scaled, interaction::no_virial>
	  (args, input);
      }
    case interaction::atomic_virial:
      {
	return do_md_virial<do_grid, do_scaled, interaction::atomic_virial>
	  (args, input);
      }
    case interaction::molecular_virial:
      {
	return do_md_virial<do_grid, do_scaled, interaction::molecular_virial>
	  (args, input);
      }
    default:
      {
	io::messages.add("Wrong virial specified", "md_global",
			 io::message::error);
	return 2;
      }
  }
  return 2;
}

template<bool do_grid, bool do_scaled, interaction::virial_enum do_virial>
int algorithm::do_md_virial(io::Argument &args, io::InInput &input)
{
  // PERTURBATION?
  int ntg, nlam;
  double rlam, dlamt;
  bool scaling;
  
  input.read_PERTURB(ntg, rlam, dlamt, nlam, scaling);
  bool do_perturb = false;
  if (ntg) do_perturb = true;
  if (do_perturb && (args.count("pert") != 1)){
    io::messages.add("PERTURB requested from input file but no "
		     "perturbation topology specified (@pert)!",
		     "md_global::do_md",
		     io::message::error);
  }

  if (do_perturb){
    algorithm::Perturbation_MD<
      algorithm::perturbed_MD_spec<do_virial>,
      algorithm::Interaction_spec<
      typename algorithm::perturbed_MD_spec<do_virial>::simulation_type,
      // perturbation
      true,
      // virial
      do_virial,
      // atomic cutoff
      false,
      // scaling
      do_scaled,
      // bekker
      do_grid
      >
      > 
      the_MD;
    
    return the_MD.do_md(args, input);
  }
  else{
    algorithm::MD<
      algorithm::MD_spec<do_virial>,
      algorithm::Interaction_spec<
      typename algorithm::MD_spec<do_virial>::simulation_type,
      // perturbation
      false,
      // virial
      do_virial,
      // atomic cutoff
      false,
      // scaling
      do_scaled,
      // bekker
      do_grid
      >
      > 
      the_MD;
    
    return the_MD.do_md(args, input);
  }

  return 2;

}
