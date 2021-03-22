/**
 * @file create_forcefield.cc
 */

#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"

#include "../../interaction/interaction.h"
#include "../../interaction/interaction_types.h"
#include "../../interaction/forcefield/forcefield.h"

#include "../../interaction/molecular_virial_interaction.h"

#include "../../io/ifp.h"

#include "create_forcefield.h"
#include "../../interaction/special/qmmm/qm_storage.h"
#include "../../interaction/special/qmmm/mm_atom.h"
#include "../../interaction/special/qmmm_interaction.h"

#include "../../interaction/bonded/create_bonded.h"
#include "../../interaction/nonbonded/create_nonbonded.h"
#include "../../interaction/special/create_special.h"

#ifdef XXMPI
#include <mpi.h>
#endif

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE forcefield

/**
 * create a Gromos96 (like) forcefield.
 */
int interaction::create_g96_forcefield(interaction::Forcefield & ff,
				       topology::Topology const & topo,
				       simulation::Simulation & sim,
				       io::IFP & it,
				       std::ostream & os,
				       bool quiet)
{
  if (!quiet)
    os << "FORCEFIELD\n";
  
  // the bonded
  DEBUG(8, "creating the bonded terms");
  if (create_g96_bonded(ff, topo, sim.param(), it, os, quiet))
    return 1;

  // the nonbonded
  DEBUG(8, "creating the nonbonded terms");
  if (create_g96_nonbonded(ff, topo, sim, it, os, quiet))
    return 1;

  // correct the virial (if molecular virial is required)
  if (sim.param().pcouple.virial == math::molecular_virial){
    
    Molecular_Virial_Interaction * mvi = new Molecular_Virial_Interaction;
    ff.push_back(mvi);
  }

  // the special
  DEBUG(8, "creating the special terms");
  if(create_special(ff, topo, sim.param(), os, quiet))
    return 1;
  
  DEBUG(8, "creating the QM/MM terms");
  bool add_qmmm = true;
#ifdef XXMPI
  if (sim.mpi && MPI::COMM_WORLD.Get_rank() != 0) {
    add_qmmm = false;
  }
#endif  
  if (sim.param().qmmm.qmmm != simulation::qmmm_off && add_qmmm) {
    interaction::QMMM_Interaction * qmmm = new interaction::QMMM_Interaction;
    ff.push_back(qmmm);
    sim.param().qmmm.interaction = qmmm;
  }

  if (!quiet){
  
    if (sim.param().perturbation.perturbation){
      os << "\n\t" << std::setw(20) << std::left << "perturbation" 
	 << std::setw(30) << "on" << std::right << "\n"
	 << "\t\tlambda         : " << sim.param().perturbation.lambda << "\n"
	 << "\t\texponent       : " << sim.param().perturbation.lambda_exponent << "\n"
	 << "\t\tdlambda        : " << sim.param().perturbation.dlamt << "\n"
	 << "\t\tscaling        : ";
      
      if (sim.param().perturbation.scaling){
	if (sim.param().perturbation.scaled_only)
	  os << "perturbing only scaled interactions\n";
	else
	  os << "on\n";
      }
      else
	os << "off\n";

      if (topo.perturbed_solute().atoms().size() == 0)
	os << "\t\t" << "using unperturbed nonbonded routines as no atoms are perturbed\n";
      else os << "\t\t" << "with " << topo.perturbed_solute().atoms().size() << " perturbed atoms\n";

      if(sim.param().lambdas.individual_lambdas){

	// say what we are doing
	os << "\n\t\tIndividual lambdas according to l_int = ALI * l^4 + BLI * l^3 + CLI * l^2 + DLI * l + ELI\n"
	   << "\t\tfor interactions between different chargegroups NILG1 and NILG2"
	   <<     "\n\t\t" << std::setw(20) << "Interaction         "
	   << std::setw(6) << "NILG1"
	   << std::setw(6) << "NILG2"
	   << std::setw(8) << "ALI"
	   << std::setw(8) << "BLI"
	   << std::setw(8) << "CLI"
	   << std::setw(8) << "DLI"
	   << std::setw(8) << "ELI"
	   << "\n";
	for(int i=0; i < simulation::last_interaction_lambda; i++){
	  for(unsigned int n1=0; n1 < sim.param().lambdas.a[i].size(); n1++){
	    for(unsigned int n2=n1; n2 < sim.param().lambdas.a[i].size(); n2++){
	      if(sim.param().lambdas.a[i][n1][n2] != 0 ||
		 sim.param().lambdas.b[i][n1][n2] != 0 ||
		 sim.param().lambdas.c[i][n1][n2] != 0 ||
		 sim.param().lambdas.d[i][n1][n2] != 1 ||
		 sim.param().lambdas.e[i][n1][n2] != 0){
		  
		os << "\t\t";
		switch(i) {
		  case simulation::bond_lambda :
		    os << "bonds              :";
		    break;
		  case simulation::angle_lambda :
		    os << "angles             :";
		    break;
		  case simulation::dihedral_lambda :
		    os << "dihedrals          :";
		    break;
		  case simulation::improper_lambda :
		    os << "impropers          :";
		    break;
		  case simulation::lj_lambda :
		    os << "vdw interaction    :";
		    break;
		  case simulation::lj_softness_lambda :
		    os << "vdw softness       :";
		    break;
		  case simulation::crf_lambda :
		    os << "crf interaction    :";
		    break;
		  case simulation::crf_softness_lambda :
		    os << "crf softness       :";
		    break;
		  case simulation::disres_lambda :
		    os << "distance restraint :";
		    break;
		  case simulation::dihres_lambda :
		    os << "dihedral restraint :";
		    break;
		  case simulation::mass_lambda :
		    os << "mass scaling       :";
		    break;
		}
		os << std::setw(6) << n1+1 << std::setw(6) << n2+1
		   << std::setw(8) << sim.param().lambdas.a[i][n1][n2]
		   << std::setw(8) << sim.param().lambdas.b[i][n1][n2]
		   << std::setw(8) << sim.param().lambdas.c[i][n1][n2]
		   << std::setw(8) << sim.param().lambdas.d[i][n1][n2]
		   << std::setw(8) << sim.param().lambdas.e[i][n1][n2]
		   << "\n";
	      }
	    }
	  }
	}
	os << "\n\t\tall other interactions / interactions for other energy groups are"
	   << " scaled using overall lambda\n";
      }
      else{
	os << "\t\tall interactions are scaled using overall lambda\n"
	   << "\t\t(no individual lambdas)\n";
      }
    }
    else{
      os << "\t" << std::setw(20) << std::left << "perturbation" 
	 << std::setw(30) << std::left << "off" << std::right << "\n";
    }
    
    os << "END\n";
  }
  
  return 0;

}

