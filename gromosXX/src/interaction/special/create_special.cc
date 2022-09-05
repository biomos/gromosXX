/**
 * @file create_special.cc
 * create the special terms.
 */

#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"
#include "../../interaction/interaction.h"
#include "../../interaction/forcefield/forcefield.h"

// special interactions
#include "../../interaction/interaction_types.h"

#include "../../interaction/special/position_restraint_interaction.h"
#include "../../interaction/special/distance_restraint_interaction.h"
#include "../../interaction/special/distance_field_interaction.h"
#include "../../interaction/special/angle_restraint_interaction.h"
#include "../../interaction/special/dihedral_restraint_interaction.h"
#include "../../interaction/special/dfunct_interaction.h"
#include "../../interaction/special/perturbed_distance_restraint_interaction.h"
#include "../../interaction/special/perturbed_distance_field_interaction.h"
#include "../../interaction/special/eds_distance_restraint_interaction.h"
#include "../../interaction/special/perturbed_angle_restraint_interaction.h"
#include "../../interaction/special/perturbed_dihedral_restraint_interaction.h"
#include "../../interaction/special/jvalue_restraint_interaction.h"
#include "../../util/umbrella_weight.h"
#include "../../interaction/special/xray_restraint_interaction.h"
#include "../../interaction/special/adde_reweighting.h"
#include "../../interaction/special/nemd.h"
#include "../../interaction/special/local_elevation_interaction.h"
#include "../../interaction/special/electric_field_interaction.h"
#include "../../interaction/special/order_parameter_restraint_interaction.h"
#include "../../interaction/special/rdc_restraint_interaction.h"
#include "../../interaction/special/symmetry_restraint_interaction.h"

#include "../../interaction/bonded/dihedral_interaction.h"
#include "../../interaction/bonded/dihedral_new_interaction.h"
#include "../../interaction/special/pscale.h"

#include "../../io/instream.h"
#include "../../io/topology/in_topology.h"

#include "create_special.h"
#include "bs_interaction.h"

int interaction::create_special(interaction::Forcefield & ff,
				topology::Topology const & topo,
				simulation::Parameter const & param,
				std::ostream & os,
				bool quiet)
{
  // if (!quiet)
  // os << "SPECIAL\n";
  if (param.dfunct.dfunct > simulation::dfunct_off) {
    if (!quiet) {
      os << "\tDFUNCT\n";
    }
    interaction::DFunct_Interaction* dfunct = new interaction::DFunct_Interaction();
    ff.push_back(dfunct);
  }

  // Position restraints / constraints
  if (param.posrest.posrest == simulation::posrest_on || 
      param.posrest.posrest == simulation::posrest_bfactor) {

    if(!quiet)
      os <<"\tPosition restraints\n";

    interaction::Position_Restraint_Interaction *pr =
      new interaction::Position_Restraint_Interaction;

    ff.push_back(pr);
    
    if (param.pcouple.virial == math::atomic_virial)
      io::messages.add("Position restraints with atomic virial ill defined",
		       "create_special", io::message::warning);
    else if (param.pcouple.virial == math::molecular_virial)
      io::messages.add("Position restraint forces not added to molecular virial",
		       "create_special", io::message::warning);
  } else if (param.posrest.posrest != simulation::posrest_off &&
             param.posrest.posrest != simulation::posrest_const){
    io::messages.add("Wrong value for position restraints",
		     "create_special", io::message::error);
  }
  // Distance restraints 
  // as this function is called before read_special, we do not know if we 
  // really need the (perturbed) interactions
  // so we generate them, but need to check inside if we really do something
  if (abs(param.distanceres.distanceres) == 1 || 
      abs(param.distanceres.distanceres) == 2){

    if(!quiet)
      os <<"\tDistance restraints\n";

    interaction::Distance_Restraint_Interaction *dr =
      new interaction::Distance_Restraint_Interaction();

    ff.push_back(dr);
    
    if(param.perturbation.perturbation){
      if(!quiet)
	os <<"\tPerturbed distance restraints\n";
      
      interaction::Perturbed_Distance_Restraint_Interaction *pdr =
	new interaction::Perturbed_Distance_Restraint_Interaction;

      ff.push_back(pdr); 
    }
    if(param.eds.eds){
      if(!quiet)
        os << "\tEDS perturbed distance restraints\n";
      interaction::Eds_Distance_Restraint_Interaction *edr =
        new interaction::Eds_Distance_Restraint_Interaction;
      
      ff.push_back(edr);
    }
  }
  
  // Distance field restraint
  // as this function is called before read_special, we do not know if we 
  // really need the (perturbed) interactions
  // so we generate them, but need to check inside if we really do something
  if (param.distancefield.distancefield == 1){
    if(!quiet)
      os << "\tDistancefield restraint\n";
    interaction::Distance_Field_Interaction *df = 
      new interaction::Distance_Field_Interaction();
    
    ff.push_back(df);
    DEBUG(10, "Distancefield interaction pushed back");
    
    if (param.perturbation.perturbation){
      if(!quiet)
	os <<"\tPerturbed distancefield restraints \n";
      
      interaction::Perturbed_Distance_Field_Interaction *pdf =
	new interaction::Perturbed_Distance_Field_Interaction;
      
      ff.push_back(pdf);
    }
  }
  
  // Angle restraints 
  if (param.angrest.angrest == simulation::angle_restr_inst ||
      param.angrest.angrest == simulation::angle_restr_inst_weighted){

    if(!quiet)
      os <<"\tAngle restraints\n";

    interaction::Angle_Restraint_Interaction *dr =
      new interaction::Angle_Restraint_Interaction();

    ff.push_back(dr);
    
    if(param.perturbation.perturbation){
      if(!quiet)
	os <<"\tPerturbed angle restraints\n";
      
      interaction::Perturbed_Angle_Restraint_Interaction *pdr =
	new interaction::Perturbed_Angle_Restraint_Interaction;
      
      ff.push_back(pdr); 
    }
  }
  
  // Dihedral restraints 
  if (param.dihrest.dihrest == simulation::dihedral_restr_inst ||
      param.dihrest.dihrest == simulation::dihedral_restr_inst_weighted){

    if(!quiet)
      os <<"\tDihedral restraints\n";

    interaction::Dihedral_Restraint_Interaction *dr =
      new interaction::Dihedral_Restraint_Interaction();

    ff.push_back(dr);
    
    if(param.perturbation.perturbation){
      if(!quiet)
	os <<"\tPerturbed dihedral restraints\n";
      
      interaction::Perturbed_Dihedral_Restraint_Interaction *pdr =
	new interaction::Perturbed_Dihedral_Restraint_Interaction;
      
      ff.push_back(pdr); 
    }
  }

  // J-Value restraints
  if (param.jvalue.mode != simulation::jvalue_restr_off){
    if(!quiet){
      os << "\tJ-Value restraints (";
      switch(param.jvalue.mode){
	case simulation::jvalue_restr_inst :
	  os << "instantaneous";
	  break;
	case simulation::jvalue_restr_av :
	  os << "time averaged";
	  break;
 	case simulation::jvalue_restr_inst_weighted :
	  os << "instantaneous, weighted";
	  break;
	case simulation::jvalue_restr_av_weighted :
	  os << "time averaged, weighted";
	  break;
	case simulation::jvalue_restr_biq_weighted :
	  os << "biquadratic, weighted";
	  break;
	default:
	  os << "unknown mode!";
	  break;
      }
      os << ")\n";
    }

    interaction::Jvalue_Restraint_Interaction *jr =
      new interaction::Jvalue_Restraint_Interaction;
    
    ff.push_back(jr);
  }

  // Xray restraints
  if (param.xrayrest.xrayrest != simulation::xrayrest_off) {
    if (!quiet) {
      os << "\tXray restraints \n";
      /*switch (param.xry) {
        case simulation::restr_inst :
                  os << "instantaneous";
          break;
        case simulation::restr_av :
                  os << "time averaged";
          break;
        case simulation::restr_biq :
                  os << "biquadratic";
          break;
        default:
          os << "unknown mode!";
          break;
      }*/
    }

    interaction::Xray_Restraint_Interaction *xr =
            new interaction::Xray_Restraint_Interaction;

    ff.push_back(xr);
  }

  // order parameter restraints
  if (param.orderparamrest.orderparamrest != simulation::oparam_restr_off){
    if(!quiet){
      os << "\tOrder-parameter restraints (";
      switch(param.orderparamrest.orderparamrest){
	case simulation::oparam_restr_av :
	  os << "time averaged";
	  break;
	case simulation::oparam_restr_av_weighted :
	  os << "time averaged, weighted";
	  break;
	default:
	  os << "unknown mode!";
	  break;
      }
      os << ")\n";
    }

    interaction::Order_Parameter_Restraint_Interaction *opr =
      new interaction::Order_Parameter_Restraint_Interaction;

    ff.push_back(opr);
  }

  // RDC restraints
  if (param.rdc.mode != simulation::rdc_restr_off){
    if(!quiet){
      os << "\tRDC restraints (";
      switch(param.rdc.mode){
	    case simulation::rdc_restr_inst :
	      os << "instantaneous";
	      break;
 	    case simulation::rdc_restr_inst_weighted :
	      os << "instantaneous, weighted";
	      break;
	    case simulation::rdc_restr_av :
	      os << "time averaged";
	      break;
	    case simulation::rdc_restr_av_weighted :
	      os << "time averaged, weighted";
	      break;
	    case simulation::rdc_restr_biq :
	      os << "biquadratic";
	      break;
	    case simulation::rdc_restr_biq_weighted :
	      os << "biquadratic, weighted";
	      break;
	    default:
	      os << "unknown mode!";
	      break;
      }
      os << ")\n";
    }

    interaction::RDC_Restraint_Interaction *rdcr =
      new interaction::RDC_Restraint_Interaction;

    ff.push_back(rdcr);
  }

  // symmetry restraints
  if (param.symrest.symrest != simulation::xray_symrest_off){
    if(!quiet){
      os << "\tSymmetry restraints (";
      switch(param.symrest.symrest){
	case simulation::xray_symrest_ind :
	  os << "harmonic";
	  break;
	case simulation::xray_symrest_constr :
	  os << "constrained";
	  break;
	default:
	  os << "unknown mode!";
	  break;
      }
      os << ")\n";
    }

    interaction::Symmetry_Restraint_Interaction *symri =
      new interaction::Symmetry_Restraint_Interaction;

    ff.push_back(symri);
  }

  if (param.localelev.localelev != simulation::localelev_off) {
    if (!quiet) {
      os << "\tlocal elevation \n";
    }
    interaction::Local_Elevation_Interaction *le =
            new interaction::Local_Elevation_Interaction();
    ff.push_back(le);
  }
  
  if (param.bsleus.bsleus != simulation::bsleus_off){
    if (!quiet) {
      os << "\tB&S-LEUS\n";
    }
    interaction::BS_Interaction *bsleus_ia =
            new interaction::BS_Interaction();
    ff.push_back(bsleus_ia);
  }

  if (param.electric.electric != simulation::electric_off) {
    if (!quiet) {
      os << "\telectric field \n";
    }
    interaction::Electric_Field_Interaction *efield =
            new interaction::Electric_Field_Interaction();
    ff.push_back(efield);
  }



  // Periodic Scaling
  // right now this only works if no dihedral angles with j-value restraints on are perturbed...
  if (param.pscale.jrest){

    if(!quiet){
      os << "\tscaling based on J-Value restraints\n";
    }
    
    interaction::Periodic_Scaling * ps = 
      new interaction::Periodic_Scaling(ff, param);

    ff.push_back(ps);
  }

  if(param.addecouple.adgr>0){
    if(!quiet){
      os << "\tadding adiabatic decoupling reweighting\n";
    }
    interaction::Adde_Reweighting * ad =
            new interaction::Adde_Reweighting();
    ff.push_back(ad);
  }
  
  if (param.nemd.nemd != simulation::nemd_off) {
    if (!quiet) {
      os << "\tnemd \n";
    }
    interaction::NEMD_Interaction *nemd =
            new interaction::NEMD_Interaction();
    ff.push_back(nemd);
  }
  
  return 0;
}
