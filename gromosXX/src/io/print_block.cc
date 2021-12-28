/**
 * @file print_block.cc
 * routines to print out the various blocks.
 */

#include <stdheader.h>
#include <fstream>


#include "../algorithm/algorithm.h"
#include "../topology/topology.h"
#include "../simulation/simulation.h"
#include "../configuration/configuration.h"
#include "../interaction/interaction.h"
#include "../interaction/interaction_types.h"

#include "../io/argument.h"
#include "../io/blockinput.h"
#include "../io/instream.h"
#include "../io/configuration/inframe.h"
#include "../io/configuration/in_configuration.h"
#include "../io/topology/in_topology.h"
#include "../io/topology/in_perturbation.h"
#include "../io/parameter/in_parameter.h"

#include "../algorithm/algorithm/algorithm_sequence.h"
#include "../algorithm/create_md_sequence.h"

#include "../interaction/forcefield/forcefield.h"

#include "../math/volume.h"

#include "print_block.h"

namespace io
{

  /* 
   * Print the DOF COUPLING table of MULTIBATH block.
   */
  void print_MULTIBATH_COUPLING(std::ostream &os,
				simulation::Multibath const &bath)
  {
    os << "MULTIBATHCOUPLING\n";
    os.precision(2);
    os.setf(std::ios_base::fixed, std::ios_base::floatfield);
    
    os << std::setw(12) << "LAST-ATOM"
       << std::setw(12) << "LAST-MOL"
       << std::setw(12) << "COM-BATH"
       << std::setw(12) << "IR-BATH"
       << "\n";
    std::vector<simulation::bath_index_struct>::const_iterator
      it = bath.bath_index().begin(),
      to = bath.bath_index().end();
    
    for(; it!=to; ++it){
      os << std::setw(12) << it->last_atom + 1
	 << std::setw(12) << it->last_temperature_group + 1
	 << std::setw(12) << it->com_bath + 1
	 << std::setw(12) << it->ir_bath + 1
	 << "\n";
    }
    
    os << "END\n";
    
  }

  /*
   * Print DEGREESOFFREEDOM block.
   */
  void print_DEGREESOFFREEDOM(std::ostream &os,
			      simulation::Multibath const &bath)
  {
    os << "DEGREES OF FREEDOM\n";
  
    os.precision(2);
    os.setf(std::ios_base::fixed, std::ios_base::floatfield);
    
    os << std::setw(10) << "BATH";
    os << std::setw( 8) << "TEMP0"
       << std::setw( 8) << "TAU"
       << std::setw(10) << "DOF"
       << std::setw(10) << "MOL-DOF"
       << std::setw(10) << "IR-DOF"
       << std::setw(10) << "SOLUC"
       << std::setw(10) << "SOLVC"
       << "\n";
  
    double avg_temp0 = 0, avg_tau = 0, sum_dof = 0, sum_soluc = 0,
      sum_solvc = 0, tau_dof = 0,
      sum_ir_dof = 0, sum_com_dof = 0;

    std::vector<simulation::bath_struct>::const_iterator
      it = bath.begin(),
      to = bath.end();
  
    for(unsigned int i=0; it != to; ++it, ++i){
      
      os << std::setw(10) << i
	 << std::setw( 8) << it->temperature
	 << std::setw( 8) << it->tau
	 << std::setw(10) << it->dof
	 << std::setw(10) << it->com_dof
	 << std::setw(10) << it->ir_dof
	 << std::setw(10) << it->solute_constr_dof
	 << std::setw(10) << it->solvent_constr_dof
	 << "\n";

      if (it->tau != -1){
	tau_dof += it->dof;
	avg_tau += it->tau * it->dof;
	avg_temp0 += it->temperature * it->dof;
      }
      sum_dof += it->dof;
      sum_ir_dof += it->ir_dof;
      sum_com_dof += it->com_dof;

      sum_soluc += it->solute_constr_dof;
      sum_solvc += it->solvent_constr_dof;

    }

    os << "    -------------------------------------------"
       << "-------------------------------\n";

    os << std::setw(10) << "Total";
    if (tau_dof)
      os << std::setw( 8) << avg_temp0 / tau_dof
	 << std::setw( 8) << avg_tau / tau_dof;
    else
      os << std::setw( 8) << "-"
	 << std::setw( 8) << "-";
    
    os << std::setw(10) << sum_dof
       << std::setw(10) << sum_com_dof
       << std::setw(10) << sum_ir_dof
       << std::setw(10) << sum_soluc
       << std::setw(10) << sum_solvc
       << "\n";
    
    os << "END\n";
    
  }
  
  /* 
   * Print the MULTIBATH block.
   */
  void print_MULTIBATH(std::ostream &os,
		       simulation::Multibath const &bath,
		       configuration::Energy const &energy,
		       std::string title)
  {
    os << title << "\n";
  
    os.precision(2);
    os.setf(std::ios_base::fixed, std::ios_base::floatfield);
  
    os << std::setw(10) << "BATH"
       << std::setw(13) << "EKIN"
       << std::setw(13) << "EKIN-MOL-TR"
       << std::setw(13) << "EKIN-MOL-IR"
       << std::setw(10) << "T"
       << std::setw(10) << "T-MOL-TR"
       << std::setw(10) << "T-MOL-IR"
       << std::setw(12) << "SCALE"
       << "\n";
  
    double avg_temp0 = 0, avg_tau = 0, sum_dof = 0, sum_soluc = 0,
      sum_solvc = 0, sum_ekin = 0, tau_dof = 0,
      sum_com_ekin = 0, sum_ir_ekin = 0,
      sum_ir_dof = 0, sum_com_dof = 0,
      avg_scale = 0;

    std::vector<simulation::bath_struct>::const_iterator
      it = bath.begin(),
      to = bath.end();
  
    for(unsigned int i=0; it != to; ++it, ++i){
      
      const double e_kin = energy.kinetic_energy[i];
      const double e_kin_com = energy.com_kinetic_energy[i];
      const double e_kin_ir = energy.ir_kinetic_energy[i];

      os << std::setw(10) << i
	 << std::setw(13) << std::setprecision(4) << std::scientific 
	 << e_kin
	 << std::setw(13) 
	 << e_kin_com
	 << std::setw(13) 
	 << e_kin_ir
	 << std::setprecision(2) << std::fixed;
      if (it->dof == 0){
	os << std::setw(10) << 0;
      }
      else{
	os << std::setw(10) 
	   << 2 * e_kin / (math::k_Boltzmann * it->dof);
      }
      if (it->com_dof == 0){
	os << std::setw(10) << "-";
      }
      else{
	os << std::setw(10) 
	   << 2 * e_kin_com / 
	  (math::k_Boltzmann * it->com_dof);
      }
      if (it->ir_dof == 0) {
	os << std::setw(10) << "-";
      }
      else{
	os << std::setw(10) 
	   << 2 * e_kin_ir / 
	  (math::k_Boltzmann * it->ir_dof);
      }
      if (it->tau != -1){
	os << std::setw(12)
	   << std::setprecision(7)
	   << it->scale;
      }
      else{
	os << std::setw(12)
	   << "-";
      }
	  
      if (it->tau != -1){
	tau_dof += it->dof;
	avg_temp0 += it->temperature * it->dof;
	avg_tau += it->tau * it->dof;
	avg_scale += it->scale * it->dof;
      }

      sum_dof += it->dof;
      sum_ir_dof += it->ir_dof;
      sum_com_dof += it->com_dof;

      sum_soluc += it->solute_constr_dof;
      sum_solvc += it->solvent_constr_dof;
      sum_ekin += e_kin;
      sum_com_ekin += e_kin_com;
      sum_ir_ekin += e_kin_ir;

      os << "\n";

    }

    os << "    ---------------------------------------------"
       << "------------------------------------------\n";

    os << std::setw(10) << "T_avg"
       << std::setw(13) << std::setprecision(4) << std::scientific
       << sum_ekin
       << std::setw(13) << sum_com_ekin
       << std::setw(13) << sum_ir_ekin
	   << std::setprecision(2) << std::fixed;
	if (sum_dof)
       os << std::setw(10) << 2 * sum_ekin / (math::k_Boltzmann * sum_dof);
	else os << std::setw(10) << "-";
	if (sum_com_dof)
		os << std::setw(10) << 2 * sum_com_ekin / (math::k_Boltzmann * sum_com_dof);
	else os << std::setw(10) << "-";
	if (sum_ir_dof)
		os << std::setw(10) << 2 * sum_ir_ekin / (math::k_Boltzmann * sum_ir_dof);
	else os << std::setw(10) << "-";
    if (tau_dof)
      os << std::setw(12) << std::setprecision(7) << avg_scale / tau_dof;
    else
      os << std::setw(12) << std::right << "-";

    os << "\n";
    
    os << "END\n";

  }

  /*
   * Print the PCOUPLE block.
   */
  void print_PCOUPLE(std::ostream &os, bool calc, 
		     math::pressure_scale_enum scale,
		     math::Matrix pres0, double comp, 
		     double tau, math::virial_enum vir,
                     int x_semi, int y_semi, int z_semi)
  {
    os << "PRESSURE SCALING\n";
    os.precision(5);
    os.setf(std::ios_base::fixed, std::ios_base::floatfield);

    os << std::setw(12) << "COUPLE"
       << std::setw(12) << "SCALE"
       << std::setw(12) << "COMP"
       << std::setw(12) << "TAU"
       << std::setw(12) << "VIRIAL"
       << "\n";
    
    if (calc && scale != math::pcouple_off)
      os << std::setw(12) << "scale";
    else if (calc)
      os << std::setw(12) << "calc";
    else
      os << std::setw(12) << "none";

    switch(scale){
      case math::pcouple_off:
	os << std::setw(12) << "off";
	break;
      case math::pcouple_isotropic:
	os << std::setw(12) << "iso";
	break;
      case math::pcouple_anisotropic:
	os << std::setw(12) << "aniso";
	break;
      case math::pcouple_full_anisotropic :
                os << std::setw(12) << "full";
        break;
      case math::pcouple_semi_anisotropic :
                os << std::setw(12) << "semianiso";
        break;
      default:
        os << std::setw(12) << "unknown";
    }
    
    os << std::setw(12) << comp
       << std::setw(12) << tau;
    
    switch(vir){
      case math::no_virial:
	os << std::setw(12) << "none";
	break;
      case math::atomic_virial:
	os << std::setw(12) << "atomic";
	break;
      case math::molecular_virial:
	os << std::setw(12) << "molecular";
	break;
      default:
	os << std::setw(12) << "unknown";
    }

    if(scale == math::pcouple_semi_anisotropic){
      os << "\n" << std::setw(39) << "SEMIANISOTROPIC COUPLINGS(X, Y, Z)" << "\n";
      os << std::setw(12) << x_semi
         << std::setw(12) << y_semi
         << std::setw(12) << z_semi;
    }
    else{
      os << "\n" << std::setw(39) << "SEMIANISOTROPIC COUPLINGS(X, Y, Z)" << "\n";
      os << std::setw(6) << x_semi << " (ignored)"
         << std::setw(6) << y_semi << " (ignored)"
         << std::setw(6) << z_semi << " (ignored)";
    }



    os << "\n" << std::setw(23) << "REFERENCE PRESSURE" << "\n";
    
    for(int i=0; i<3; ++i){
      for(int j=0; j<3; ++j){
	os << std::setw(12) << pres0(i,j);
      }
      os << "\n";
    }
    
    os << "END\n";

  }

  /*
   * Print the PRESSURE block.
   */
  void print_PRESSURE(std::ostream &os,
		      configuration::Configuration const & conf)
  {
    os << "PRESSURE\n";
    os.precision(5);
    os.setf(std::ios_base::scientific, std::ios_base::floatfield);
    
    os << "\tmolecular kinetic energy:\n\t";
    for(int i=0; i<3; ++i){
      for(int j=0; j<3; ++j)
	os << std::setw(15) << conf.old().kinetic_energy_tensor(i,j);
      os << "\n\t";
    }
    
    os << "\n\tvirial\n\t";
    // the virial is stored internally as just the outer product of positions and forces
    // so without the -0.5 prefactor.

    for(int i=0; i<3; ++i){
      for(int j=0; j<3; ++j)
	os << std::setw(15) << -0.5 * conf.old().virial_tensor(i,j);
      os << "\n\t";
    }

    os << "\n\tpressure tensor\n\t";
    for(int i=0; i<3; ++i){
      for(int j=0; j<3; ++j)
	os << std::setw(15) << conf.old().pressure_tensor(i,j);
      os << "\n\t";
    }
    os << "\n\tpressure: "
       << std::setw(15)
       << (conf.old().pressure_tensor(0,0) + 
	   conf.old().pressure_tensor(1,1) +
	   conf.old().pressure_tensor(2,2)) / 3 
       << "\n";

    os << "\tvolume:   "
       << std::setw(15)
       << math::volume(conf.old().box, conf.boundary_type)
       << "\n";
    
    os << "\nEND\n";

    os.precision(5);
    os.setf(std::ios_base::fixed, std::ios_base::floatfield);
    
  }

  /*
   * Print the ENERGY block.
   */
  void print_ENERGY(std::ostream &os,
		    configuration::Energy const &e,
		    std::vector<unsigned int> const &energy_groups,
		    std::string const title,
		    std::string const type)
  {

    unsigned int numenergygroups = unsigned(e.bond_energy.size());
 
    std::vector<std::string> energroup;
   
    int b=1;
    
    for(unsigned int i=0; i<numenergygroups; i++){

      std::ostringstream ostring;
      ostring << b << "-" << energy_groups[i]+1;
      energroup.push_back(ostring.str());
      b = energy_groups[i]+2;
    }
        
    os << title << "\n";

    os.precision(4);
    os.setf(std::ios_base::scientific, std::ios_base::floatfield);
    os << type << "Total                : " << std::setw(12) << e.total << "\n";
    os << type << "Kinetic              : " << std::setw(21) << e.kinetic_total << "\n";
    os << type << "Potential+Special    : " << std::setw(21) << e.potential_total + e.special_total << "\n";
    os << type << "Potential            : " << std::setw(21) << e.potential_total << "\n";
    os << type << "Covalent             : " << std::setw(30) << e.bonded_total << "\n";
    os << type << "Bonds                : " << std::setw(39) << e.bond_total << "\n";
    os << type << "Angles               : " << std::setw(39) << e.angle_total << "\n";
    os << type << "Improper             : " << std::setw(39) << e.improper_total << "\n";
    os << type << "Dihedral             : " << std::setw(39) << e.dihedral_total << "\n";
    os << type << "Crossdihedral        : " << std::setw(39) << e.crossdihedral_total << "\n";
    os << type << "Non-bonded           : " << std::setw(30) << e.nonbonded_total  << "\n";
    os << type << "Vdw                  : " << std::setw(39) << e.lj_total << "\n";
    os << type << "El (RF)              : " << std::setw(39) << e.crf_total  << "\n";
    os << type << "El (LS)              : " << std::setw(39) << e.ls_total << "\n";
    os << type << "El (pair)            : " << std::setw(48) << e.ls_pair_total << "\n";
    os << type << "El (real-space)      : " << std::setw(57) << e.ls_realspace_total << "\n";
    os << type << "El (k-space)         : " << std::setw(57) << e.ls_kspace_total << "\n";
    os << type << "El (A term)          : " << std::setw(57) << e.ls_a_term_total << "\n";
    os << type << "El (lattice sum self): " << std::setw(48) << e.ls_self_total << "\n";
    os << type << "El (surface term)    : " << std::setw(48) << e.ls_surface_total << "\n";
    os << type << "Polarisation self    : " << std::setw(39) << e.self_total << "\n";
    os << type << "QM                   : " << std::setw(30) << e.qm_total << "\n";
    os << type << "Special              : " << std::setw(21) << e.special_total << "\n";
    os << type << "SASA                 : " << std::setw(30) << e.sasa_total << "\n";
    os << type << "SASA Volume          : " << std::setw(30) << e.sasa_volume_total << "\n";
    os << type << "Constraints          : " << std::setw(30) << e.constraints_total << "\n";
    os << type << "Distanceres          : " << std::setw(30) << e.distanceres_total << "\n";
    os << type << "Disfieldres          : " << std::setw(30) << e.disfieldres_total << "\n";
    os << type << "Angrest              : " << std::setw(30) << e.angrest_total << "\n";
    os << type << "Dihrest              : " << std::setw(30) << e.dihrest_total << "\n";
    os << type << "Posrest              : " << std::setw(30) << e.posrest_total << "\n";
    os << type << "Jrest                : " << std::setw(30) << e.jvalue_total << "\n";
    os << type << "X-ray restraints     : " << std::setw(30) << e.xray_total << "\n";
    os << type << "Local elevation      : " << std::setw(30) << e.leus_total << "\n";
    os << type << "Order-parameter rest.: " << std::setw(30) << e.oparam_total << "\n";
    os << type << "RDCrest              : " << std::setw(30) << e.rdc_total << "\n";    
    os << type << "Symmetry restraints  : " << std::setw(30) << e.symrest_total << "\n";
    os << type << "EDS reference        : " << std::setw(30) << e.eds_vr << "\n";
    os << type << "Entropy              : " << std::setw(30) << e.entropy_term << "\n";
    os << type << "NN validation        : " << std::setw(30) << e.nn_valid << "\n";
    os << "\n";

    os << std::setw(20) << "COV";
    
    for(unsigned int i=0; i < numenergygroups; i++) os << std::setw(12) << energroup[i];
    os << "\n" << std::setw(20) << type + "bonds";
    for(unsigned int i=0; i < numenergygroups; i++) os << std::setw(12) << e.bond_energy[i];
    os << "\n" << std::setw(20) << type + "angles";
    for(unsigned int i=0; i < numenergygroups; i++) os << std::setw(12) << e.angle_energy[i];
    os << "\n" << std::setw(20) << type + "impropers";
    for(unsigned int i=0; i < numenergygroups; i++) os << std::setw(12) << e.improper_energy[i];
    os << "\n" << std::setw(20) << type + "dihedrals";
    for(unsigned int i=0; i < numenergygroups; i++) os << std::setw(12) << e.dihedral_energy[i];
    os << "\n" << std::setw(20) << type + "crossdihedrals";
    for(unsigned int i=0; i < numenergygroups; i++) os << std::setw(12) << e.crossdihedral_energy[i];

    os << "\n" << "\n";
    os << std::setw(20) << type + "VDW";
    
    for(unsigned int i=0; i < numenergygroups; i++) os << std::setw(12) << energroup[i];
    os << "\n";
    for(unsigned int j=0; j < numenergygroups; j++) {
      os << std::setw(20) << energroup[j];
      for(unsigned int i=0; i<j; i++) os << std::setw(12) << " ";
      for(unsigned int i=j; i < numenergygroups; i++){
	// now in calculate_totals
	// if(i==j)
	os << std::setw(12) << e.lj_energy[i][j];
	// else 
	// os << std::setw(12) << e.lj_energy[i][j] + e.lj_energy[j][i] ;
      }
      os << "\n";
    }
    os << "\n" << std::setw(20) << type + "CRF";
    
    for(unsigned int i=0; i < numenergygroups; i++) os << std::setw(12) << energroup[i];
    os << "\n";
    for(unsigned int j=0; j < numenergygroups; j++) {
      os << std::setw(20) << energroup[j];
      for(unsigned int i=0; i<j; i++) os << std::setw(12) << " ";
      for(unsigned int i=j; i < numenergygroups; i++){
	// now in calculate totals
	// if(i==j)
	os << std::setw(12) << e.crf_energy[i][j];
	// else
	// os << std::setw(12) << e.crf_energy[i][j] +  e.crf_energy[j][i];
      }
      os << "\n";
    }

    // ANITA 
    // print A_e_lj, B_e_lj, A_e_crf, B_e_crf for the lambda at which
    // we are simulating
/*    double lambda_step = (simulation::Parameter().precalclam.max_lam -
             simulation::Parameter().precalclam.min_lam) /
             (simulation::Parameter().precalclam.nr_lambdas-1);
    os << std::setw(20) << "lambda_step: " << lambda_step << "\n";
    double lambda_s = simulation::Parameter().perturbation.lambda;
    unsigned int lambda_index = (lambda_s - simulation::Parameter().precalclam.min_lam) /
                    lambda_step;
    os << "\n" << "\n";
    os << std::setw(20) << "lambda_s: " << lambda_s << "\n";
    os << std::setw(20) << "lambda_index: " << lambda_index << "\n";
    os << std::setw(20) << type + "A_e_lj (lambda_s)";
    for(unsigned int i=0; i < numenergygroups; i++) os << std::setw(12) << energroup[i];
    os << "\n";
    for(unsigned int j=0; j < numenergygroups; j++) {
      os << std::setw(20) << energroup[j];
      for(unsigned int i=0; i<j; i++) os << std::setw(12) << " ";
      for(unsigned int i=j; i < numenergygroups; i++){
        os << std::setw(12) << e.A_lj_energy[lambda_index][i][j];
      }
      os << "\n";
    }

    os << "\n" << "\n";
    os << std::setw(20) << type + "B_e_lj (lambda_s)";
    for(unsigned int i=0; i < numenergygroups; i++) os << std::setw(12) << energroup[i];
    os << "\n";
    for(unsigned int j=0; j < numenergygroups; j++) {
      os << std::setw(20) << energroup[j];
      for(unsigned int i=0; i<j; i++) os << std::setw(12) << " ";
      for(unsigned int i=j; i < numenergygroups; i++){
        os << std::setw(12) << e.B_lj_energy[lambda_index][i][j];
      }
      os << "\n";
    }

    os << "\n" << "\n";
    os << std::setw(20) << type + "A_e_crf (lambda_s)";
    for(unsigned int i=0; i < numenergygroups; i++) os << std::setw(12) << energroup[i];
    os << "\n";
    for(unsigned int j=0; j < numenergygroups; j++) {
      os << std::setw(20) << energroup[j];
      for(unsigned int i=0; i<j; i++) os << std::setw(12) << " ";
      for(unsigned int i=j; i < numenergygroups; i++){
        os << std::setw(12) << e.A_crf_energy[lambda_index][i][j];
      }
      os << "\n";
    }
    
    os << "\n" << "\n";
    os << std::setw(20) << type + "B_e_crf (lambda_s)";
    for(unsigned int i=0; i < numenergygroups; i++) os << std::setw(12) << energroup[i];
    os << "\n";
    for(unsigned int j=0; j < numenergygroups; j++) {
      os << std::setw(20) << energroup[j];
      for(unsigned int i=0; i<j; i++) os << std::setw(12) << " ";
      for(unsigned int i=j; i < numenergygroups; i++){
        os << std::setw(12) << e.B_crf_energy[lambda_index][i][j];
      }
      os << "\n";
    }//ANITA
*/    
    os << "\n" << std::setw(20) << type + "LS (real)";
    
    for(unsigned int i=0; i < numenergygroups; i++) os << std::setw(12) << energroup[i];
    os << "\n";
    for(unsigned int j=0; j < numenergygroups; j++) {
      os << std::setw(20) << energroup[j];
      for(unsigned int i=0; i<j; i++) os << std::setw(12) << " ";
      for(unsigned int i=j; i < numenergygroups; i++){
	os << std::setw(12) << e.ls_real_energy[i][j];
      }
      os << "\n";
    }    

    /*
    energy groups not implemented for k-space energy.
    os << "\n" << std::setw(20) << type + "LS (k-space)";
    
    for(unsigned int i=0; i < numenergygroups; i++) os << std::setw(12) << energroup[i];
    os << "\n";
    for(unsigned int j=0; j < numenergygroups; j++) {
      os << std::setw(20) << energroup[j];
      for(unsigned int i=0; i<j; i++) os << std::setw(12) << " ";
      for(unsigned int i=j; i < numenergygroups; i++){
	os << std::setw(12) << e.ls_k_energy[i][j];
      }
      os << "\n";
    }        
    */
    os << "\n" << std::setw(20) << type + "Self-energy";
    for(unsigned int i=0; i < numenergygroups; i++) os << std::setw(12) << e.self_energy[i];

    os << "\n" << "\n";
    os << std::setw(20) << type + "SPECIAL";

    for(unsigned int i=0; i < numenergygroups; i++) os << std::setw(12) << energroup[i];

    os << "\n" << std::setw(20) << type + "SASA";
    for(unsigned int i=0; i < numenergygroups; i++) os << std::setw(12) << e.sasa_energy[i];

    os << "\n" << std::setw(20) << type + "SASA Volume";
    for (unsigned int i = 0; i < numenergygroups; i++) os << std::setw(12) << e.sasa_volume_energy[i];

    os << "\n" << std::setw(20) << type + "Constraints";
    for(unsigned int i=0; i < numenergygroups; i++) os << std::setw(12) << e.constraints_energy[i];

    os << "\n" << std::setw(20) << type + "Distanceres";
    for(unsigned int i=0; i < numenergygroups; i++) os << std::setw(12) << e.distanceres_energy[i];

    os << "\n" << std::setw(20) << type + "Disfieldres";
    for(unsigned int i=0; i < numenergygroups; i++) os << std::setw(12) << e.disfieldres_energy[i];

    os << "\n" << std::setw(20) << type + "Angrest";
    for(unsigned int i=0; i < numenergygroups; i++) os << std::setw(12) << e.angrest_energy[i];

    os << "\n" << std::setw(20) << type + "Dihrest";
    for(unsigned int i=0; i < numenergygroups; i++) os << std::setw(12) << e.dihrest_energy[i];

    os << "\n" << std::setw(20) << type + "Posrest";
    for(unsigned int i=0; i < numenergygroups; i++) os << std::setw(12) << e.posrest_energy[i];

    os << "\n" << std::setw(20) << type + "Jrest";
    for(unsigned int i=0; i < numenergygroups; i++) os << std::setw(12) << e.jvalue_energy[i];

    os << "\n" << std::setw(20) << type + "RDCrest";
    for(unsigned int i=0; i < numenergygroups; i++) os << std::setw(12) << e.rdc_energy[i];


    os << "\nEND\n";
    
  }

  /*
   * Print a matrix.
   */
  void print_MATRIX(std::ostream &os, math::Matrix const &m,
		    std::string const title)
  {
    os.precision(5);
    os.setf(std::ios::scientific, std::ios::floatfield);

    os << title << "\n";
    for(int a=0; a<3; ++a){
      os << "\t";
      for(int b=0; b<3; ++b){
	os << std::setw(15) << m(a,b);
      }
      os << "\n";
    }
    os << "END\n";

    os.setf(std::ios::fixed, std::ios::floatfield);

  }

  void print_TIMESTEP(std::ostream &os, double const steps, double const time)
  {
    os.setf(std::ios::fixed, std::ios::floatfield);
    os.precision(5);
    
    os << "TIMESTEP"
       << std::setw(15) << steps
       << " "
       << std::setw(14) << time
       << "\n";

  }

  /**
   * @todo don't print rotation in periodic system
   */
  void print_CENTREOFMASS(std::ostream &os, 
			  double const ekin_trans, 
			  double const ekin_rot)
  {
    os.setf(std::ios_base::scientific, std::ios_base::floatfield);
    os.precision(5);
    
    os << "CENTREOFMASS\n"
       << std::setw(15) << "E-KIN trans " << std::setw(15) << ekin_trans <<"\n"
       << std::setw(15) << "E-KIN rot "   << std::setw(15) << ekin_rot <<"\n"
       << std::setw(15) << "E-KIN COM " << std::setw(15) << ekin_trans + ekin_rot
       << "\nEND\n";

  }
  
  void print_sasa(std::ostream &os, topology::Topology const & topo,
                    configuration::Configuration const & conf,
                    simulation::Simulation const & sim,
		    std::string const title)
  {
    os << "\n";
    os.precision(7);
    //std::setw(13) << std::setprecision(4) << std::scientific
    if (sim.param().sasa.switch_volume) {
      os << title << "\n";
      os << "ATOM" << "\t" << "SASA" << "\t\t" << "VOLUME" << "\n";
      for (unsigned int i = 0; i < topo.sasa_parameter().size(); ++i) {
        os << topo.sasa_parameter(i).atom + 1 << "\t"
         << std::scientific << conf.old().sasa_area[i]
        << "\t" << conf.old().sasa_buriedvol[i] << "\n";
      }
      os << "\n";
      os << "TOTAL SASA\tTOTAL BURIED VOLUME\tTOTAL VOLUME\n";
      os << conf.old().sasa_tot << "\t" << conf.old().sasa_buriedvol_tot <<
      "\t\t" << topo.sasa_volume_tot() << "\n\n";
    } else {
      os << "SOLVENT ACCESSIBLE SURFACE AREAS" << "\n";
      os << "ATOM" << "\t" << "SASA" << "\n";
      for (unsigned int i = 0; i < topo.sasa_parameter().size(); ++i) {
        os << topo.sasa_parameter(i).atom + 1  << "\t"
        << std::scientific << conf.old().sasa_area[i] << "\n";
      }
      os << "\n";
      os << "TOTAL SASA" << "\n";
      os << conf.old().sasa_tot << "\n\n";
    }
  }

  void print_sasa_avg(std::ostream &os, topology::Topology const & topo,
                      std::vector<double> const &s,
		      std::vector<double> const &v, double & stot, double & vtot,
                      std::string const title, int volume) {
    os.precision(7);

    if (volume)
      os << title << "\n";
    else os << "SASA AVERAGES" << "\n";

    if (volume) {
      os << "ATOM" << "\t" << "<SASA>" << "\t\t" << "<BURIED VOLUME>" << "\n";
      for (unsigned int i = 0; i < s.size(); ++i) {
        os << topo.sasa_parameter(i).atom + 1 << "\t"
        << std::scientific << s[i] << "\t" << v[i] << "\n";
      }
    } else {
      os << "ATOM" << "\t" << "<SASA>" << "\n";
      for (unsigned int i = 0; i < s.size(); ++i) {
        os << topo.sasa_parameter(i).atom + 1 << "\t"
        << std::scientific << s[i] << "\n";
      }
    }
    if (volume){
    os << "TOTAL AVERAGE SASA\tTOTAL AVERAGE BURIED VOLUME\n";
    os << stot << "\t\t" << vtot << "\n";
    os << "END\n";
    }
    else {
      os << "TOTAL AVERAGE SASA\n";
      os << stot << "\n";
      os << "END\n";
    }
    os << "\n";
}

  void print_sasa_fluct(std::ostream &os, topology::Topology const & topo,
                        std::vector<double> const &s,
		        std::vector<double> const &v, double & stot, double & vtot,
                        std::string const title, int volume)
{
    os.precision(7);

    if (volume)
      os << title << "\n";
    else os << "SASA FLUCTUATIONS" << "\n";

    if (volume) {
      os << "ATOM" << "\t" << "<<SASA>>" << "\t" << "<<BURIED VOLUME>>" << "\n";
      for (unsigned int i = 0; i < s.size(); ++i) {
        os << topo.sasa_parameter(i).atom + 1 << "\t"
        << std::scientific << s[i] << "\t" << v[i] << "\n";
      }
    } else {
      os << "ATOM" << "\t" << "<<SASA>>" << "\n";
      for (unsigned int i = 0; i < s.size(); ++i) {
        os << topo.sasa_parameter(i).atom + 1 << "\t"
        << std::scientific << s[i] << "\n";
      }
    }
    if (volume){
    os << "TOTAL SASA FLUCTUATIONS" << "\t" << "TOTAL BURIED VOLUME FLUCTUATIONS" << "\n";
    os << stot << "\t\t" << vtot << "\n";
    os << "END\n";
    }
    else {
      os << "TOTAL SASA FLUCTUATIONS" << "\n";
      os << stot << "\n";
      os << "END\n";
    }
    os << "\n";
}

  /*
   * Print the MINIMISATION block.
   */
  /*void print_MINIMISATION(std::ostream &os,
		    configuration::Forces const &f, // Here we will read force RMSD and MAX force
		    ???, // Here we will read total number of iterations, doubling and interpolations
		    bool last,
		    bool minimised)
  {
    os << "MINIMISATION\n";
    if (last) {
      if (minimised) {
        os << "MINIMUM REACHED\n";
      }
      else {
        os << "MINIMUM CRITERION NOT MET\n";
      }

      //os << "Total interaction calculations : " << total_iterations + 1 << "\n";
      //os << "Total search interval doublings: " << total_doubling << "\n";
      //os << "Total interpolations           : " << total_interpolations << "\n";
      os << "Final RMS force                : " << f.rms_force << "\n";
      os << "Final MAX force                : " << f.max_force << "\n";
    }
    else {
      os << "RMS force:\t" << rms_force << "\n";
      os << "MAX force:\t" << max_force << "\n";
    }
    os << "END\n";

  }*/


} // io
