/**
 * @file print_block.cc
 * routines to print out the various blocks.
 */

#include <util/stdheader.h>
#include <fstream>

#include <topology/core/core.h>
#include <topology/topology.h>
#include <simulation/multibath.h>
#include <simulation/parameter.h>
#include <simulation/simulation.h>
#include <configuration/energy.h>
#include <configuration/energy_average.h>
#include <configuration/configuration.h>
#include <interaction/interaction.h>
#include <interaction/interaction_types.h>

#include <io/argument.h>
#include <io/blockinput.h>
#include <io/instream.h>
#include <io/configuration/inframe.h>
#include <io/configuration/in_configuration.h>
#include <io/topology/in_topology.h>
#include <io/topology/in_perturbation.h>
#include <io/parameter/in_parameter.h>

#include <algorithm/algorithm.h>
#include <algorithm/algorithm/algorithm_sequence.h>
#include <algorithm/create_md_sequence.h>

#include <interaction/forcefield/forcefield.h>

#include "print_block.h"

namespace io
{

  /** 
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
	 << std::setw(12) << it->last_molecule + 1
	 << std::setw(12) << it->com_bath
	 << std::setw(12) << it->ir_bath
	 << "\n";
    }
    
    os << "END\n";
    
  }

  /**
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
  
    for(size_t i=0; it != to; ++it, ++i){
      
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
  
  /** 
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
       << std::setw(12) << "EKIN"
       << std::setw(12) << "EKIN-MOL"
       << std::setw(12) << "EKIN-IR"
       << std::setw(10) << "TEMP"
       << std::setw(10) << "TEMP-MOL"
       << std::setw(10) << "TEMP-IR"
       << "\n";
  
    double avg_temp0 = 0, avg_tau = 0, sum_dof = 0, sum_soluc = 0,
      sum_solvc = 0, sum_ekin = 0, tau_dof = 0,
      sum_com_ekin = 0, sum_ir_ekin = 0,
      sum_ir_dof = 0, sum_com_dof = 0;

    std::vector<simulation::bath_struct>::const_iterator
      it = bath.begin(),
      to = bath.end();
  
    for(size_t i=0; it != to; ++it, ++i){
      
      const double e_kin = energy.kinetic_energy[i];
      const double e_kin_com = energy.com_kinetic_energy[i];
      const double e_kin_ir = energy.ir_kinetic_energy[i];

      os << std::setw(10) << i
	 << std::setw(12) << std::setprecision(4) << std::scientific 
	 << e_kin
	 << std::setw(12) 
	 << e_kin_com
	 << std::setw(12) 
	 << e_kin_ir
	 << std::setprecision(2) << std::fixed;
      if (e_kin == 0){
	os << std::setw(10) << 0;
      }
      else{
	os << std::setw(10) 
	   << 2 * e_kin / (math::k_Boltzmann * it->dof);
      }
      if (e_kin_com == 0){
	os << std::setw(10) << 0;
      }
      else{
	os << std::setw(10) 
	   << 2 * e_kin_com / 
	  (math::k_Boltzmann * it->com_dof);
      }
      if (e_kin_ir == 0){
	os << std::setw(10) << 0;
      }
      else{
	os << std::setw(10) 
	   << 2 * e_kin_ir / 
	  (math::k_Boltzmann * it->ir_dof);
      }

      if (it->tau != -1){
	tau_dof += it->dof;
	avg_temp0 += it->temperature * it->dof;
	avg_tau += it->tau * it->dof;
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
       << "-----------------------------------------\n";

    os << std::setw(10) << "Avg"
       << std::setw(12) << std::setprecision(4) << std::scientific
       << sum_ekin
       << std::setw(12) << sum_com_ekin
       << std::setw(12) << sum_ir_ekin
       << std::setprecision(2) << std::fixed
       << std::setw(10) << 2 * sum_ekin / (math::k_Boltzmann * sum_dof)
       << std::setw(10) << 2 * sum_com_ekin / (math::k_Boltzmann * sum_com_dof)
       << std::setw(10) << 2 * sum_ir_ekin / (math::k_Boltzmann * sum_ir_dof)
       << "\n";
    
    os << "END\n";

  }

  /**
   * Print the PCOUPLE block.
   */
  void print_PCOUPLE(std::ostream &os, bool calc, 
		     math::pressure_scale_enum scale,
		     math::Matrix pres0, double comp, 
		     double tau, math::virial_enum vir)
  {
    os << "PCOUPLE\n";
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
      case math::pcouple_full_anisotropic:
	os << std::setw(12) << "full";
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
    
    os << "\n" << std::setw(23) << "REFERENCE PRESSURE" << "\n";
    
    for(int i=0; i<3; ++i){
      for(int j=0; j<3; ++j){
	os << std::setw(12) << pres0(i,j);
      }
      os << "\n";
    }
    
    os << "END\n";

  }

  /**
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
    for(int i=0; i<3; ++i){
      for(int j=0; j<3; ++j)
	os << std::setw(15) << conf.old().virial_tensor(i,j);
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
       << 0.0 // sys.periodicity().volume()
       << "\n";
    
    os << "\nEND\n";

    os.precision(5);
    os.setf(std::ios_base::fixed, std::ios_base::floatfield);
    
  }

  /**
   * Print the ENERGY block.
   */
  void print_ENERGY(std::ostream &os,
		    configuration::Energy const &e,
		    std::vector<size_t> const &energy_groups,
		    std::string const title)
  {

    int numenergygroups=e.bond_energy.size();
 
    std::vector<std::string> energroup;
   
    int b=1;
    
    for(int i=0; i<numenergygroups; i++){

      std::ostringstream ostring;
      ostring << b << "-" << energy_groups[i]+1;
      energroup.push_back(ostring.str());
      b = energy_groups[i]+2;
    }
        
    os << title << "\n";

    os.precision(4);
    os.setf(std::ios_base::scientific, std::ios_base::floatfield);
    os << "Total      : " << std::setw(12) << e.total << "\n";
    os << "Kinetic    : " << std::setw(21) << e.kinetic_total << "\n";
    os << "Potential  : " << std::setw(21) << e.potential_total << "\n";
    os << "Covalent   : " << std::setw(30) << e.bonded_total << "\n";
    os << "Bonds      : " << std::setw(39) << e.bond_total << "\n";
    os << "Angles     : " << std::setw(39) << e.angle_total << "\n";
    os << "Improper   : " << std::setw(39) << e.improper_total << "\n";
    os << "Dihedral   : " << std::setw(39) << e.dihedral_total << "\n";
    os << "Non-bonded : " << std::setw(30) << e.nonbonded_total  << "\n";
    os << "Vdw        : " << std::setw(39) << e.lj_total << "\n";
    os << "El (RF)    : " << std::setw(39) << e.crf_total  << "\n";
    os << "Special    : " << std::setw(21) << e.special_total << "\n";
    os << "\n";

    os << std::setw(10) << "COV";
    
    for(int i=0; i < numenergygroups; i++) os << std::setw(12) << energroup[i];
    os << "\n" << std::setw(12) << "bonds";
    for(int i=0; i < numenergygroups; i++) os << std::setw(12) << e.bond_energy[i];
    os << "\n" << std::setw(12) << "angles";
    for(int i=0; i < numenergygroups; i++) os << std::setw(12) << e.angle_energy[i];
    os << "\n" << std::setw(12) << "dihedrals";
    for(int i=0; i < numenergygroups; i++) os << std::setw(12) << e.dihedral_energy[i];
    os << "\n" << std::setw(12) << "impropers";
    for(int i=0; i < numenergygroups; i++) os << std::setw(12) << e.improper_energy[i];
    os << "\n" << "\n";
    os << std::setw(10) << "VDW";
    
    for(int i=0; i < numenergygroups; i++) os << std::setw(12) << energroup[i];
    os << "\n";
    for(int j=0; j < numenergygroups; j++) {
      os << std::setw(12) << energroup[j];
      for(int i=0; i<j; i++) os << std::setw(12) << " ";
      for(int i=j; i < numenergygroups; i++){
	if(i==j)
	  os << std::setw(12) << e.lj_energy[i][j];
	else 
	  os << std::setw(12) << e.lj_energy[i][j] + e.lj_energy[j][i] ;
      }
      os << "\n";
    }
    os << "\n" << std::setw(10) << "CRF";
    
    for(int i=0; i < numenergygroups; i++) os << std::setw(12) << energroup[i];
    os << "\n";
    for(int j=0; j < numenergygroups; j++) {
      os << std::setw(12) << energroup[j];
      for(int i=0; i<j; i++) os << std::setw(12) << " ";
      for(int i=j; i < numenergygroups; i++){
	if(i==j)
	  os << std::setw(12) << e.crf_energy[i][j];
	else
	  os << std::setw(12) << e.crf_energy[i][j] +  e.crf_energy[j][i];
      }
      os << "\n";
    }

    os << "\n" << "\n";
    os << std::setw(10) << "SPECIAL";

    os << "\n" << std::setw(12) << "Posrest";
    // for(int i=0; i < numenergygroups; i++) os << std::setw(12) << e.energy_posrest[i];


    os << "\nEND\n";
    
  }

  /**
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
    os.precision(9);
    
    os << "TIMESTEP\n"
       << std::setw(15) << steps
       << std::setw(15) << time
       << "\nEND\n";

  }

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
  
} // io
