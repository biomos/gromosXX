/**
 * @file print_block.tcc
 * routines to print out the various blocks.
 */

namespace io
{

  /** 
   * Print the DOF COUPLING table of MULTIBATH block.
   */
  inline std::ostream & 
  print_MULTIBATH_COUPLING(std::ostream &os,
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
    return os;
    
  }

  /**
   * Print DEGREESOFFREEDOM block.
   */
  inline std::ostream & print_DEGREESOFFREEDOM(std::ostream &os,
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

    os << "    --------------------------------------------------------------------------\n";
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
    return os;
    
  }
  
  /** 
   * Print the MULTIBATH block.
   */
  inline std::ostream & print_MULTIBATH(std::ostream &os,
					simulation::Multibath const &bath,
					simulation::Energy const &energy)
  {
    os << "MULTIBATH\n";
  
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
      
      os << std::setw(10) << i
	 << std::setw(12) << std::setprecision(4) << std::scientific 
	 << energy.kinetic_energy[i]
	 << std::setw(12) << energy.com_kinetic_energy[i]
	 << std::setw(12) << energy.ir_kinetic_energy[i]
	 << std::setprecision(2) << std::fixed;
      if (energy.kinetic_energy[i] == 0){
	os << std::setw(10) << 0;
      }
      else{
	os << std::setw(10) 
	   << 2 * energy.kinetic_energy[i] / (math::k_Boltzmann * it->dof);
      }
      if (energy.com_kinetic_energy[i] == 0){
	os << std::setw(10) << 0;
      }
      else{
	os << std::setw(10) 
	   << 2 * energy.com_kinetic_energy[i] / 
	  (math::k_Boltzmann * it->com_dof);
      }
      if (energy.ir_kinetic_energy[i] == 0){
	os << std::setw(10) << 0;
      }
      else{
	os << std::setw(10) 
	   << 2 * energy.ir_kinetic_energy[i] / 
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
      sum_ekin += energy.kinetic_energy[i];
      sum_com_ekin += energy.com_kinetic_energy[i];
      sum_ir_ekin += energy.ir_kinetic_energy[i];
      
      os << "\n";

    }

    os << "    --------------------------------------------------------------------------\n";
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
    return os;
  }

  /**
   * Print the PRESSURE block.
   */
  template<math::boundary_enum b>
  inline std::ostream &print_PRESSURE(std::ostream &os,
				      simulation::System<b> const & sys)
  {
    os << "PRESSURE\n";
    os.precision(5);
    os.setf(std::ios_base::fixed, std::ios_base::floatfield);
    
    os << "\tmolecular kinetic energy:\n\t";
    for(int i=0; i<3; ++i){
      for(int j=0; j<3; ++j)
	os << std::setw(12) << sys.molecular_kinetic_energy()(i,j);
      os << "\n\t";
    }
    
    os << "\n\tvirial\n\t";
    for(int i=0; i<3; ++i){
      for(int j=0; j<3; ++j)
	os << std::setw(12) << sys.virial()(i,j);
      os << "\n\t";
    }

    os << "\n\tpressure tensor\n\t";
    for(int i=0; i<3; ++i){
      for(int j=0; j<3; ++j)
	os << std::setw(12) << sys.pressure()(i,j);
      os << "\n\t";
    }
    os << "\n\tpressure: " 
       << (sys.pressure()(0,0)+sys.pressure()(1,1)+sys.pressure()(2,2))/3 
       << "\n";

    os << "\nEND\n";
    
    return os;
  }

  /**
   * Print the ENERGY block.
   */
  inline std::ostream &print_ENERGY(std::ostream &os,
				    simulation::Energy const &e,
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
    // os << "Temperature: " << std::setw(21) << v_kin/(0.5*k_Boltzmann*Ndf) << "\n";
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
    os << "Constraint : " << std::setw(30) << e.constraint_energy << "\n";
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
    
  
    return os;
  }

  inline std::ostream & print_TIMESTEP(std::ostream &os, double const steps, double const time)
  {
    os.setf(std::ios::fixed, std::ios::floatfield);
    os.precision(9);
    
    os << "TIMESTEP\n"
       << std::setw(15) << steps
       << std::setw(15) << time
       << "\nEND\n";

    return os;
  }

  
} // io
