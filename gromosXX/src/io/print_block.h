/**
 * @file print_block.h
 * routines to print out the various blocks.
 */

/**
 * Print the multibath block.
 */
std::ostream & operator<<(std::ostream &os, simulation::Multibath &bath);


namespace io
{
  
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
  template<typename t_simulation>
  inline std::ostream &print_ENERGY(std::ostream &os,
				    t_simulation const &sim)
  {

    simulation::Energy const & e = sim.system().energies();

    int numenergygroups=e.bond_energy.size();
    double tot_vdw=0.0, tot_es=0.0, tot_bond=0.0, tot_angle=0.0,
      tot_dih=0.0, tot_imp=0.0, tot_nb=0.0, tot_b=0.0, tot_pot=0.0, 
      tot_posrest=0.0, tot_special=0.0, tot=0.0;
    double tot_kin=0.0;

    for(unsigned int i=0; i<sim.multibath().size(); i++)
      tot_kin += sim.multibath()[i].kinetic_energy;
 
    std::vector<std::string> energroup;
   
    int b=1;
    
    for(int i=0; i<numenergygroups; i++){
      for(int j=0; j<numenergygroups; j++){
	tot_vdw+=e.lj_energy[i][j];
	tot_es +=e.crf_energy[i][j];
      }
      tot_bond +=e.bond_energy[i];
      tot_angle+=e.angle_energy[i];
      tot_imp  +=e.improper_energy[i];
      tot_dih  +=e.dihedral_energy[i];
      // tot_posrest += e.posrest_energy[i];
      
      std::ostringstream ostring;
      ostring << b << "-" << sim.topology().energy_groups()[i]+1;
      energroup.push_back(ostring.str());
      b=sim.topology().energy_groups()[i]+2;
    }
    tot_nb = tot_vdw + tot_es;
    tot_b  = tot_bond + tot_angle + tot_dih + tot_imp;
    tot_pot= tot_nb + tot_b;
    tot_special = tot_posrest;

    tot    = tot_pot + tot_kin + tot_special;
        
    os << "ENERGIES\n";

    os.precision(4);
    os.setf(std::ios_base::scientific, std::ios_base::floatfield);
    os << "Total      : " << setw(12) << tot << endl;
    os << "Kinetic    : " << setw(21) << tot_kin << endl;
    // os << "Temperature: " << setw(21) << v_kin/(0.5*k_Boltzmann*Ndf) << endl;
    os << "Potential  : " << setw(21) << tot_pot << endl;
    os << "Covalent   : " << setw(30) << tot_b << endl;
    os << "Bonds      : " << setw(39) << tot_bond << endl;
    os << "Angles     : " << setw(39) << tot_angle << endl;
    os << "Improper   : " << setw(39) << tot_imp << endl;
    os << "Dihedral   : " << setw(39) << tot_dih << endl;
    os << "Non-bonded : " << setw(30) << tot_nb  << endl;
    os << "Vdw        : " << setw(39) << tot_vdw << endl;
    os << "El (RF)    : " << setw(39) << tot_es  << endl;
    os << "Special    : " << setw(21) << tot_special << endl;
    os << endl;

    os << setw(10) << "COV";
    
    for(int i=0; i < numenergygroups; i++) os << setw(12) << energroup[i];
    os << endl << setw(12) << "bonds";
    for(int i=0; i < numenergygroups; i++) os << setw(12) << e.bond_energy[i];
    os << endl << setw(12) << "angles";
    for(int i=0; i < numenergygroups; i++) os << setw(12) << e.angle_energy[i];
    os << endl << setw(12) << "dihedrals";
    for(int i=0; i < numenergygroups; i++) os << setw(12) << e.dihedral_energy[i];
    os << endl << setw(12) << "impropers";
    for(int i=0; i < numenergygroups; i++) os << setw(12) << e.improper_energy[i];
    os << endl << endl;
    os << setw(10) << "VDW";
    
    for(int i=0; i < numenergygroups; i++) os << setw(12) << energroup[i];
    os << endl;
    for(int j=0; j < numenergygroups; j++) {
      os << setw(12) << energroup[j];
      for(int i=0; i<j; i++) os << setw(12) << " ";
      for(int i=j; i < numenergygroups; i++){
	if(i==j)
	  os << setw(12) << e.lj_energy[i][j];
	else 
	  os << setw(12) << e.lj_energy[i][j] + e.lj_energy[j][i] ;
      }
      os << endl;
    }
    os << endl << setw(10) << "CRF";
    
    for(int i=0; i < numenergygroups; i++) os << setw(12) << energroup[i];
    os << endl;
    for(int j=0; j < numenergygroups; j++) {
      os << setw(12) << energroup[j];
      for(int i=0; i<j; i++) os << setw(12) << " ";
      for(int i=j; i < numenergygroups; i++){
	if(i==j)
	  os << setw(12) << e.crf_energy[i][j];
	else
	  os << setw(12) << e.crf_energy[i][j] +  e.crf_energy[j][i];
      }
      os << endl;
    }

    os << endl << endl;
    os << setw(10) << "SPECIAL";

    os << endl << setw(12) << "Posrest";
    // for(int i=0; i < numenergygroups; i++) os << setw(12) << e.energy_posrest[i];


    os << "\nEND\n";
    
  
    return os;
  }
  
} // io


  
