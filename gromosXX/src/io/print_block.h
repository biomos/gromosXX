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

} // io


  
