/**
 * @file InInput.tcc
 * implements methods of InInput.
 */

#undef MODULE
#define MODULE io
#undef SUBMODULE
#define SUBMODULE input

#include "../../debug.h"

/**
 * read the stream into blocks.
 */
inline void io::InInput::read_stream()
{
  std::vector<std::string> buffer;
  
  while(!stream().eof()){

    try{
      io::getblock(stream(), buffer);
    }
    catch(std::runtime_error e){
      break;
    }
    
    m_block[buffer[0]] = buffer;    
    buffer.clear();
    
  }
}

/**
 * Store standard parameters in the simulation.
 */
template<typename t_topology, typename t_system>
inline io::InInput & io::InInput
::operator>>(simulation::Simulation<t_topology, t_system> &sim)
{
  // let's for now only do the 'numbers only' stuff
  std::vector<std::string> buffer;
  std::vector<std::string>::const_iterator it;
  
  { // PLIST
    buffer = m_block["PLIST"];
    
    it = buffer.begin() + 1;
    _lineStream.clear();
    _lineStream.str(*it);
    
    int i, update_step;
    double rcutp, rcutl;
    
    _lineStream >> i >> update_step
		>> rcutp >> rcutl;

    if (_lineStream.fail())
      throw std::runtime_error("bad line in PLIST block");

    DEBUG(7, "pairlist update=" << update_step);
    DEBUG(7, "setting short cutoff=" << rcutp << " long cutoff=" << rcutl);
    
    sim.nonbonded().update(update_step);
    sim.nonbonded().cutoff_short(rcutp);
    sim.nonbonded().cutoff_long(rcutl);
  }
  
  { // LONGRANGE
    buffer = m_block["LONGRANGE"];

    it = buffer.begin() + 1;
    _lineStream.clear();
    _lineStream.str(*it);
    
    double epsilon, kappa, cutoff;
    
    _lineStream >> epsilon >> kappa >> cutoff;

    if (_lineStream.fail())
      io::messages.add("bad line in LONGRANGE block",
		       "InInput", io::message::error);

    if (!_lineStream.eof())
      io::messages.add("End of line not reached, but should have been: \n" + *it +  "\n",
		       "InInput", io::message::warning);


    sim.nonbonded().RF_constant(epsilon, kappa, cutoff);

    DEBUG(7, "calculating Crf: epsilon= " << epsilon << " kappa= "
	  << kappa << " cutoff= " << cutoff << endl << "RF_constant= " 
	  << sim.nonbonded().RF_constant());
    
  } // LONGRANGE
  
  { // SUBMOLECULE
    buffer = m_block["SUBMOLECULES"];
    
    // std::cerr << "reading SUBMOLECULES" << std::endl;

    if (buffer.begin() == buffer.end()){
      io::messages.add("empty SUBMOLECULES block",
		       "InInput", io::message::error);
    }
    else{
      
      std::string submol;
      concatenate(buffer.begin()+1, buffer.end()-1, submol);
      
      _lineStream.clear();
      _lineStream.str(submol);
    
      size_t num;
      _lineStream >> num;
      // std::cerr << "num: " << num << std::endl;
      
      size_t m;

      for(size_t i=0; i<num; ++i){
	_lineStream >> m;
	sim.topology().molecules().push_back(m);
	// std::cerr << "\t" << m << std::endl;
      }
    
      if (_lineStream.fail())
	io::messages.add("bad line in SUBMOLECULES block",
			 "InInput", io::message::error);
    }
    
  } // SUBMOLECULE

  { // TEMPERATURE COUPLING
    
    // is there a MULTIBATH block
    buffer = m_block["MULTIBATH"];
    
    it = buffer.begin();
    if (it != buffer.end()){

      io::messages.add("using MULTIBATH block",
		       "InInput", io::message::notice);
      ++it;
      _lineStream.clear();
      _lineStream.str(*it);
      
      int num;
      int last;
      double temp, tau;

      _lineStream >> num;
      ++it;
      
      for(int i=0; i<num; ++i, ++it){
	_lineStream.clear();
	_lineStream.str(*it);
	
	_lineStream >> last >> temp >> tau;

	sim.multibath().add_bath(last, temp, tau);
	
      }
      
      if (_lineStream.fail()){
	io::messages.add("bad line in MULTIBATH block",
			 "InInput", io::message::error);
      }

    }
    else{
      // try a TCOUPLE block
      
      buffer = m_block["TCOUPLE"];
      if (buffer.size()){

	int ntt[3];
	double temp[3];
	double tau[3];
	
	it = buffer.begin()+1;

	for(int i=0; i<3; ++i, ++it){
	  _lineStream.clear();
	  _lineStream.str(*it);
	  _lineStream >> ntt[i] >> temp[i] >> tau[i];
	}

	if (_lineStream.fail()){
	  io::messages.add("bad line in TCOUPLE block",
			   "InInput", io::message::error);
	  return *this;
	}	

	// need the number of solvents...
	// has to be adapted in case of multiple solvents...
	int nsm;
	read_SYSTEM(nsm);

	// the last solute atom (and really the last index)
	int last_solute = sim.topology().num_solute_atoms() - 1;
	// the last solvent atom
	int last_solvent = sim.topology().solvent(0).num_atoms() * nsm + last_solute;

	// cases to handle
	// 0 0 0
	if (ntt[0] == 0 && ntt[1] == 0 && ntt[2] == 0){
	  // nothing
	}
	// 2 -2 0
	else if (ntt[0] == 2 && ntt[1] == -2 && ntt[2] == 0){
	  sim.multibath().add_bath(last_solute, temp[0], tau[0]);
	}
	// 2 -2 1
	else if (ntt[0] == 2 && ntt[1] == -2 && ntt[2] == 1){
	  sim.multibath().add_bath(last_solute, temp[0], tau[0]);
	  sim.multibath().add_bath(last_solvent, temp[2], tau[2]);
	}
	// 3 3 3
	else if (ntt[0] == 3 && ntt[1] == 3 && ntt[2] == 3){
	  sim.multibath().add_bath(last_solvent, temp[0], tau[0]);
	}
	// rest is not handled!
	else{
	  io::messages.add("TCOUPLE ntt combination not handled",
			   "InInput", io::message::error);
	  return *this;
	}

      }
      else{
	// no TCOUPLE block
      }
      
    }
    
  }

  return *this;
}

/**
 * read the SYSTEM block.
 */
inline void io::InInput::read_SYSTEM(int &nsm)
{
  std::vector<std::string> buffer;
  buffer = m_block["SYSTEM"];
  _lineStream.clear();
  _lineStream.str(buffer[1]);
  
  int npm;
  _lineStream >> npm >> nsm;
 
  if (_lineStream.fail())
    io::messages.add("bad line in SYSTEM block",
		       "InInput", io::message::error);

  if (!_lineStream.eof())
    io::messages.add("End of line not reached in SYSTEM block, but should have been: \n" + _lineStream.str() +  "\n",
		       "InInput", io::message::warning);

  
  if (npm != 1)
    io::messages.add("SYSTEM: only NPM=1 allowed",
		     "io::InInput::read_SYSTEM",
		     io::message::error);
  
} 

/**
 * read the STEP block.
 */
inline void io::InInput::read_STEP(int &num_steps, double &t0, double &dt)
{
  std::vector<std::string> buffer;
  std::vector<std::string>::const_iterator it;
  
  buffer = m_block["STEP"];
  
  it = buffer.begin() + 1;
  _lineStream.clear();
  _lineStream.str(*it);
  
  _lineStream >> num_steps >> t0 >> dt;
  
  if (_lineStream.fail())
    io::messages.add("bad line in STEP block",
		       "InInput", io::message::error);

  if (!_lineStream.eof())
    io::messages.add("End of line not reached in STEP block, but should have been: \n" + *it +  "\n",
		       "InInput", io::message::warning);
}

/**
 * the SHAKE block.
 */
inline void io::InInput::read_SHAKE(int &ntc, double &tolerance)
{
  std::vector<std::string> buffer;
  std::vector<std::string>::const_iterator it;
  
  buffer = m_block["SHAKE"];
  
  it = buffer.begin() + 1;
  _lineStream.clear();
  _lineStream.str(*it);
  
  _lineStream >> ntc >> tolerance;
  
  if (_lineStream.fail())
    io::messages.add("bad line in SHAKE block",
		       "InInput", io::message::error);

  if (!_lineStream.eof())
    io::messages.add("End of line not reached in SHAKE block, but should have been: \n" + *it +  "\n",
		       "InInput", io::message::warning);
}

/**
 * read the PRINT and WRITE block.
 */
inline void io::InInput::read_PRINT(int &print_trajectory, 
				    int &print_velocity,
				    int &print_energy)
{
  std::vector<std::string> buffer;
  std::vector<std::string>::const_iterator it;
  
  buffer = m_block["PRINT"];
  it = buffer.begin() + 1;
  _lineStream.clear();
  _lineStream.str(*it);
  
  int com, dih_monitoring;
  _lineStream >> print_energy >> com >> dih_monitoring;

  if (_lineStream.fail())
    io::messages.add("bad line in PRINT block",
		       "InInput", io::message::error);

  if (!_lineStream.eof())
    io::messages.add("End of line not reached in PRINT, but should have been: \n" + *it +  "\n",
		       "InInput", io::message::warning);

  buffer = m_block["WRITE"];
  it = buffer.begin() + 1;
  _lineStream.clear();
  _lineStream.str(*it);
  
  int format, selection, energy, free_energy;
  _lineStream >> print_trajectory >> selection >> print_velocity
	      >> energy >> free_energy >> format;
  
  if (_lineStream.fail())
    io::messages.add("bad line in WRITE block",
		       "InInput", io::message::error);

  if (!_lineStream.eof())
    io::messages.add("End of line not reached in WRITE block, but should have been: \n" + *it +  "\n",
		       "InInput", io::message::warning);
    
}

/**
 * read the PCOUPLE block.
 */
inline void io::InInput::read_PCOUPLE(int &ntp, double &pres0,
				      double &comp, double &tau)
{
  std::vector<std::string> buffer;
  
  buffer = m_block["PCOUPLE"];
  _lineStream.clear();
  _lineStream.str(buffer[1]);
  
  _lineStream >> ntp >> pres0 >> comp >> tau;
  
  if (_lineStream.fail())
    io::messages.add("bad line in PCOUPLE block",
		     "InInput", io::message::error);
}

/**
 * read the BOUNDARY block.
 */
inline void io::InInput::read_BOUNDARY(int &ntb, int &nrdbox)
{
  std::vector<std::string> buffer;
  buffer = m_block["BOUNDARY"];
  _lineStream.clear();
  _lineStream.str(buffer[1]);

  double b1, b2, b3, beta;

  _lineStream >> ntb >> b1 >> b2 >> b3 >> beta >> nrdbox;

  if (_lineStream.fail())
    io::messages.add("bad line in BOUNDARY block",
		     "InInput", io::message::error);
}

/**
 * read FORCE block.
 */
inline void io::InInput::read_FORCE(int &do_bond, int &do_angle,
				    int &do_improper, int &do_dihedral,
				    int &do_nonbonded)
{
  std::vector<std::string> buffer;
  std::vector<std::string>::const_iterator it;
  
  buffer = m_block["FORCE"];
  
  it = buffer.begin()+1;
  _lineStream.clear();
  _lineStream.str(*it);
  
  int bondH, angleH, impH, dihedralH, charge;
  _lineStream >> bondH >> do_bond >> angleH >> do_angle
	      >> impH >> do_improper >> dihedralH >> do_dihedral
	      >> charge >> do_nonbonded;
  
  if (bondH ^ do_bond)
    io::messages.add("Force switch for bond and bond H has to be equal",
		     "InInput", io::message::error);

  if (angleH ^ do_angle)
    io::messages.add("Force switch for angle and angle H has to be equal",
		     "InInput", io::message::error);

  if (impH ^ do_improper)
    io::messages.add("Force switch for improper and improper H has to be equal",
		     "InInput", io::message::error);

  if (dihedralH ^ do_dihedral)
    io::messages.add("Force switch for dihedral and dihedral H has to be equal",
		     "InInput", io::message::error);

  if (charge ^ do_nonbonded)
    io::messages.add("Force switch for lj and charge has to be equal",
		     "InInput", io::message::error);

  if (_lineStream.fail())
    io::messages.add("bad line in FORCE block",
		       "InInput", io::message::error);

  if (!_lineStream.eof())
    io::messages.add("End of line not reached in FORCE block, but should have been: \n" + *it +  "\n",
		       "InInput", io::message::warning);
}

