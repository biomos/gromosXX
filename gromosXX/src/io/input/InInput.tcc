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
 * Store standard parameters in the simulation.
 */
template<typename t_topology, typename t_system>
inline io::InInput & io::InInput
::operator>>(simulation::Simulation<t_topology, t_system> &sim)
{
  // let's for now only do the 'numbers only' stuff
  std::vector<std::string> buffer;
  std::vector<std::string>::const_iterator it;

  DEBUG(7, "reading input (to simulation)");
  
  { // PLIST
    DEBUG(10, "pairlist block");
    buffer = m_block["PLIST"];

    if (!buffer.size())
      io::messages.add("no PLIST block in input","InInput",io::message::error);
    else{
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
  }
  
  { // LONGRANGE
    DEBUG(10, "longrange block");
    buffer = m_block["LONGRANGE"];

    if (!buffer.size()){
      io::messages.add("no LONGRANGE block in input","InInput",io::message::error);
    }
    else{

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
      
      if(cutoff < 0) {
	sim.nonbonded().RF_exclusion(true);
	io::messages.add("Reaction field contribution of excluded atoms is " 
			 "taken into account", "InInput", 
			 io::message::notice);
	cutoff *= -1;
      }
      else{
	sim.nonbonded().RF_exclusion(false);
	io::messages.add("Reaction field contribution of excluded atoms is " 
			 "NOT taken into account", "InInput", 
			 io::message::notice);
      }
      
      sim.nonbonded().RF_constant(epsilon, kappa, cutoff);
      
      DEBUG(7, "calculating Crf: epsilon= " << epsilon << " kappa= "
	    << kappa << " cutoff= " << cutoff << endl << "RF_constant= " 
	    << sim.nonbonded().RF_constant());
    }
    
  } // LONGRANGE
  
  { // SUBMOLECULES
    DEBUG(10, "submolecules block");
    buffer = m_block["SUBMOLECULES"];
    
    if (!buffer.size()){
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
	DEBUG(10, "add submol " << m);
      }
    
      if (_lineStream.fail())
	io::messages.add("bad line in SUBMOLECULES block",
			 "InInput", io::message::error);
    }
    
  } // SUBMOLECULES

  { // TEMPERATURE COUPLING
    DEBUG(10, "TEMPERATURE COUPLING block");

    // is there a MULTIBATH block
    buffer = m_block["MULTIBATH"];
    
    it = buffer.begin();
    if (it != buffer.end()){
      DEBUG(11, "MULTIBATH present");
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

	sim.multibath().add_bath(last-1, temp, tau);
	
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
	DEBUG(11, "TCOUPLE present");
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
	else if (ntt[0] == 3 && ntt[1] == -3 && ntt[2] == -3){
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
	DEBUG(11, "no TEMPERATURE COUPLING block");
	// no TCOUPLE block
	// that's fine, same as 0,0,0
      }
      
    }
    
  } // TEMPERATURE coupling
  
  { // ENERGY groups

    buffer = m_block["FORCE"];
    DEBUG(10, "FORCE block (energy groups)");
    if (buffer.size() < 4){
      DEBUG(7, "FORCE block wrong");
      io::messages.add("wrong FORCE block (energy groups)",
		       "InInput", io::message::error);
    }
    else{
      DEBUG(10, "ENERGY groups:");
      std::string egroup;
      concatenate(buffer.begin()+2, buffer.end()-1, egroup);
      _lineStream.clear();
      _lineStream.str(egroup);
      
      size_t num, e, atom = 0;
      _lineStream >> num;

      DEBUG(10, "num: " << num);
      
      for(size_t i=0; i<num; ++i){
	_lineStream >> e;
	sim.topology().energy_groups().push_back(e-1);
	for( ; atom < e; ++atom){
	  sim.topology().atom_energy_group().push_back(i);
	  DEBUG(11, "atom " << atom << ": " << i);
	}
      }

      if (_lineStream.fail())
	io::messages.add("bad line in ENERGYGROUP (FORCE) block",
			 "InInput", io::message::error);

      // and resize
      sim.system().energies().resize(num);
      
    }

  } // ENERGY groups

  DEBUG(7, "input read...");

  return *this;
}

/**
 * read the SYSTEM block.
 */
inline void io::InInput::read_SYSTEM(int &nsm)
{
  std::vector<std::string> buffer;
  buffer = m_block["SYSTEM"];
  if (!buffer.size()){
    io::messages.add("no SYSTEM block in input", "InInput", io::message::error);
    nsm = 0;
    return;
  }
  
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

  
  // we might need to also allow for 0...
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
  
  if (!buffer.size()){
    io::messages.add("no STEP block in input", "InInput", io::message::error);
    num_steps = 0;
    t0 = 0;
    dt = 0;
    return;
  }

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

  if (!buffer.size()){
    io::messages.add("no SHAKE block", "InInput",io::messages.notice);
    ntc = 1;
    tolerance = 0;
    return;
  }
  
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

  if (!buffer.size()){
    io::messages.add("no PRINT block", "InInput", io::message::notice);
    print_trajectory = 1;
    print_velocity = 1;
    print_energy = 1;
    return;
  }

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

  if (!buffer.size()){
    ntp = 0;
    pres0 = 0;
    comp = 0;
    tau = -1;
    return;
  }

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

  if (!buffer.size()){
    io::messages.add("no BOUNDARY block", "InInput", io::message::error);
    ntb = 0;
    nrdbox = 1;
    return;
  }

  _lineStream.clear();
  _lineStream.str(buffer[1]);

  double b1, b2, b3, beta;

  _lineStream >> ntb >> b1 >> b2 >> b3 >> beta >> nrdbox;

  if (_lineStream.fail())
    io::messages.add("bad line in BOUNDARY block",
		     "InInput", io::message::error);
}

/**
 * read the PERTURB block.
 */
inline void io::InInput::read_PERTURB(int &ntg, double &rlam, double &dlamt,
				      double &alphlj, double &alphc,int &nlam)
{
  std::vector<std::string> buffer;
  std::vector<std::string>::const_iterator it;
  
  buffer = m_block["PERTURB"];
  if (!buffer.size()){
    ntg = 0;
    rlam = 0;
    dlamt = 0;
    alphlj= 0;
    alphc = 0;
    nlam = 1;
    return;
  }
  
  it = buffer.begin()+1;
  _lineStream.clear();
  _lineStream.str(*it);
  
  int nrdgl;
  double dmu, dmut;
  double mmu;
  
  _lineStream >> ntg >> nrdgl >> rlam >> dlamt >> dmu >> dmut;
  _lineStream.clear();
  _lineStream.str(*(++it));
  _lineStream >> alphlj >> alphc >> nlam >> mmu;
  
  if (_lineStream.fail())
    io::messages.add("bad line in PERTURB block",
		     "InInput", io::message::error);

  if (nrdgl)
    io::messages.add("PERTURB: nrdgl != 0 not allowed",
		     "InInput", io::message::error);
  
  if (ntg != 0 && ntg != 1)
    io::messages.add("PERTURB: only ntg = 0 or ntg = 1 allowed",
		     "InInput", io::message::error);
  
  if (alphlj || alphc){
    std::cerr << "alphlj: " << alphlj << "alphc: " << alphc << std::endl;
    io::messages.add("PERTURB: softness not implemented",
		     "InInput", io::message::error);
  }

  if (nlam<=0){
    std::cerr << "nlam: " << nlam << std::endl;
    io::messages.add("PERTURB: nlam > 0",
		     "InInput", io::message::error);
  }
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

  if (!buffer.size()){
    io::messages.add("no FORCE block", "InInput", io::message::error);
    do_bond = 0;
    do_angle = 0;
    do_improper = 0;
    do_dihedral = 0;
    do_nonbonded = 0;
    return;
  }
  
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

