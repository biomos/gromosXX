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

    bool grid;
    int nsnb;
    double rcutp, rcutl, size;

    read_PLIST(grid, nsnb, rcutp, rcutl, size);
	
    sim.nonbonded().update(nsnb);
    sim.nonbonded().cutoff_short(rcutp);
    sim.nonbonded().cutoff_long(rcutl);
    sim.nonbonded().grid_cell_size(size);

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
      size_t last;
      size_t com_bath, ir_bath;
      double temp, tau;

      // the baths
      _lineStream >> num;
      ++it;
      
      for(int i=0; i<num; ++i, ++it){
	_lineStream.clear();
	_lineStream.str(*it);
	
	_lineStream >> temp >> tau;

	sim.multibath().add_bath(temp, tau);
	
      }
      
      if (_lineStream.fail()){
	io::messages.add("bad line in MULTIBATH block",
			 "InInput", io::message::error);
      }

      // now the ranges
      _lineStream.clear();
      _lineStream.str(*it);
      _lineStream >> num;
      ++it;
      
      for(int i=0; i<num; ++i, ++it){
	_lineStream.clear();
	_lineStream.str(*it);
	_lineStream >> last >> com_bath >> ir_bath;
	sim.multibath().add_bath_index(last - 1, com_bath, ir_bath);
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
	// 1 0 0
	else if (ntt[0] == 1 && ntt[1] == 0 && ntt[2] == 0){
	  // the baths
	  sim.multibath().add_bath(temp[0], tau[0]);
	  sim.multibath().add_bath(0, -1);
	  // the atoms in the baths
	  sim.multibath().add_bath_index(last_solute, 0, 1);
	  // and an uncoupled one for the solvent...
	  if (last_solvent != last_solute)
	    sim.multibath().add_bath_index(last_solvent, 1, 1);

	}
	// 0 1 0
	else if (ntt[0] == 0 && ntt[1] == 1 && ntt[2] == 0){
	  // the baths
	  sim.multibath().add_bath(temp[1], tau[1]);
	  sim.multibath().add_bath(0, -1);
	  // the atoms in the baths
	  sim.multibath().add_bath_index(last_solute, 1, 0);
	  // and an uncoupled one for the solvent...
	  if (last_solvent != last_solute)
	    sim.multibath().add_bath_index(last_solvent, 1, 1);

	}
	// 0 0 1
	else if (ntt[0] == 0 && ntt[1] == 0 && ntt[2] == 1){
	  // the baths
	  sim.multibath().add_bath(temp[2], tau[2]);
	  sim.multibath().add_bath(0, -1);
	  // the atoms in the baths
	  sim.multibath().add_bath_index(last_solute, 1, 1);
	  sim.multibath().add_bath_index(last_solvent, 0, 0);
	}
	// 1 1 0
	else if (ntt[0] == 1 && ntt[1] == 1 && ntt[2] == 0){
	  // the baths
	  sim.multibath().add_bath(temp[0], tau[0]);
	  sim.multibath().add_bath(temp[1], tau[1]);
	  // the atoms in the baths
	  sim.multibath().add_bath_index(last_solute, 0, 1);

	  // and an uncoupled one for the solvent...
	  if (last_solvent != last_solute){
	    sim.multibath().add_bath(0, -1);
	    sim.multibath().add_bath_index(last_solvent, 2, 2);
	  }
	  
	}
	// 1 1 1
	else if (ntt[0] == 1 && ntt[1] == 1 && ntt[2] == 1){
	  // the baths
	  sim.multibath().add_bath(temp[0], tau[0]);
	  sim.multibath().add_bath(temp[1], tau[1]);
	  sim.multibath().add_bath(temp[2], tau[2]);
	  // the atoms in the baths
	  sim.multibath().add_bath_index(last_solute, 0, 1);
	  sim.multibath().add_bath_index(last_solvent, 2, 2);
	}
	// 2 -2 0
	else if (ntt[0] == 2 && ntt[1] == -2 && ntt[2] == 0){
	  // the bath
	  sim.multibath().add_bath(temp[0], tau[0]);
	  // the atoms in the bath
	  sim.multibath().add_bath_index(last_solute, 0, 0);
	  // and an uncoupled one for the solvent...
	  if (last_solvent != last_solute){
	    sim.multibath().add_bath(0, -1);
	    sim.multibath().add_bath_index(last_solvent, 1, 1);
	  }

	}
	// 2 -2 1
	else if (ntt[0] == 2 && ntt[1] == -2 && ntt[2] == 1){
	  // the baths
	  sim.multibath().add_bath(temp[0], tau[0]);
	  sim.multibath().add_bath(temp[2], tau[2]);
	  // the atoms in the baths
	  sim.multibath().add_bath_index(last_solute, 0, 0);
	  sim.multibath().add_bath_index(last_solvent, 1, 1);
	}
	// 3 3 3
	else if (ntt[0] == 3 && ntt[1] == -3 && ntt[2] == -3){
	  // the bath
	  sim.multibath().add_bath(temp[0], tau[0]);
	  // the atoms in the bath
	  sim.multibath().add_bath_index(last_solvent, 0, 0);
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
 * the FLEXCON block.
 */
inline void io::InInput::read_FLEXCON(int &lfcon)
{
  std::vector<std::string> buffer;
  std::vector<std::string>::const_iterator it;
  
  buffer = m_block["FLEXCON"];

  if (!buffer.size()){
    io::messages.add("no FLEXCON block", "InInput",io::messages.notice);
    lfcon = 0;
    return;
  }
  
  it = buffer.begin() + 1;
  _lineStream.clear();
  _lineStream.str(*it);
  
  _lineStream >> lfcon;
  
  if (_lineStream.fail())
    io::messages.add("bad line in FLEXCON block",
		       "InInput", io::message::error);

  if (!_lineStream.eof())
    io::messages.add("end of line not reached in FLEXCON block,"
		     " but should have been: \n" + *it +  "\n",
		     "InInput", io::message::warning);
}

/**
 * read the PRINT
 */
inline void io::InInput::read_PRINT(int &NTPR, 
				    int &NTPL,
				    int &NTPP)
{
  std::vector<std::string> buffer;
  std::vector<std::string>::const_iterator it;
  
  buffer = m_block["PRINT"];

  if (!buffer.size()){
    io::messages.add("no PRINT block", "InInput", io::message::notice);
    NTPR = 1;
    NTPL = 1;
    NTPP = 1;
    return;
  }

  it = buffer.begin() + 1;
  _lineStream.clear();
  _lineStream.str(*it);
  
  _lineStream >> NTPR >> NTPL >> NTPP;

  if (_lineStream.fail())
    io::messages.add("bad line in PRINT block",
		       "InInput", io::message::error);

  if (!_lineStream.eof())
    io::messages.add("End of line not reached in PRINT, but should have been: \n" + *it +  "\n",
		       "InInput", io::message::warning);
    
}

/**
 * read the WRITE
 */
inline void io::InInput::read_WRITE(int &NTWX,
				    int &NTWSE,
				    int &NTWV,
				    int &NTWE,
				    int &NTWG)
{
  std::vector<std::string> buffer;
  std::vector<std::string>::const_iterator it;
  
  buffer = m_block["WRITE"];

  if (!buffer.size()){
    io::messages.add("no WRITE block", "InInput", io::message::notice);
    NTWX = 1;
    NTWSE = 0;
    NTWV = 1;
    NTWE = 1;
    NTWG = 1;
    return;
  }
  int NTWP;
  
  it = buffer.begin() + 1;
  _lineStream.clear();
  _lineStream.str(*it);
  
  _lineStream >> NTWX >> NTWSE >> NTWV >> NTWE >> NTWG >> NTWP;
  
  if (_lineStream.fail())
    io::messages.add("bad line in WRITE block",
		       "InInput", io::message::error);

  if (!_lineStream.eof())
    io::messages.add("End of line not reached in WRITE block, "
		     "but should have been: \n" + *it +  "\n",
		     "InInput", io::message::warning);
    
}

/**
 * read the PCOUPLE block.
 */
inline bool io::InInput::read_PCOUPLE(bool &calc, int &ntp, 
				      math::Matrix &pres0,
				      double &comp, double &tau,
				      interaction::virial_enum &vir)
{
  std::vector<std::string> buffer;

  // read the boundary block
  int ntb, nrdbox;
  read_BOUNDARY(ntb, nrdbox);
  
  // first try for a PCOUPLE03 block
  buffer = m_block["PCOUPLE03"];

  if (buffer.size()){
  
    std::string s1, s2, s3, sblock;
    
    concatenate(buffer.begin()+1, buffer.end()-1, sblock);

    _lineStream.clear();
    _lineStream.str(sblock);
    
    _lineStream >> s1 >> s2 >> comp >> tau >> s3;

    for(int i=0; i<3; ++i){
      for(int j=0; j<3; ++j){
	_lineStream >> pres0(i, j);
      }
    }

    if (_lineStream.fail())
      io::messages.add("bad line in PCOUPLE03 block",
		       "InInput", io::message::error);
    
    std::transform(s1.begin(), s1.end(), s1.begin(), tolower);
    std::transform(s2.begin(), s2.end(), s2.begin(), tolower);
    std::transform(s3.begin(), s3.end(), s3.begin(), tolower);

    bool scale = false;
    
    if (s1 == "off")
      calc = false;
    else if (s1 == "calc")
      calc = true;
    else if (s1 == "scale"){
      calc = true;
      scale = true;
    }
    else{
      io::messages.add("bad value for calc switch in PCOUPLE03 block\n"
		       "(off,calc,scale)",
		       "InInput", io::message::error);
      calc = false;
    }
  
    if (scale){
      if (s2 == "iso")
	ntp = 1;
      else if (s2 == "aniso")
	ntp = 2;
      else if (s2 == "full")
	ntp = 3;
      else{
	io::messages.add("bad value for ntp switch in PCOUPLE03 block\n"
			 "(iso,aniso,full)",
			 "InInput", io::message::error);
	ntp = 0;
      }
    }
    else ntp = 0;
    
    if (calc){
      if (s3 == "none")
	vir = interaction::no_virial;
      else if (s3 == "atomic")
	vir = interaction::atomic_virial;
      else if (s3 == "molecular")
	vir = interaction::molecular_virial;
      else{
	io::messages.add("bad value for virial switch in PCOUPLE03 block\n"
			 "(none,atomic,molecular)",
			 "InInput", io::message::error);
	vir = interaction::no_virial;
      }
    }
    else vir = interaction::no_virial;
    
  } // PCOUPLE03 block
  else{

    buffer = m_block["PCOUPLE"];
    if (!buffer.size()){

      if (abs(ntb) == 2){
	calc = true;
	ntp = 0;
	vir = interaction::molecular_virial;
      }
      else{
	calc = false;
	ntp = 0;
	vir = interaction::no_virial;
      }
      
      for(int i=0; i<3; ++i)
	for(int j=0; j<3; ++j)
	  pres0(i,j) = 0;
      comp = 0;
      tau = -1;

      return false;
    }
  
    _lineStream.clear();
    _lineStream.str(buffer[1]);

    double p0;
    _lineStream >> ntp >> p0 >> comp >> tau;
  
    if (_lineStream.fail())
      io::messages.add("bad line in PCOUPLE block",
		       "InInput", io::message::error);

    for(int i=0; i<3; ++i){
      for(int j=0; j<3; ++j){
	if (i == j)
	  pres0(i,j) = p0;
	else pres0(i,j) = 0.0;
      }
    }

    if (abs(ntb) == 2){
      // pressure calculation
      calc = true;
      vir = interaction::molecular_virial;
    }
    else{
      vir = interaction::no_virial;
      calc = false;
    }
    
  }

  if (ntp && (!calc))
    io::messages.add("pressure coupling activated but "
		     "not calculating pressure",
		     "InInput",
		     io::message::error);


  if (calc && (vir == interaction::no_virial)){
    io::messages.add("PCOUPLE03 block: pressure calculation requested but"
		     " no virial specified!", "InInput",
		     io::message::error);
  }

  if (ntp < 0 || ntp > 3)
    io::messages.add("wrong value for NTP in PCOUPLE block",
		     "InInput",
		     io::message::error);
  
  return true;
  
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
inline void io::InInput::read_PERTURB(int &ntg, double &rlam, double &dlamt, int &nlam, bool &scaling)
{
  std::vector<std::string> buffer;
  std::vector<std::string>::const_iterator it;
  
  // try the new PERTURB03 block
  buffer = m_block["PERTURB03"];
  if (buffer.size()){
    
    std::string b, s1, s2;
    concatenate(buffer.begin()+1, buffer.end()-1, b);

    _lineStream.clear();
    _lineStream.str(b);
  
    _lineStream >> s1 >> rlam >> dlamt >> nlam >> s2;
    
    if (_lineStream.fail())
      io::messages.add("bad line in PERTURB block",
		       "InInput", io::message::error);

    std::transform(s1.begin(), s1.end(), s1.begin(), tolower);
    std::transform(s2.begin(), s2.end(), s2.begin(), tolower);
    
    if (s1 == "on") ntg = 1;
    else if (s1 == "off") ntg = 0;
    else ntg = atoi(s1.c_str());
    
    if (s2 == "on") scaling = true;
    else if (s2 == "off") scaling = false;
    else{
      io::messages.add("wrong value for scaling in PERTURB03 block (on / off): " + s2,
		       "InInput", io::message::error);
      scaling = false;
    }

  }
  else{
    // a common PERTURB block...
    buffer = m_block["PERTURB"];
    if (!buffer.size()){
      // no block at all???
      // std::cerr << "PERTURB block not found" << std::endl;
      ntg = 0;
      rlam = 0;
      dlamt = 0;
      nlam = 1;
      scaling = false;
      return;
    }
  
    it = buffer.begin()+1;
    _lineStream.clear();
    _lineStream.str(*it);
    
    int nrdgl;
    double dmu, dmut;
    double mmu;
    double alpha_lj, alpha_crf;
    
    _lineStream >> ntg >> nrdgl >> rlam >> dlamt >> dmu >> dmut;
    _lineStream.clear();
    _lineStream.str(*(++it));
    _lineStream >> alpha_lj >> alpha_crf >> nlam >> mmu;

    scaling = false;
    
    if (_lineStream.fail())
      io::messages.add("bad line in PERTURB block",
		       "InInput", io::message::error);
    
    if (nrdgl)
      io::messages.add("PERTURB: nrdgl != 0 not allowed",
		       "InInput", io::message::error);

    if (alpha_lj || alpha_crf){
      io::messages.add("PERTURB: softness constants taken from topology!",
		       "InInput", io::message::notice);
    }

  }
  
  if (ntg != 0 && ntg != 1)
    io::messages.add("PERTURB: only ntg = 0 or ntg = 1 allowed",
		     "InInput", io::message::error);
  
  if (nlam<=0){
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

/**
 * read FORCEFIELD block.
 */
inline bool io::InInput::read_FORCEFIELD(int &bond_term,
					 int &angle_term)
{
  std::vector<std::string> buffer;
  std::vector<std::string>::const_iterator it;
  
  buffer = m_block["FORCEFIELD"];

  if (!buffer.size()){
    bond_term = 0;
    angle_term = 0;
    return false;
  }
  
  io::messages.add("using FORCEFIELD block to determine bond term",
		   "InInput", io::message::notice);

  it = buffer.begin()+1;
  _lineStream.clear();
  _lineStream.str(*it);
  
  _lineStream >> bond_term >> angle_term;
  
  if (_lineStream.fail())
    io::messages.add("bad line in FORCEFIELD block",
		     "InInput", io::message::error);

  if (!_lineStream.eof())
    io::messages.add("end of line not reached in FORCEFIELD block,"
		     " but should have been: \n" + *it +  "\n",
		     "InInput", io::message::warning);

  return true;

}

/**
 * read START block.
 */
inline void io::InInput::read_START(int &ntx, int &init, double &tempi, unsigned int &ig)
{
  std::vector<std::string> buffer;
  std::vector<std::string>::const_iterator it;
  
  DEBUG(10, "reading START block");
  buffer = m_block["START"];

  if (!buffer.size()){
    ntx = 1;
    init = 1;
    tempi = 0;
    ig = 0;
    return;
  }
  
  it = buffer.begin()+1;
  _lineStream.clear();
  _lineStream.str(*it);
  
  double heat;
  int ntxo;
  double boltz;
  
  _lineStream >> ntx >> init >> ig >> tempi >> heat >> ntxo >> boltz;
  _lineStream.clear();
  
  if (_lineStream.fail())
    io::messages.add("bad line in START block",
		     "InInput", io::message::error);
  
}
/**
 * read CENTREOFMASS block.
 */
inline bool io::InInput::read_CENTREOFMASS(int &ndfmin, int &ntcm, int &nscm)
{
  std::vector<std::string> buffer;
  std::vector<std::string>::const_iterator it;
  
  DEBUG(10, "reading CENTREOFMASS block");
  buffer = m_block["CENTREOFMASS"];

  if (!buffer.size()){
    ndfmin = 0;
    ntcm = 0;
    nscm = 0;
    return false;
  }
  
  it = buffer.begin()+1;
  _lineStream.clear();
  _lineStream.str(*it);
  
  _lineStream >> ndfmin >> ntcm >> nscm;
  
  _lineStream.clear();
  
  if (_lineStream.fail())
    io::messages.add("bad line in CENTREOFMASS block",
		     "InInput", io::message::error);
  return true;
  
}

/**
 * read PLIST block.
 */
inline bool io::InInput::read_PLIST(bool &grid, int &nsnb, double &rcutp, double &rcutl, double &size)
{
  std::vector<std::string> buffer;
  std::vector<std::string>::const_iterator it;

  DEBUG(10, "pairlist block");
  
  // try a PLIST03
  buffer = m_block["PLIST03"];
  if (buffer.size()){
    // we got a NEW block!!!
    std::string s1, s2;
    
    it = buffer.begin() + 1;
    _lineStream.clear();
    _lineStream.str(*it);
    
    _lineStream >> s1 >> nsnb >> rcutp >> rcutl >> s2;
    
    if (_lineStream.fail())
      throw std::runtime_error("bad line in PLIST03 block");
    
    std::transform(s1.begin(), s1.end(), s1.begin(), tolower);
    std::transform(s2.begin(), s2.end(), s2.begin(), tolower);
    
    if (s1 == "grid") grid = true;
    else if (s1 == "standard") grid = false;
    else{
      io::messages.add("wrong pairlist algorithm chosen (allowed: standard, grid) in PLIST03 block",
		       "InInput", io::message::error);
      grid = false;
      }
    
    if (grid){
      if (s2 == "auto") size = 0.5 * rcutp;
      else{
	size = atof(s2.c_str());
	if (errno){
	  io::messages.add("wrong pairlist grid size chosen (allowed: auto, [size]) in PLIST03 block",
			   "InInput", io::message::error);
	  size = 0.5 * rcutp;
	}
      }
    }
    else size = 0;
  }
  else{
    buffer = m_block["PLIST"];
    
    if (!buffer.size()){
      io::messages.add("no PLIST block in input","InInput",io::message::error);
      size = 0;
      grid = false;
      rcutp = 0;
      rcutl = 0;
      nsnb = 1;
      return false;
    }
    else{
      
      it = buffer.begin() + 1;
      _lineStream.clear();
      _lineStream.str(*it);
	
      int i;
      
      _lineStream >> i >> nsnb >> rcutp >> rcutl;
      size = 0;
      grid = false;
      
      if (_lineStream.fail()){
	io::messages.add("bad line in PLIST block",
			 "InInput", io::message::error);
      }
      
      DEBUG(7, "pairlist update=" << nsnb);
      DEBUG(7, "setting short cutoff=" << rcutp << " long cutoff=" << rcutl);
      
    }
  }
  return true;
}
