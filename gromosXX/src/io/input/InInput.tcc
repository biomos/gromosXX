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

    if (_lineStream.fail() || ! _lineStream.eof())
      throw std::runtime_error("bad line in PLIST block");

    DEBUG(7, "setting short cutoff=" << rcutp << " long cutoff=" << rcutl);
    
    sim.nonbonded_update(update_step);
    sim.nonbonded_cutoff_short(rcutp);
    sim.nonbonded_cutoff_long(rcutl);
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
 
  if (_lineStream.fail() || ! _lineStream.eof())
    throw std::runtime_error("bad line in SYSTEM block");
  
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
  
  if (_lineStream.fail() || ! _lineStream.eof())
    throw std::runtime_error("bad line in STEP block");
  
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
  
  if (_lineStream.fail() || ! _lineStream.eof())
    throw std::runtime_error("bad line in SHAKE block");

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
  
  if (_lineStream.fail() || ! _lineStream.eof())
    throw std::runtime_error("bad line in PRINT block");

  buffer = m_block["WRITE"];
  it = buffer.begin() + 1;
  _lineStream.clear();
  _lineStream.str(*it);
  
  int format, selection, energy, free_energy;
  _lineStream >> print_trajectory >> selection >> print_velocity
	      >> energy >> free_energy >> format;
  
  if (_lineStream.fail() || ! _lineStream.eof())
    throw std::runtime_error("bad line in WRITE block");
    
}

/**
 * read FORCE block.
 */
inline void io::InInput::read_FORCE(bool &do_bond, bool &do_angle,
				    bool &do_improper, bool &do_dihedral,
				    bool &do_nonbonded)
{
  std::vector<std::string> buffer;
  std::vector<std::string>::const_iterator it;
  
  buffer = m_block["FORCE"];
  
  it = buffer.begin()+1;
  _lineStream.clear();
  _lineStream.str(*it);
  
  bool bondH, angleH, impH, dihedralH, charge;
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

  if (_lineStream.fail() || ! _lineStream.eof())
    throw std::runtime_error("bad line in FORCE block");

}

