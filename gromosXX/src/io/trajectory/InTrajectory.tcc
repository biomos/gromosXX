/**
 * @file InTrajectory.tcc
 * implements methods of InTrajectory.
 */

#undef MODULE
#undef SUBMODULE

#define MODULE io
#define SUBMODULE trajectory

#include "../../debug.h"

/**
 * Default Constructor.
 */
inline io::InTrajectory::InTrajectory()
  : GInStream(),
    read_position(true),
    read_velocity(true),
    read_box(true),
    read_boxindices(false)
{
}

/**
 * Constructor.
 */
inline io::InTrajectory::InTrajectory(std::istream &is)
  : GInStream(is),
    read_position(true),
    read_velocity(true),
    read_box(true),
    read_boxindices(false)
{
}

/**
 * read in a trajectory.
 * this function is kind of specialized to reading initial position and velocities.
 * if a block that should be read is not present, it will go ahead through the trajectory,
 * if it encounters one of the other blocks a second time, this is considered an error.
 */
template<math::boundary_enum b>
inline io::InTrajectory &io::InTrajectory::operator>>(simulation::System<b>& sys){

  std::vector<std::string> buffer;
  std::vector<std::string>::const_iterator it;

  bool pos_read = false;
  bool vel_read = false;
  bool box_read = false;
  bool boxindices_read = false;
  
  while (true){

    bool required = ((!pos_read && read_position) ||
		     (!vel_read && read_velocity) ||
		     (!box_read && read_box) ||
		     (!boxindices_read && read_boxindices));

    // std::cerr << "required: " << required << std::endl;
    
    try{
      io::getblock(stream(), buffer);
    }
    catch(std::runtime_error e){
      // std::cerr << "no more blocks" << std::endl;
      
      // check whether we already have all required blocks...
      if (!required)
	return *this;
      
      std::cout << "failed to read block from trajectory\n"
		<< e.what() << std::endl;
      std::cout << "position: " << pos_read << " velocity: "
		<< vel_read << " box: " << box_read << std::endl;
      throw;
    }

    // std::cerr << "block " << buffer[0] << std::endl;

    DEBUG(10, "\tblock " << buffer[0]);
    
    if (buffer[0] == "POSITION"){
      if (pos_read){
	// check whether we have all the required ones but the boxindices
	if (required)
	  throw std::runtime_error("second POSITION/POSITIONRED block before "
				   "all other required blocks have been read in");
	else return *this;
      }
      
      DEBUG(10, "resize: " << buffer.size()-2);
      sys.resize(buffer.size()-2);
      pos_read = _read_position(sys.pos(), buffer);
      DEBUG(10, "\tPOSITION read");

    }
    
    if (buffer[0] == "POSITIONRED"){
      if (pos_read){
	if (required)
	  throw std::runtime_error("second POSITION/POSITIONRED block before "
				   "all other required blocks have been read in");
	else return *this;
      }
	  
      sys.resize(buffer.size()-2);
      pos_read = _read_positionred(sys.pos(), buffer);
      DEBUG(10, "\tPOSITIONRED read");
    }

    if (buffer[0] == "VELOCITY"){
      if (vel_read){
	if (required)
	  throw std::runtime_error("second VELOCITY/VELOCITYRED block before "
				   "all other required blocks have been read in");
	else return *this;
      }
      
      sys.resize(buffer.size()-2);
      vel_read = _read_velocity(sys.vel(), buffer);
      sys.old_vel() = sys.vel();

      DEBUG(10, "\tVELOCITY read");
    }
    
    if (buffer[0] == "VELOCITYRED"){
      if (vel_read){
	if (required)
	  throw std::runtime_error("second VELOCITY/VELOCITYRED block before "
				   "all other required blocks have been read in");
	else return *this;
      }
      
      sys.resize(buffer.size()-2);
      vel_read = _read_velocityred(sys.vel(), buffer);
      sys.old_vel() = sys.vel();
      
      DEBUG(10, "\tVELOCITYRED read");
    }

    if (buffer[0] == "TRICLINICBOX"){
      if (box_read){
	if (required)
	  throw std::runtime_error("second TRICLINICBOX block before "
				   "all other required blocks have been read in");
	else return *this;
      }
      
      box_read = _read_box(sys, buffer);
      DEBUG(10, "\tTRICLINICBOX read");
    }

    if (buffer[0] == "BOXINDICES"){
      if (boxindices_read){
	if (required)
	  throw std::runtime_error("second BOXINDICES block before "
				   "all other required blocks have been read in");
	else return *this;
      }
      
      sys.resize(buffer.size()-2);
      boxindices_read = _read_boxindices(sys, buffer);
      DEBUG(10, "\tBOXINDICES read");
    }
	
    
  }

  return *this;
}

inline bool io::InTrajectory::_read_positionred(math::VArray &pos, std::vector<std::string> &buffer)
{
  std::vector<std::string>::const_iterator it = buffer.begin()+1,
    to = buffer.end()-1;
  
  for(int i=0; it != to; ++i, ++it){
    _lineStream.clear();
    _lineStream.str(*it);
    _lineStream >> pos(i)(0) >> pos(i)(1) >> pos(i)(2);
    
    if(_lineStream.fail())
      throw std::runtime_error("bad line in POSITIONRED block");
    if (!_lineStream.eof()) {
      std::string msg = "Warning, end of line not reached, but should have been: \n" + *it +  "\n";
      DEBUG(10, msg);
    }
  }
  
  return true;
  
}

inline bool io::InTrajectory::_read_position(math::VArray &pos, std::vector<std::string> &buffer)
{
  std::vector<std::string>::const_iterator it = buffer.begin()+1,
    to = buffer.end()-1;
  
  std::string s1, s2;
  int n, nr;

  for(int i=0; it != to; ++i, ++it){

    _lineStream.clear();
    _lineStream.str(*it);
    // ignore first 4 fields
    _lineStream >> n >> s1 >> s2 >> nr;
    _lineStream >> pos(i)(0) >> pos(i)(1) >> pos(i)(2);
    
    if(_lineStream.fail()){
      io::messages.add("bad line in POSITION block",
		       "InTrajectory",
		       io::message::critical);
      throw std::runtime_error("bad line in POSITION block");
    }
  }

  return true;
  
}

inline bool io::InTrajectory::_read_velocityred(math::VArray &vel, std::vector<std::string> &buffer)
{
  std::vector<std::string>::const_iterator it = buffer.begin()+1,
    to = buffer.end()-1;
  
  for(int i=0; it != to; ++i, ++it){
    _lineStream.clear();
    _lineStream.str(*it);
    _lineStream >> vel(i)(0) >> vel(i)(1) >> vel(i)(2);
    
    if(_lineStream.fail())
      throw std::runtime_error("bad line in VELOCITYRED block");
    if (!_lineStream.eof()) {
      std::string msg = "Warning, end of line not reached, but should have been: \n" + *it +  "\n";
      DEBUG(10, msg);
    }
  }
  
  return true;
  
}

inline bool io::InTrajectory::_read_velocity(math::VArray &vel, std::vector<std::string> &buffer)
{
  std::vector<std::string>::const_iterator it = buffer.begin()+1,
    to = buffer.end()-1;
  
  std::string s1, s2;
  int n, nr;

  for(int i=0; it != to; ++i, ++it){

    _lineStream.clear();
    _lineStream.str(*it);
    // ignore first 4 fields
    _lineStream >> n >> s1 >> s2 >> nr;
    _lineStream >> vel(i)(0) >> vel(i)(1) >> vel(i)(2);
    
    if(_lineStream.fail()){
      io::messages.add("bad line in VELOCITY block",
		       "InTrajectory",
		       io::message::critical);
      
      throw std::runtime_error("bad line in VELOCITY block");
    }
  }

  return true;
  
}

template<math::boundary_enum b>
inline bool io::InTrajectory
::_read_boxindices(simulation::System<b> &sys,
		   std::vector<std::string> &buffer)
{
  std::vector<typename simulation::System<b>::index_struct> &ind = sys.box_indices();

  std::vector<std::string>::const_iterator it = buffer.begin()+1,
    to = buffer.end()-1;
  
  for(int i=0; it != to; ++i, ++it){

    _lineStream.clear();
    _lineStream.str(*it);
  
    _lineStream >> ind[i].k >> ind[i].l >> ind[i].m;
    
    if(_lineStream.fail()){
      io::messages.add("bad line in BOXINDICES block",
		       "InTrajectory",
		       io::message::critical);
      
      throw std::runtime_error("bad line in BOXINDICES block");
    }
  }

  return true;
}

template<math::boundary_enum b>
inline bool io::InTrajectory::_read_box(simulation::System<b> &sys, std::vector<std::string> &buffer)
{

  std::vector<std::string>::const_iterator it = buffer.begin()+1,
    to = buffer.end()-1;
  
  math::Box box;

  int boundary;
  _lineStream.clear();
  _lineStream.str(*it);
  _lineStream >> boundary;

  ++it;
  
  for(int i=0; it != to; ++i, ++it){

    assert(i<3);
    
    _lineStream.clear();
    _lineStream.str(*it);

    _lineStream >> box(0)(i) >> box(1)(i) >> box(2)(i);
    
    if(_lineStream.fail())
      throw std::runtime_error("bad line in TRICLINICBOX block");
    if (!_lineStream.eof()) {
      std::string msg = "Warning, end of line not reached, but should have been: \n" + *it +  "\n";
      DEBUG(10, msg);
    }
  }

  // set the box...
  sys.periodicity().box(box);
  
  // and the boundary condition...
  switch(boundary){
    case 0:
      sys.periodicity().boundary_condition(math::vacuum);
      io::messages.add("boundary conditions set to VACUUM", "InTrajectory", io::message::notice);
      break;
    case 1:
      sys.periodicity().boundary_condition(math::triclinic);
      io::messages.add("boundary conditions set to TRICLINIC", "InTrajectory", io::message::notice);
      break;
    default:
      throw std::runtime_error("bad boundary conditions in TRICLINICBOX block");
  }

  return true;
  
}
