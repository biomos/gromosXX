/**
 * @file InTrajectory.tcc
 * implements methods of InTrajectory.
 */

/**
 * Constructor.
 */
inline io::InTrajectory::InTrajectory(std::istream &is)
  : GInStream(is),
    read_position(true),
    read_velocity(true),
    read_box(true)
{
}

/**
 * read in a trajectory.
 * this function is kind of specialized to reading initial position and velocities.
 * if a block that should be read is not present, it will go ahead through the trajectory,
 * if it encounters one of the other blocks a second time, this is considered an error.
 */
inline io::InTrajectory &io::InTrajectory::operator>>(simulation::system& sys){

  std::vector<std::string> buffer;
  std::vector<std::string>::const_iterator it;

  bool pos_read = false;
  bool vel_read = false;
  bool box_read = false;
  
  while ((!pos_read && read_position) ||
	 (!vel_read && read_velocity) ||
	 (!box_read && read_box)){
    
    io::getblock(stream(), buffer);
    
    if (buffer[0] == "POSITION"){
      if (pos_read)
	throw std::runtime_error("second POSITION/POSITIONRED block before "
				 "all other required blocks have been read in");
      sys.resize(buffer.size()-2);
      pos_read = _read_position(sys.pos(), buffer);
    }
    
    if (buffer[0] == "POSITIONRED"){
      if (pos_read)
	throw std::runtime_error("second POSITION/POSITIONRED block before "
				 "all other required blocks have been read in");
      sys.resize(buffer.size()-2);
      pos_read = _read_positionred(sys.pos(), buffer);
    }

    if (buffer[0] == "VELOCITY"){
      if (vel_read)
	throw std::runtime_error("second VELOCITY/VELOCITYRED block before "
				 "all other required blocks have been read in");
      sys.resize(buffer.size()-2);
      vel_read = _read_velocity(sys.vel(), buffer);
    }
    
    if (buffer[0] == "VELOCITYRED"){
      if (vel_read)
	throw std::runtime_error("second VELOCITY/VELOCITYRED block before "
				 "all other required blocks have been read in");
      sys.resize(buffer.size()-2);
      vel_read = _read_velocityred(sys.vel(), buffer);
    }

    if (buffer[0] == "TRICLINICBOX"){
      if (box_read)
	throw std::runtime_error("second TRICLINICBOX block before "
				 "all other required blocks have been read in");
      box_read = _read_box(sys, buffer);
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
    
    if(_lineStream.fail() || !_lineStream.eof())
      throw std::runtime_error("bad line in POSITIONRED block");
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
    
    if(_lineStream.fail() || !_lineStream.eof())
      throw std::runtime_error("bad line in POSITION block");
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
    
    if(_lineStream.fail() || !_lineStream.eof())
      throw std::runtime_error("bad line in VELOCITYRED block");
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
    
    if(_lineStream.fail() || !_lineStream.eof())
      throw std::runtime_error("bad line in VELOCITY block");
  }

  return true;
  
}

inline bool io::InTrajectory::_read_box(simulation::system &sys, std::vector<std::string> &buffer)
{
  std::vector<std::string>::const_iterator it = buffer.begin()+1,
    to = buffer.end()-1;
  
  math::Matrix &box = sys.box();

  int i;
  _lineStream.clear();
  _lineStream.str(*it);
  _lineStream >> i;
  switch(i){
    case 0:
      sys.boundary_condition(math::vacuum);
      break;
    case 1:
      sys.boundary_condition(math::triclinic);
      break;
    default:
      throw std::runtime_error("bad boundary conditions in TRICLINICBOX block");
  }

  ++it;
  
  for(int i=0; it != to; ++i, ++it){

    assert(i<3);
    
    _lineStream.clear();
    _lineStream.str(*it);

    _lineStream >> box(0)(i) >> box(1)(i) >> box(2)(i);
    
    if(_lineStream.fail() || !_lineStream.eof())
      throw std::runtime_error("bad line in TRICLINICBOX block");
  }

  return true;
  
}
