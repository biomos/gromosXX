/**
 * @file OutTrajectory.tcc
 * definition of the OutTrajectory template methods.
 */

template<typename t_simulation>
inline io::OutTrajectory<t_simulation>::OutTrajectory(std::ostream &os,
						      std::ostream &final)
  : m_format(reduced),
    m_old_format(reduced),
    m_pos_traj(&os),
    m_final_traj(&final),
    m_vel_traj(NULL),
    m_force_traj(NULL),
    m_pos(true),
    m_vel(false),
    m_force(false)
{
}

template<typename t_simulation>
inline io::OutTrajectory<t_simulation> & io::OutTrajectory<t_simulation>
::operator<<(t_simulation &sim)
{
  
  if (m_format == reduced){

    _print_timestep(sim, *m_pos_traj);
    _print_positionred(sim.system(), *m_pos_traj);
    _print_box(sim.system(), *m_pos_traj);
    
    if(m_vel){
      _print_timestep(sim, *m_vel_traj);
      _print_velocityred(sim.system(), *m_vel_traj);
    }
    
    if(m_force){
      _print_timestep(sim, *m_force_traj);
      _print_forcered(sim.system(), *m_force_traj);
    }

  }
  else if(m_format == final){
    _print_timestep(sim, *m_final_traj);
    _print_position(sim.system(), sim.topology(), *m_final_traj);
    _print_velocity(sim.system(), sim.topology(), *m_final_traj);
    _print_box(sim.system(), *m_final_traj);
    
    // reset the format after one output (compare std::setw)
    m_format = m_old_format;    
  }
  else{
    _print_timestep(sim, *m_pos_traj);
    _print_position(sim.system(), sim.topology(), *m_pos_traj);
    _print_box(sim.system(), *m_pos_traj);
    
    if(m_vel){
      _print_timestep(sim, *m_vel_traj);
      _print_velocity(sim.system(), sim.topology(), *m_vel_traj);
    }
    
    if(m_force){
      _print_timestep(sim, *m_force_traj);
      _print_force(sim.system(), sim.topology(), *m_force_traj);
    }

    // reset the format after one output (compare std::setw)
    m_format = m_old_format;        
  }
  
  return *this;
}

template<typename t_simulation>
inline io::OutTrajectory<t_simulation> & io::OutTrajectory<t_simulation>
::operator<<(output_format f)
{
  m_old_format = m_format;
  m_format = f;
  return *this;
}

template<typename t_simulation>
inline void io::OutTrajectory<t_simulation>
::velocity_trajectory(std::ostream &os)
{
  m_vel_traj = &os;
  m_vel = true;
}

template<typename t_simulation>
inline void io::OutTrajectory<t_simulation>
::force_trajectory(std::ostream &os)
{
  m_force_traj = &os;
  m_force = true;
}

template<typename t_simulation>
inline void io::OutTrajectory<t_simulation>
::_print_timestep(t_simulation &sim, std::ostream &os)
{
  os.setf(std::ios::fixed, std::ios::floatfield);
  os.precision(9);
  
  os << "TIMESTEP\n"
     << std::setw(15) << sim.steps()
     << std::setw(15) << sim.time()
     << "\nEND\n";
  
}

template<typename t_simulation>
inline void io::OutTrajectory<t_simulation>
::_print_position(simulation::system &sys, simulation::topology &topo, std::ostream &os)
{
  os.setf(std::ios::fixed, std::ios::floatfield);
  os.precision(9);
  
  os << "POSITION\n";
  
  math::VArray &pos = sys.pos();
  simulation::soluteatom &soluteatom = topo.soluteatoms();
  std::vector<std::string> &residue_name = topo.residue_name();

  os << "# first 24 chars ignored\n";
  
  for(int i=0,to = pos.size(); i<to; ++i){

    os << std::setw(6)  << soluteatom(i).residue_nr+1
       << std::setw(5)  << residue_name[soluteatom(i).residue_nr]
       << std::setw(6)  << soluteatom(i).name
       << std::setw(8)  << i+1
       << std::setw(15) << pos(i)(0)
       << std::setw(15) << pos(i)(1)
       << std::setw(15) << pos(i)(2)
       << "\n";
  }
  
  os << "END\n";
  
}

template<typename t_simulation>
inline void io::OutTrajectory<t_simulation>
::_print_positionred(simulation::system &sys, std::ostream &os)
{
  os.setf(std::ios::fixed, std::ios::floatfield);
  os.precision(9);
  
  os << "POSITIONRED\n";
  
  math::VArray &pos = sys.pos();
  
  for(int i=0,to = pos.size(); i<to; ++i){
    if(i% 10 == 0 && i) os << '#' << std::setw(15) << i << "\n";
    
    os << std::setw(15) << pos(i)(0)
       << std::setw(15) << pos(i)(1)
       << std::setw(15) << pos(i)(2)
       << "\n";
  }
  
  os << "END\n";
  
}

template<typename t_simulation>
inline void io::OutTrajectory<t_simulation>
::_print_velocity(simulation::system &sys, simulation::topology &topo, std::ostream &os)
{
  os.setf(std::ios::fixed, std::ios::floatfield);
  os.precision(9);
  
  os << "VELOCITY\n";
  
  math::VArray &vel = sys.vel();
  simulation::soluteatom &soluteatom = topo.soluteatoms();
  std::vector<std::string> &residue_name = topo.residue_name();

  os << "# first 24 chars ignored\n";
  
  for(int i=0,to = vel.size(); i<to; ++i){

    os << std::setw(6)  << soluteatom(i).residue_nr+1
       << std::setw(5)  << residue_name[soluteatom(i).residue_nr]
       << std::setw(6)  << soluteatom(i).name
       << std::setw(8)  << i+1
       << std::setw(15) << vel(i)(0)
       << std::setw(15) << vel(i)(1)
       << std::setw(15) << vel(i)(2)
       << "\n";
  }
  
  os << "END\n";
  
}

template<typename t_simulation>
inline void io::OutTrajectory<t_simulation>
::_print_velocityred(simulation::system &sys, std::ostream &os)
{
  os.setf(std::ios::fixed, std::ios::floatfield);
  os.precision(9);
  
  os << "VELOCITYRED\n";
  
  math::VArray &vel = sys.vel();
  
  for(int i=0,to = vel.size(); i<to; ++i){
    if(i% 10 == 0) os << '#' << std::setw(15) << i << "\n";
    
    os << std::setw(15) << vel(i)(0)
       << std::setw(15) << vel(i)(1)
       << std::setw(15) << vel(i)(2)
       << "\n";
  }
  
  os << "END\n";
  
}

template<typename t_simulation>
inline void io::OutTrajectory<t_simulation>
::_print_force(simulation::system &sys, simulation::topology &topo, std::ostream &os)
{
  os.setf(std::ios::fixed, std::ios::floatfield);
  os.precision(9);
  
  os << "FORCE\n";
  
  math::VArray &force = sys.force();
  simulation::soluteatom &soluteatom = topo.soluteatoms();
  std::vector<std::string> &residue_name = topo.residue_name();

  os << "# first 24 chars ignored\n";
  
  for(int i=0,to = force.size(); i<to; ++i){

    os << std::setw(6)  << soluteatom(i).residue_nr+1
       << std::setw(5)  << residue_name[soluteatom(i).residue_nr]
       << std::setw(6)  << soluteatom(i).name
       << std::setw(8)  << i+1
       << std::setw(15) << force(i)(0)
       << std::setw(15) << force(i)(1)
       << std::setw(15) << force(i)(2)
       << "\n";
  }
  
  os << "END\n";
  
}

template<typename t_simulation>
inline void io::OutTrajectory<t_simulation>
::_print_forcered(simulation::system &sys, std::ostream &os)
{
  os.setf(std::ios::fixed, std::ios::floatfield);
  os.precision(9);
  
  os << "FORCERED\n";
  
  math::VArray &force = sys.force();
  
  for(int i=0,to = force.size(); i<to; ++i){
    if(i% 10 == 0) os << '#' << std::setw(15) << i << "\n";
    
    os << std::setw(15) << force(i)(0)
       << std::setw(15) << force(i)(1)
       << std::setw(15) << force(i)(2)
       << "\n";
  }
  
  os << "END\n";
  
}

template<typename t_simulation>
inline void io::OutTrajectory<t_simulation>
::_print_box(simulation::system &sys, std::ostream &os)
{
  os.setf(std::ios::fixed, std::ios::floatfield);
  os.precision(9);
  
  os << "TRICLINICBOX\n";
  
  math::Matrix &box = sys.box();
  
  os << std::setw(5) << sys.boundary_condition() << "\n";
  
  for(int i=0,to = 3; i<to; ++i){
    
    os << std::setw(15) << box(0)(i)
       << std::setw(15) << box(1)(i)
       << std::setw(15) << box(2)(i)
       << "\n";
  }
  
  os << "END\n";
  
}

