/**
 * @file OutTrajectory.tcc
 * definition of the OutTrajectory template methods.
 */

#undef MODULE
#undef SUBMODULE
#define MODULE io
#define SUBMODULE trajectory
#include "../../debug.h"

template<typename t_simulation>
inline io::OutTrajectory<t_simulation>::OutTrajectory(std::ostream &os,
						      std::ostream &final, int every)
  : m_format(reduced),
    m_old_format(reduced),
    m_pos_traj(&os),
    m_final_traj(&final),
    m_vel_traj(NULL),
    m_force_traj(NULL),
    m_pos(true),
    m_vel(false),
    m_force(false),
    m_every_pos(every),
    m_every_vel(0),
    m_every_force(0),
    m_precision(9),
    m_force_precision(9),
    m_width(15),
    m_force_width(18)
{
  assert(m_every_pos > 0);
}

template<typename t_simulation>
inline void io::OutTrajectory<t_simulation>::print_title(std::string title)
{
  if (m_pos){
    *m_pos_traj << "TITLE\n"
		<< title << "\n"
		<< "\tposition trajectory\n"
		<< "END\n";
  }
  
  *m_final_traj << "TITLE\n"
		<< title << "\n"
		<< "\tfinal structure\n"
		<< "END\n";

  if (m_vel){
    *m_vel_traj << "TITLE\n"
		<< title << "\n"
		<< "\tvelocity trajectory\n"
		<< "END\n";
  }

  if (m_force){
    *m_force_traj << "TITLE\n"
		  << title << "\n"
		  << "\tforce trajectory\n"
		  << "END\n";
  }

}

template<typename t_simulation>
inline io::OutTrajectory<t_simulation> & io::OutTrajectory<t_simulation>
::operator<<(t_simulation &sim)
{
  
  if (m_format == reduced){

    if((sim.steps() % m_every_pos) == 0){
      _print_timestep(sim, *m_pos_traj);
      _print_positionred(sim.system(), *m_pos_traj);
      if (sim.system().boundary_condition() != math::vacuum)
	_print_box(sim.system(), *m_pos_traj);
    }
    
    if (m_vel && (sim.steps() % m_every_vel) == 0){
      _print_timestep(sim, *m_vel_traj);
      _print_velocityred(sim.system(), *m_vel_traj);
    }
    
    if(m_force && (sim.steps() % m_every_force) == 0){
      if(sim.steps()){
	_print_timestep(sim, *m_force_traj);
	_print_forcered(sim.system(), *m_force_traj);
      }
    }

  }
  else if(m_format == final){
    _print_timestep(sim, *m_final_traj);
    _print_position(sim.system(), sim.topology(), *m_final_traj);
    _print_velocity(sim.system(), sim.topology(), *m_final_traj);
    _print_box(sim.system(), *m_final_traj);
    // forces still go to the force trajectory
    if (m_force){
      _print_timestep(sim, *m_force_traj);
      _print_forcered(sim.system(), *m_force_traj);
    }
    
    // reset the format after one output (compare std::setw)
    m_format = m_old_format;    
  }
  else{

    if((sim.steps() % m_every_pos) == 0){
      _print_timestep(sim, *m_pos_traj);
      _print_position(sim.system(), sim.topology(), *m_pos_traj);
      _print_box(sim.system(), *m_pos_traj);
    }
    
    if (m_vel && (sim.steps() % m_every_vel) == 0){
      _print_timestep(sim, *m_vel_traj);
      _print_velocity(sim.system(), sim.topology(), *m_vel_traj);
    }
    
    if(m_force && (sim.steps() % m_every_force) == 0){
      if (sim.steps()){
	_print_timestep(sim, *m_force_traj);	
	_print_force(sim.system(), sim.topology(), *m_force_traj);
      }
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
::velocity_trajectory(std::ostream &os, int every)
{
  m_vel_traj = &os;
  m_vel = true;
  m_every_vel = every;
  assert(m_every_vel > 0);
}

template<typename t_simulation>
inline void io::OutTrajectory<t_simulation>
::force_trajectory(std::ostream &os, int every)
{
  m_force_traj = &os;
  m_force = true;
  m_every_force = every;
  assert(m_every_force > 0);
}

template<typename t_simulation>
inline void io::OutTrajectory<t_simulation>
::_print_timestep(t_simulation &sim, std::ostream &os)
{
  os.setf(std::ios::fixed, std::ios::floatfield);
  os.precision(m_precision);
  
  os << "TIMESTEP\n"
     << std::setw(m_width) << sim.steps()
     << std::setw(m_width) << sim.time()
     << "\nEND\n";
  
}

template<typename t_simulation>
inline void io::OutTrajectory<t_simulation>
::_print_position(typename t_simulation::system_type &sys,
		  typename t_simulation::topology_type &topo, std::ostream &os)
{
  os.setf(std::ios::fixed, std::ios::floatfield);
  os.precision(m_precision);
  
  os << "POSITION\n";
  
  math::VArray &pos = sys.pos();
  simulation::Solute &solute = topo.solute();
  std::vector<std::string> &residue_name = topo.residue_name();

  math::Vec v;
  
  os << "# first 24 chars ignored\n";
  
  for(int i=0,to = topo.num_solute_atoms(); i<to; ++i){

    v = pos(i);
    sys.periodicity().put_into_box(v);

    os << std::setw(5)  << solute.atom(i).residue_nr+1 << " "
       << std::setw(5)  << std::left << residue_name[solute.atom(i).residue_nr] << " "
       << std::setw(6)  << std::left << solute.atom(i).name << std::right
       << std::setw(6)  << i+1
       << std::setw(m_width) << v(0)
       << std::setw(m_width) << v(1)
       << std::setw(m_width) << v(2)
       << "\n";
  }

  int index = topo.num_solute_atoms();
  int res_nr = 1;

  for(size_t s=0; s < topo.num_solvents(); ++s){

    for(size_t m=0; m < topo.num_solvent_molecules(s); ++m, ++res_nr){
      
      for(size_t a=0; a < topo.solvent(s).num_atoms(); ++a, ++index){
	
	v = pos(index);
	sys.periodicity().put_into_box(v);
	
	os << std::setw(5)  << res_nr
	   << ' ' << std::setw(5)  << std::left
	   << residue_name[topo.solvent(s).atom(a).residue_nr] << " "
	   << std::setw(6)  << std::left << topo.solvent(s).atom(a).name << std::right
	   << std::setw(6)  << index + 1
	   << std::setw(m_width) << v(0)
	   << std::setw(m_width) << v(1)
	   << std::setw(m_width) << v(2)
	   << "\n";
      }
    }
  }

  os << "END\n";
  
}

template<typename t_simulation>
inline void io::OutTrajectory<t_simulation>
::_print_positionred(typename t_simulation::system_type &sys, std::ostream &os)
{
  os.setf(std::ios::fixed, std::ios::floatfield);
  os.precision(m_precision);
  
  os << "POSITIONRED\n";
  
  math::VArray &pos = sys.pos();
  math::Vec v;

  DEBUG(10, "writing POSITIONRED " << pos.size() );
  
  for(int i=0,to = pos.size(); i<to; ++i){

    v = pos(i);
    sys.periodicity().put_into_box(v);

    os << std::setw(m_width) << v(0)
       << std::setw(m_width) << v(1)
       << std::setw(m_width) << v(2)
       << "\n";

    if((i+1)% 10 == 0) os << '#' << std::setw(10) << i+1 << "\n";
  }
  
  os << "END\n";

}

template<typename t_simulation>
inline void io::OutTrajectory<t_simulation>
::_print_velocity(typename t_simulation::system_type &sys,
		  typename t_simulation::topology_type &topo, std::ostream &os)
{
  os.setf(std::ios::fixed, std::ios::floatfield);
  os.precision(m_precision);
  
  os << "VELOCITY\n";
  
  math::VArray &vel = sys.vel();
  simulation::Solute &solute = topo.solute();
  std::vector<std::string> &residue_name = topo.residue_name();

  os << "# first 24 chars ignored\n";
  
  for(int i=0,to = topo.num_solute_atoms(); i<to; ++i){

    os << std::setw(5)  << solute.atom(i).residue_nr+1 << " "
       << std::setw(5)  << std::left << residue_name[solute.atom(i).residue_nr] << " "
       << std::setw(6)  << std::left << solute.atom(i).name << std::right
       << std::setw(6)  << i+1
       << std::setw(m_width) << vel(i)(0)
       << std::setw(m_width) << vel(i)(1)
       << std::setw(m_width) << vel(i)(2)
       << "\n";
  }
  
  int index = topo.num_solute_atoms();
  int res_num = 1;
  
  for(size_t s=0; s < topo.num_solvents(); ++s){

    for(size_t m=0; m < topo.num_solvent_molecules(s); ++m, ++res_num){
      
      for(size_t a=0; a < topo.solvent(s).num_atoms(); ++a, ++index){
	
	os << std::setw(5)  << res_num << " "
	   << std::setw(5)  << std::left 
	   << residue_name[topo.solvent(s).atom(a).residue_nr] << " "
	   << std::setw(6)  << std::left << topo.solvent(s).atom(a).name << std::right
	   << std::setw(6)  << index + 1
	   << std::setw(m_width) << vel(index)(0)
	   << std::setw(m_width) << vel(index)(1)
	   << std::setw(m_width) << vel(index)(2)
	   << "\n";
      }
    }
  }

  os << "END\n";

}

template<typename t_simulation>
inline void io::OutTrajectory<t_simulation>
::_print_velocityred(typename t_simulation::system_type &sys, std::ostream &os)
{
  os.setf(std::ios::fixed, std::ios::floatfield);
  os.precision(m_precision);
  
  os << "VELOCITYRED\n";
  
  math::VArray &vel = sys.vel();
  
  for(int i=0,to = vel.size(); i<to; ++i){

    os << std::setw(m_width) << vel(i)(0)
       << std::setw(m_width) << vel(i)(1)
       << std::setw(m_width) << vel(i)(2)
       << "\n";

    if((i+1)% 10 == 0) os << '#' << std::setw(10) << i+1 << "\n";
  }
  
  os << "END\n";
  
}

template<typename t_simulation>
inline void io::OutTrajectory<t_simulation>
::_print_force(typename t_simulation::system_type &sys,
	       typename t_simulation::topology_type &topo, std::ostream &os)
{
  os.setf(std::ios::fixed, std::ios::floatfield);
  os.precision(m_force_precision);
  
  os << "FORCE\n";
  
  math::VArray &force = sys.force();
  simulation::Solute &solute = topo.solute();
  std::vector<std::string> &residue_name = topo.residue_name();

  os << "# first 24 chars ignored\n";
  
  for(int i=0,to = topo.num_solute_atoms(); i<to; ++i){

    os << std::setw(6)  << solute.atom(i).residue_nr+1
       << std::setw(5)  << residue_name[solute.atom(i).residue_nr]
       << std::setw(6)  << solute.atom(i).name
       << std::setw(8)  << i+1
       << std::setw(m_force_width) << force(i)(0)
       << std::setw(m_force_width) << force(i)(1)
       << std::setw(m_force_width) << force(i)(2)
       << "\n";
  }
  
  int index = topo.num_solute_atoms();
  
  for(size_t s=0; s < topo.num_solvents(); ++s){

    for(size_t m=0; m < topo.num_solvent_molecules(s); ++m){
      
      for(size_t a=0; a < topo.solvent(s).num_atoms(); ++a, ++index){
	
	os << std::setw(6)  << topo.solvent(s).atom(a).residue_nr+1
	   << std::setw(5)  << residue_name[topo.solvent(s).atom(a).residue_nr]
	   << std::setw(6)  << topo.solvent(s).atom(a).name
	   << std::setw(8)  << index + 1
	   << std::setw(m_force_width) << force(index)(0)
	   << std::setw(m_force_width) << force(index)(1)
	   << std::setw(m_force_width) << force(index)(2)
	   << "\n";
      }
    }
  }

  os << "END\n";

}

template<typename t_simulation>
inline void io::OutTrajectory<t_simulation>
::_print_forcered(typename t_simulation::system_type &sys, std::ostream &os)
{
  os.setf(std::ios::fixed, std::ios::floatfield);
  os.precision(m_force_precision);
  
  os << "FORCERED\n";
  
  math::VArray &force = sys.force();
  
  for(int i=0,to = force.size(); i<to; ++i){
    
    os << std::setw(m_force_width) << force(i)(0)
       << std::setw(m_force_width) << force(i)(1)
       << std::setw(m_force_width) << force(i)(2)
       << "\n";

    if((i+1)% 10 == 0) os << '#' << std::setw(10) << i+1 << "\n";
  }
  
  os << "END\n";
  
}

template<typename t_simulation>
inline void io::OutTrajectory<t_simulation>
::_print_box(typename t_simulation::system_type &sys, std::ostream &os)
{
  os.setf(std::ios::fixed, std::ios::floatfield);
  os.precision(m_precision);
  
  os << "TRICLINICBOX\n";
  
  math::Matrix const &box = sys.periodicity().box();
  
  os << std::setw(5) << sys.periodicity().boundary_condition() << "\n";
  
  for(int i=0,to = 3; i<to; ++i){
    
    os << std::setw(m_width) << box(0)(i)
       << std::setw(m_width) << box(1)(i)
       << std::setw(m_width) << box(2)(i)
       << "\n";
  }
  
  os << "END\n";
  
}

template<typename t_simulation>
inline void io::OutTrajectory<t_simulation>
::precision(int prec, int add)
{
  m_precision = prec;
  m_width = prec + add;
}

template<typename t_simulation>
inline void io::OutTrajectory<t_simulation>
::force_precision(int prec, int add)
{
  m_force_precision = prec;
  m_force_width = prec + add;
}

template<typename t_simulation>
inline int io::OutTrajectory<t_simulation>
::precision()
{
  return m_precision;
}

template<typename t_simulation>
inline int io::OutTrajectory<t_simulation>
::force_precision()
{
  return m_force_precision;
}
