/**
 * @file OutG96Trajectory.tcc
 * definition of the OutG96Trajectory template methods.
 */

#undef MODULE
#undef SUBMODULE
#define MODULE io
#define SUBMODULE trajectory
#include "../../debug.h"

template<typename t_simulation>
inline io::OutG96Trajectory<t_simulation>::OutG96Trajectory(std::ostream &os,
						      std::ostream &final,
						      int every)
  : OutTrajectory<t_simulation>(os, final, every)
{
}

/**
 * @override because functions are not virtual...
 */
template<typename t_simulation>
inline io::OutG96Trajectory<t_simulation> & io::OutG96Trajectory<t_simulation>
::operator<<(t_simulation &sim)
{
  
  if (m_format == reduced){

    if((sim.steps() % m_every_pos) == 0){
      _print_timestep(sim, *m_pos_traj);
      _print_positionred(sim.system(), sim.topology(), *m_pos_traj);
      if (sim.system().periodicity().boundary_condition() != math::vacuum)
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
inline void io::OutG96Trajectory<t_simulation>
::_print_position(typename t_simulation::system_type &sys,
		  typename t_simulation::topology_type &topo, std::ostream &os)
{
  os.setf(std::ios::fixed, std::ios::floatfield);
  os.precision(9);
  
  os << "POSITION\n";
  
  math::VArray &pos = sys.pos();
  simulation::Solute &solute = topo.solute();
  std::vector<std::string> &residue_name = topo.residue_name();

  math::Vec v, v_box, trans, r;
  
  os << "# first 24 chars ignored\n";

  // put chargegroups into the box (on the fly)
  simulation::chargegroup_iterator cg_it = topo.chargegroup_begin(),
    cg_to = topo.chargegroup_end();

  // solute chargegroups...
  size_t i = 0;
  for( ; i < topo.num_solute_chargegroups(); ++cg_it, ++i){
    cg_it.cog(pos, v);
    v_box = v;
    sys.periodicity().put_into_positive_box(v_box);
    trans = v_box - v;
    
    // atoms in a chargegroup
    simulation::chargegroup_iterator::atom_iterator at_it = cg_it.begin(),
      at_to = cg_it.end();
    for( ; at_it != at_to; ++at_it){
      r = pos(*at_it) + trans;

      os << std::setw(5)  << solute.atom(*at_it).residue_nr+1 << " "
	 << std::setw(5)  << std::left << residue_name[solute.atom(*at_it).
						       residue_nr] << " "
	 << std::setw(6)  << std::left << solute.atom(*at_it).name 
	 << std::right
	 << std::setw(6)  << *at_it+1
	 << std::setw(15) << r(0)
	 << std::setw(15) << r(1)
	 << std::setw(15) << r(2)
	 << "\n";
    } 
  }

  // solvent chargegroups
  size_t s = 0;
  size_t mol = 0;

  for( ; cg_it != cg_to; ++cg_it, ++mol){
    v = pos(**cg_it);
    v_box = v;
    sys.periodicity().put_into_positive_box(v_box);
    trans = v_box - v;
    
    if (mol >= topo.num_solvent_molecules(s)) ++s;
    
    // loop over the atoms
    simulation::chargegroup_iterator::atom_iterator at_it = cg_it.begin(),
      at_to = cg_it.end();
    // one chargegroup per solvent
    for(size_t atom=0; at_it != at_to; ++at_it, ++atom){
      r = pos(*at_it) + trans;
	
      // std::cout << "atom it: " << *at_it << " atom: " << atom << std::endl;
      
      os << std::setw(5)  << mol+1
	 << ' ' << std::setw(5)  << std::left
	 << residue_name[topo.solvent(s).atom(atom).residue_nr] << " "
	 << std::setw(6)  << std::left << topo.solvent(s).atom(atom).name
	 << std::right
	 << std::setw(6)  << *at_it + 1
	 << std::setw(15) << r(0)
	 << std::setw(15) << r(1)
	 << std::setw(15) << r(2)
	   << "\n";
    }
  }

  os << "END\n";
  
}

template<typename t_simulation>
inline void io::OutG96Trajectory<t_simulation>
::_print_positionred(typename t_simulation::system_type &sys, 
		     typename t_simulation::topology_type &topo, 
		     std::ostream &os)
{
  os.setf(std::ios::fixed, std::ios::floatfield);
  os.precision(9);
  
  os << "POSITIONRED\n";
  
  math::VArray &pos = sys.pos();
  simulation::Solute &solute = topo.solute();
  std::vector<std::string> &residue_name = topo.residue_name();

  math::Vec v, v_box, trans, r;
  
  // put chargegroups into the box (on the fly)
  simulation::chargegroup_iterator cg_it = topo.chargegroup_begin(),
    cg_to = topo.chargegroup_end();

  // solute chargegroups...
  size_t i = 0;
  for( ; i < topo.num_solute_chargegroups(); ++cg_it, ++i){
    cg_it.cog(pos, v);
    v_box = v;
    sys.periodicity().put_into_positive_box(v_box);
    trans = v_box - v;
    
    // atoms in a chargegroup
    simulation::chargegroup_iterator::atom_iterator at_it = cg_it.begin(),
      at_to = cg_it.end();
    for( ; at_it != at_to; ++at_it){
      r = pos(*at_it) + trans;

      os << std::setw(15) << r(0)
	 << std::setw(15) << r(1)
	 << std::setw(15) << r(2)
	 << "\n";

      if((*at_it+1)% 10 == 0) os << '#' << std::setw(10) << *at_it+1 << "\n";

    } 
  }

  // solvent chargegroups
  size_t s = 0;
  size_t mol = 0;

  for(; cg_it != cg_to; ++cg_it, ++mol){
    v = pos(**cg_it);
    v_box = v;
    sys.periodicity().put_into_positive_box(v_box);
    trans = v_box - v;
    
    if (mol >= topo.num_solvent_molecules(s)) ++s;
    
    // loop over the atoms
    simulation::chargegroup_iterator::atom_iterator at_it = cg_it.begin(),
      at_to = cg_it.end();
    for( ; at_it != at_to; ++at_it){
      r = pos(*at_it) + trans;
	
      os << std::setw(15) << r(0)
	 << std::setw(15) << r(1)
	 << std::setw(15) << r(2)
	   << "\n";
      
      if((*at_it+1)% 10 == 0) os << '#' << std::setw(10) << *at_it+1 << "\n";
    }
  }

  os << "END\n";
  
}

template<typename t_simulation>
inline io::OutG96Trajectory<t_simulation> & io::OutG96Trajectory<t_simulation>
::operator<<(output_format f)
{
  m_old_format = m_format;
  m_format = f;
  return *this;
}
