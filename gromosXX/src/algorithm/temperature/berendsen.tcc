/**
 * @file temperature/berendsen.tcc
 * methods of the berendsen thermostat
 */

#undef MODULE
#undef SUBMODULE

#define MODULE algorithm
#define SUBMODULE temperature

#include "../../debug.h"

algorithm::Berendsen_Thermostat::Berendsen_Thermostat()
{
}

template<typename t_simulation>
inline void algorithm::Berendsen_Thermostat
::apply(t_simulation &sim, double const dt)
{
  
  // calculate the old kinetic energies
  sim.calculate_mol_ekin(-1);
  
  math::VArray &vel = sim.system().vel();
  simulation::Energy &e = sim.system().energies();
  
  // loop over the baths
  std::vector<simulation::bath_struct>::iterator
    b_it = sim.multibath().begin(),
    b_to = sim.multibath().end();
  
  for(size_t num=0; b_it != b_to; ++b_it, ++num){
    // do temperature coupling for that bath?
    if (b_it->tau != -1){

      assert(e.kinetic_energy.size() > num);
      DEBUG(7, "pre-scale ekin: " << e.kinetic_energy[num]);

      const double free_temp = 2 * e.kinetic_energy[num] / 
	(b_it->dof * math::k_Boltzmann);

      b_it->scale = sqrt(1.0 + dt / b_it->tau *
		       (b_it->temperature / free_temp - 1));
      DEBUG(8, "free T " << free_temp << " dof " << b_it->dof);
      DEBUG(8, "ref. T " << b_it->temperature);
      DEBUG(8, "bath " << num << " scaling: " << b_it->scale);

    } // scaling ?
  } // loop over the baths

  //--------------------------------
  // now we have the scaling factors
  //--------------------------------

  // loop over the ranges
  std::vector<simulation::bath_index_struct>::const_iterator
    r_it = sim.multibath().bath_index().begin(),
    r_to = sim.multibath().bath_index().end();

  int last_atom = -1, last_mol = -1;
  
  for(; r_it != r_to; ++r_it){

    DEBUG(8, "atoms " << last_atom + 1 << " - " << r_it->last_atom
	  << " mols " << last_mol + 1 << " - " << r_it->last_molecule);
    
    // decide whether molecular translational kinetic energy and
    // internal, rotational kinetic energies are jointly coupled
    if (r_it->com_bath == r_it->ir_bath){
      DEBUG(8, "jointly coupled, scaling with "
	    << sim.multibath()[r_it->com_bath].scale);

      for(size_t i=last_atom + 1; i <= r_it->last_atom; ++i){
	vel(i) *= sim.multibath()[r_it->com_bath].scale;
      }

    }
    else{
  
      simulation::Molecule_Iterator 
	m_it = sim.topology().molecule_begin(),
	m_to = sim.topology().molecule_begin();
      
      m_it += last_mol + 1;
      m_to += r_it->last_molecule + 1;

      DEBUG(8, "separately coupled");
      DEBUG(8, "scaling from mol " << last_mol << " to "
	    << r_it->last_molecule);

      math::Vec com_v, ir_v;
      double com_ekin, ekin;
      size_t ir_bath, com_bath;

      // which bathes?
      sim.multibath().in_bath(*(m_it.begin()), com_bath, ir_bath);  

      for( ; m_it != m_to; ++m_it){
	// new molecular translational velocities
	sim.system().
	  molecular_translational_ekin(m_it.begin(), m_it.end(),
				       sim.topology().mass(),
				       com_v, com_ekin, ekin,
				       1);

	simulation::Atom_Iterator start = m_it.begin(),
	  end = m_it.end();
    
	// loop over the atoms in the molecule
	for( ; start != end; ++start){
    
	  assert(unsigned(vel.size()) > *start);
	  ir_v = vel(*start) - com_v;
	  vel(*start) =
	    sim.multibath()[com_bath].scale * com_v +
	    sim.multibath()[ir_bath].scale * ir_v;
	} // loop over atoms
      } // loop over molecules of bath
    } // seperately coupled

    last_atom = r_it->last_atom;
    last_mol = r_it->last_molecule;
  } // loop over baths
  
}

template<typename t_simulation>
inline void algorithm::Berendsen_Thermostat
::calculate_kinetic_energy(t_simulation &sim)
{
  // get rid of previous (temperature scaling) energies
  sim.system().energies().zero(false, true);
  sim.calculate_mol_ekin();
  
}
  
template<typename t_simulation>
inline void algorithm::Berendsen_Thermostat
::calculate_kinetic_energy_lambda_derivative(t_simulation &sim)
{
  math::VArray &vel = sim.system().vel();
  math::VArray const & old_vel = sim.system().old_vel();
  math::SArray const & mass = sim.topology().mass();

  // loop over the baths
  std::vector<simulation::bath_index_struct>::iterator
    it = sim.multibath().bath_index().begin(),
    to = sim.multibath().bath_index().end();
  
  size_t last = 0;
  std::vector<double> &e_kin = sim.system().lambda_energies().kinetic_energy;
  
  assert(e_kin.size() == sim.system().energies().kinetic_energy.size());
  // e_kin.resize(sim.system().energies().kinetic_energy.size());

  for(; it != to; ++it){
    // or just put everything into the first bath...?
    size_t bath = it->com_bath;
    
    // assert(e_kin.size() > bath);
    e_kin[bath] = 0.0;
    for(size_t i=last; i<=it->last_atom; ++i){
      
      if(sim.topology().perturbed_atom()[i]){
	DEBUG(7, "\tpertrubed kinetic energy for " << i << " in bath " << bath);
	DEBUG(7, "\tA_mass: " << sim.topology().perturbed_solute().atoms()[i].A_mass() << " B_mass: " << sim.topology().perturbed_solute().atoms()[i].B_mass());
	// for some reason we take the new velocities here
	e_kin[bath] -=
	  (sim.topology().perturbed_solute().atoms()[i].B_mass() -
	   sim.topology().perturbed_solute().atoms()[i].A_mass()) 
	  * math::dot(vel(i), vel(i));
	DEBUG(7, "\tdE_kin/dl: " << e_kin[bath]);
	
      }
      
    } // atoms in bath
    
    e_kin[bath] *= 0.5;
    last = it->last_atom + 1;
    
  } // baths    
}
