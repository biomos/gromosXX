/**
 * @file perturbed_solute.tcc
 * inline methods of Perturbed_Solute
 */

#undef MODULE
#undef SUBMODULE
#define MODULE simulation
#define SUBMODULE topology

#include "../../debug.h"

inline std::vector<simulation::Perturbed_Bond> &
simulation::Perturbed_Solute::bonds()
{
  return m_bond;
}

inline std::vector<simulation::Perturbed_Bond> const &
simulation::Perturbed_Solute::bonds()const
{
  return m_bond;
}

inline std::vector<simulation::Perturbed_Angle> &
simulation::Perturbed_Solute::angles()
{
  return m_angle;
}

inline std::vector<simulation::Perturbed_Angle> const &
simulation::Perturbed_Solute::angles()const
{
  return m_angle;
}

inline std::vector<simulation::Perturbed_Improper_Dihedral> const &
simulation::Perturbed_Solute::improper_dihedrals()const
{
  return m_improper_dihedral;
}

inline std::vector<simulation::Perturbed_Improper_Dihedral> &
simulation::Perturbed_Solute::improper_dihedrals()
{
  return m_improper_dihedral;
}

inline std::vector<simulation::Perturbed_Dihedral> const &
simulation::Perturbed_Solute::dihedrals()const
{
  return m_dihedral;
}
inline std::vector<simulation::Perturbed_Dihedral> &
simulation::Perturbed_Solute::dihedrals()
{
  return m_dihedral;
}

inline std::map<size_t, simulation::Perturbed_Atom> &
simulation::Perturbed_Solute::atoms()
{
  return m_atom;
}

inline std::map<size_t, simulation::Perturbed_Atom> const &
simulation::Perturbed_Solute::atoms()const
{
  return m_atom;
}
inline simulation::Perturbed_Atom & 
simulation::Perturbed_Solute::atom(const size_t i)
{
  return m_atom[i];
}

 
inline std::vector<simulation::Perturbed_Atompair> &
simulation::Perturbed_Solute::atompairs()
{
  return m_atompair;
}

inline std::vector<simulation::Perturbed_Atompair> const &
simulation::Perturbed_Solute::atompairs()const
{
  return m_atompair;
}

inline std::vector<simulation::Perturbed_Solute::perturbed_distance_constraint_struct> const &
simulation::Perturbed_Solute::distance_constraints()const
{
  return m_distance_constraint;
}

inline std::vector<simulation::Perturbed_Solute::perturbed_distance_constraint_struct> &
simulation::Perturbed_Solute::distance_constraints()
{
  return m_distance_constraint;
}

inline void
simulation::Perturbed_Solute::
add_distance_constraint(int const i, int const j, 
			double const A_b0, double const B_b0)
{
  perturbed_distance_constraint_struct s;
  s.i = i;
  s.j = j;
  s.A_b0 = A_b0;
  s.B_b0 = B_b0;
  s.b0 = A_b0;
  m_distance_constraint.push_back(s);
}

inline void
simulation::Perturbed_Solute::
set_distance_constraints(double const lambda)
{
  std::vector<perturbed_distance_constraint_struct>::iterator
    it = m_distance_constraint.begin(),
    to = m_distance_constraint.end();
  
  DEBUG(7, "perturbed distance constraints:");
  for( ; it!=to; ++it){
    it->b0 = (1-lambda)*it->A_b0 + lambda*it->B_b0;
    DEBUG(7, "\tdistance: " << it->b0);
  }
}

