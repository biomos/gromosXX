/**
 * @file solvent.tcc
 * inline methods for solvent topology.
 */

/**
 * add a solvent atom.
 */
inline void simulation::solvent::add_atom(std::string name, int iac, double mass, double charge)
{
  solventatom_struct s;
  s.name = name;
  s.iac = iac;
  s.mass = mass;
  s.charge = charge;

  m_atom.push_back(s);
}

/**
 * accessor - atom i.
 */
inline
simulation::solvent::solventatom_struct &
simulation::solvent::atom(size_t i)
{
  assert(i < m_atom.size());
  return m_atom[i];
}

/**
 * accessor - vector.
 */
inline
std::vector<simulation::solvent::solventatom_struct> &
simulation::solvent::atoms()
{
  return m_atom;
}

/**
 * add a solvent constraint.
 */
inline void simulation::solvent
::add_constraint(int i, int j, double b0)
{
  solventconstraint_struct s;
  s.i = i;
  s.j = j;
  s.b0 = b0;

  m_constraint.push_back(s);
}

/**
 * accessor - constraint i.
 */
inline
simulation::solvent::solventconstraint_struct &
simulation::solvent::constraint(size_t i)
{
  assert(i < m_constraint.size());
  return m_constraint[i];
}

/**
 * accessor - vector.
 */
inline
std::vector<simulation::solvent::solventconstraint_struct> &
simulation::solvent::constraints()
{
  return m_constraint;
}

/**
 * accessor: number of atoms.
 */
inline
size_t simulation::solvent::num_atoms()const
{
  return m_atom.size();
}
