/**
 * @file compound.tcc
 * inline methods of compound
 */

/**
 * Constructor.
 */
inline simulation::compound::compound()
  : m_num_atoms(0)
{
}

/**
 * distance constraints accessor.
 */
inline std::vector<simulation::compound::distance_constraint_struct> &
simulation::compound::distance_constraints()
{
  return m_distance_constraint;
}

/**
 * const distance constraints accessor.
 */
inline std::vector<simulation::compound::distance_constraint_struct> const &
simulation::compound::distance_constraints()const
{
  return m_distance_constraint;
}

/**
 * a single distance constraint.
 */
inline simulation::compound::distance_constraint_struct &
simulation::compound::distance_constraint(size_t i)
{
  return m_distance_constraint[i];
}

/**
 * add a distance constraint.
 */
inline void 
simulation::compound::add_distance_constraint(int i, int j, double b0)
{
  distance_constraint_struct s;
  s.i = i;
  s.j = j;
  s.b0 = b0;
  m_distance_constraint.push_back(s);
}

/**
 * number of atoms in compound.
 */
inline size_t
simulation::compound::num_atoms()const
{
  return m_num_atoms;
}


