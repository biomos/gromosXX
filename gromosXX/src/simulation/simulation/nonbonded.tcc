/**
 * @file nonbonded.tcc
 * inline methods for nonbonded simulation parameters
 */

/**
 * constructor
 */
inline
simulation::Nonbonded::Nonbonded():
  m_update(5),
  m_cutoff_short(0.8),
  m_cutoff_long(1.4),
  m_RF_epsilon(1.0),
  m_RF_kappa(0.0),
  m_RF_cutoff(1.4), 
  m_RF_exclusion(true),
  m_RF_constant(0.0)
{
}

/**
 * pairlist update every n steps.
 */
inline void simulation::Nonbonded::update(int const update_step)
{
  m_update = update_step;
}
  
/**
 * accessor pairlist update.
 */
inline int simulation::Nonbonded::update()const
{
  return m_update;
}

/**
 * set short range cutoff.
 */
inline void simulation::Nonbonded::cutoff_short(double const cutoff_short)
{
  m_cutoff_short = cutoff_short;
}

/**
 * get short range cutoff.
 */
inline double simulation::Nonbonded::cutoff_short()const
{
  return m_cutoff_short;
}

/**
 * set long range cutoff.
 */
inline void simulation::Nonbonded::cutoff_long(double const cutoff_long)
{
  m_cutoff_long = cutoff_long;
}

/**
 * get long range cutoff.
 */
inline double simulation::Nonbonded::cutoff_long()const
{
  return m_cutoff_long;
}

/**
 * set reaction field epsilon and kappa, calculate the constant
 */
inline void simulation::Nonbonded::RF_constant(double const epsilon, 
					      double const kappa,
					      double const cutoff)
{
  m_RF_epsilon = epsilon;
  m_RF_kappa   = kappa;
  m_RF_cutoff  = cutoff; 
  m_RF_constant = (2 - 2 * epsilon) * (1 + kappa * cutoff)
                   - epsilon * kappa * kappa * cutoff * cutoff;
  m_RF_constant /= (1 + 2 * epsilon) * (1 + kappa * cutoff)
                   + epsilon * kappa * kappa * cutoff * cutoff;
}

/**
 * set reaction field exclusion flag
 */
inline void simulation::Nonbonded::RF_exclusion(bool const flag)
{
  m_RF_exclusion = flag;
}

/**
 * get reaction field epsilon
 * Crf in the book, II-54
 */
inline
double simulation::Nonbonded::RF_epsilon()const
{
  return m_RF_epsilon;
}

/**
 * get reaction field kappa
 */
inline
double simulation::Nonbonded::RF_kappa()const
{
  return m_RF_kappa;
}

/**
 * get reaction field constant Crf
 */
inline
double simulation::Nonbonded::RF_constant()const
{
  return m_RF_constant;
}

/**
 * get reaction field cutoff
 */
inline
double simulation::Nonbonded::RF_cutoff()const
{
  return m_RF_cutoff;
}

/**
 * get reaction field exclusion flag
 */
inline
bool simulation::Nonbonded::RF_exclusion()const
{
  return m_RF_exclusion;
}

/**
 * check state
 */
inline
int simulation::Nonbonded::check_state()const
{
  int result = 0;

  // check if the box is large enough...
  return result;
}
