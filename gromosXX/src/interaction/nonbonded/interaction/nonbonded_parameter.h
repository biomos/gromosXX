/**
 * @file nonbonded_parameter.h
 * the base class for the nonbonded parameter.
 */

#ifndef INCLUDED_NONBONDED_PARAMETER_H
#define INCLUDED_NONBONDED_PARAMETER_H

namespace interaction
{
  /**
   * @class Nonbonded_Parameter
   * base class of the nonbonded interaction storing the parameter.
   */
  class Nonbonded_Parameter
  {
  public:
    /**
     * Copy constructor
     */
    Nonbonded_Parameter(Nonbonded_Parameter const & nbp)
      : m_lj_parameter(nbp.m_lj_parameter)
    {
    }

    /**
     * Constructor.
     */
    Nonbonded_Parameter(){};
    
    /**
     * resize the lj_parameter matrix.
     */
    void resize(unsigned int i)
    {
      m_lj_parameter.resize(i);
      std::vector< std::vector<lj_parameter_struct> >::iterator
	it = m_lj_parameter.begin(),
	to = m_lj_parameter.end();
      
      for(; it!=to; ++it)
	it->resize(i);
    }

    /**
     * the lj parameter.
     */
    std::vector<std::vector<lj_parameter_struct> > & lj_parameter()
    {
      return m_lj_parameter;
    }

    /**
     * get the lj parameters for atom type i and j.
     */
    lj_parameter_struct const & lj_parameter(unsigned int iac_i, unsigned int iac_j){
      assert(iac_i < m_lj_parameter.size());
      assert(iac_j < m_lj_parameter[iac_i].size());      
      return m_lj_parameter[iac_i][iac_j];
    }
    
    /**
     * add the lj parameters for atom type i and j.
     */
    void add_lj_parameter(unsigned int iac_i, unsigned int iac_j,
			  lj_parameter_struct lj)
    {
      assert(iac_i < m_lj_parameter.size());
      assert(iac_j < m_lj_parameter.size());
      assert(iac_i < m_lj_parameter[iac_j].size());
      assert(iac_j < m_lj_parameter[iac_i].size());
      
      m_lj_parameter[iac_i][iac_j] = lj;
      m_lj_parameter[iac_j][iac_i] = lj;
    }
    
  protected:
    /**
     * the lj parameter.
     */
    std::vector< std::vector<lj_parameter_struct> > m_lj_parameter;

  };
  
} // interaction

#endif
