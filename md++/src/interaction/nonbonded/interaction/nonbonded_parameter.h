/*
 * This file is part of GROMOS.
 * 
 * Copyright (c) 2011, 2012, 2016, 2018, 2021, 2023 Biomos b.v.
 * See <https://www.gromos.net> for details.
 * 
 * GROMOS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

/**
 * @file nonbonded_parameter.h
 * the base class for the nonbonded parameter.
 */

#ifndef INCLUDED_NONBONDED_PARAMETER_H
#define INCLUDED_NONBONDED_PARAMETER_H

#include "../../interaction_types.h"

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
      : m_lj_parameter(nbp.m_lj_parameter),
	m_cg_parameter(nbp.m_cg_parameter)
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
     * resize the cg_parameter matrix.
     */
    void cg_resize(unsigned int i)
    {
      m_cg_parameter.resize(i);
      std::vector< std::vector<lj_parameter_struct> >::iterator
	it = m_cg_parameter.begin(),
	to = m_cg_parameter.end();
      
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
     * the cg (lj) parameter.
     */
    std::vector<std::vector<lj_parameter_struct> > & cg_parameter()
    {
      return m_cg_parameter;
    }

    /**
     * get the scaling factor for electrostatic 1,4-interactions
     */
    double const & get_coulomb_scaling() {
      return m_coulomb_scaling;
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
     * get the cg parameters for atom type i and j.
     */
    lj_parameter_struct const & cg_parameter(unsigned int iac_i, unsigned int iac_j){
      assert(iac_i < m_cg_parameter.size());
      assert(iac_j < m_cg_parameter[iac_i].size());      
      return m_cg_parameter[iac_i][iac_j];
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

    /**
     * add the cg parameters for atom type i and j.
     */
    void add_cg_parameter(unsigned int iac_i, unsigned int iac_j,
			  lj_parameter_struct lj)
    {
      assert(iac_i < m_cg_parameter.size());
      assert(iac_j < m_cg_parameter.size());
      assert(iac_i < m_cg_parameter[iac_j].size());
      assert(iac_j < m_cg_parameter[iac_i].size());
      
      m_cg_parameter[iac_i][iac_j] = lj;
      m_cg_parameter[iac_j][iac_i] = lj;
    }

    /**
     * set the scaling factor for electrostatic 1,4-interactions
     * e.g. used with Amber Forcefield
     */
    void set_coulomb_scaling(double factor) {
      m_coulomb_scaling = factor;
    }
    
  protected:
    /**
     * the lj parameter.
     */
    std::vector< std::vector<lj_parameter_struct> > m_lj_parameter;

    /**
     * the cg parameter.
     */
    std::vector< std::vector<lj_parameter_struct> > m_cg_parameter;

    /**
    * scaling factor for electrostatic 1,4-interactions
    * e.g. used with Amber Forcefield
    */
    double m_coulomb_scaling;

  };
  
} // interaction

#endif
