/**
 * @file nonbonded.h
 * the nonbonded simulation parameters class.
 */

#ifndef INCLUDED_NONBONDED_H
#define INCLUDED_NONBONDED_H

namespace simulation
{
  /**
   * @class Nonbonded
   * holds nonbonded information.
   */
  class Nonbonded
  {
  public:
    /**
     * constructor
     */
    Nonbonded();
    
    /**
     * set pairlist update.
     */
    void update(int const update_step);

    /**
     * accessor pairlist update.
     */
    int update()const;
    /**
     * set short range cutoff.
     */
    void cutoff_short(double const cutoff_short);
    /**
     * get short range cutoff.
     */
    double cutoff_short()const;
    /**
     * set long range cutoff.
     */
    void cutoff_long(double const cutoff_long);
    /**
     * get long range cutoff.
     */
    double cutoff_long()const;
    /**
     * set reaction field epsilon and kappa, calculate the constant
     */
    void RF_constant(double const epsilon, 
		     double const kappa, 
		     double const cutoff);
    /**
     * get reaction field epsilon
     */
    double RF_epsilon()const;
    /**
     * get reaction field kappa
     */
    double RF_kappa()const;
    /**
     * get reaction field constant Crf
     */
    double RF_constant()const;
    /**
     * get reaction field cutoff
     */
    double RF_cutoff()const;
    
  private:      
    /**
     * nonbonded update.
     */
    int m_update;
    
    /**
     * nonbonded short range cutoff.
     */
    double m_cutoff_short;
    
    /**
     * nonbonded long range cutoff.
     */
    double m_cutoff_long;
 
    /**
     * nonbonded reaction field epsilon
     */
    double m_RF_epsilon;
    
    /**
     * nonbonded reaction field kappa
     */
    double m_RF_kappa;
    
    /**
     * nonbonded reaction field cutoff
     */
    double m_RF_cutoff;
    
    /**
     * nonbonded reaction field constant
     */
    double m_RF_constant;
  };
	  
  
} // simulation

#include "nonbonded.tcc"

#endif
