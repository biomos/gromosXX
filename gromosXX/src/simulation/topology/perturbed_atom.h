/**
 * @file perturbed_atom.h
 * a perturbed atom.
 */

#ifndef INCLUDED_PERTURBED_ATOM_H
#define INCLUDED_PERTURBED_ATOM_H

namespace simulation
{
  /**
   * @class Perturbed_Atom
   * holds the perturbed atom information.
   */
  class Perturbed_Atom
  {
  public:
    /**
     * Default Constructor.
     */
    Perturbed_Atom()
      : m_sequence_number(0),
	m_A_IAC(0),
	m_A_mass(0),
	m_A_charge(0),
	m_B_IAC(0),
	m_B_mass(0),
	m_B_charge(0),
	m_LJ_softcore(0),
        m_crf_softcore(0),
        m_exclusion(),
        m_one_four_pair()
    {};

    /**
     * Constructor.
     * @param JLA atom sequence number.
     * @param IACA integer atom code of state A.
     * @param WMA mass of state A.
     * @param CGA charge of state A.
     * @param IACB integer atom code of state B.
     * @param WMB mass of state B.
     * @param CGB charge of state B.
     * @param SCLJ soft core van der Waals parameter.
     * @param SCC soft core electrostatic parameter.
     */
    Perturbed_Atom(size_t const JLA, size_t const IACA,
		   double const WMA, double const CGA,
		   size_t const IACB, double const WMB,
		   double const CGB,
		   double const SCLJ, double const SCC)
      : m_sequence_number(JLA),
	m_A_IAC(IACA),
	m_A_mass(WMA),
	m_A_charge(CGA),
	m_B_IAC(IACB),
	m_B_mass(WMB),
	m_B_charge(CGB),
	m_LJ_softcore(SCLJ),
        m_crf_softcore(SCC),
        m_exclusion(),
        m_one_four_pair()
    {};
  public:
    size_t sequence_number();
    size_t const sequence_number()const;
    void sequence_number(const size_t);
    size_t A_IAC();
    size_t const A_IAC()const;
    void A_IAC(const size_t);
    
    double A_mass();
    double const A_mass()const;
    void A_mass(const double);
    
    double A_charge();
    double const A_charge()const;
    void A_charge(const double);
    
    size_t B_IAC();
    size_t const B_IAC()const;
    void B_IAC(const size_t);
    
    double B_mass();
    double const B_mass()const;
    void B_mass(const double);
    
    double B_charge();
    double const B_charge()const;
    void B_charge(const double);
    
    double LJ_softcore();
    double const LJ_softcore()const;
    void LJ_softcore(const double);
    
    double crf_softcore();
    double const crf_softcore()const;
    void crf_softcore(const double);
    
    std::set<int> & exclusion();
    std::set<int> const & exclusion()const;
    
    std::set<int> & one_four_pair();
    std::set<int> const & one_four_pair()const;
    
  private:
    size_t m_sequence_number;
    size_t m_A_IAC;
    double m_A_mass;
    double m_A_charge;
    size_t m_B_IAC;
    double m_B_mass;
    double m_B_charge;
    double m_LJ_softcore;
    double m_crf_softcore;
    std::set<int> m_exclusion;
    std::set<int> m_one_four_pair;
  };
  
} // simulation

#include "perturbed_atom.tcc"
#endif

  
    
      
    
