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
      : sequence_number(0),
	A_IAC(0),
	A_mass(0),
	A_charge(0),
	B_IAC(0),
	B_mass(0),
	B_charge(0),
	LJ_softcore(0),
	crf_softcore(0)
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
      : sequence_number(JLA),
	A_IAC(IACA),
	A_mass(WMA),
	A_charge(CGA),
	B_IAC(IACB),
	B_mass(WMB),
	B_charge(CGB),
	LJ_softcore(SCLJ),
	crf_softcore(SCC)
    {};
    
    size_t sequence_number;
    size_t A_IAC;
    double A_mass;
    double A_charge;
    size_t B_IAC;
    double B_mass;
    double B_charge;
    double LJ_softcore;
    double crf_softcore;
  };
  
} // simulation

#endif

  
    
      
    
