/**
 * @file perturbed_atompair.h
 * a perturbed atom pair.
 */

#ifndef INCLUDED_PERTURBED_ATOMPAIR_H
#define INCLUDED_PERTURBED_ATOMPAIR_H

namespace simulation
{
  /**
   * @class Perturbed_Atom
   * holds the perturbed atom information.
   */
  class Perturbed_Atompair
  {
  public:
    /**
     * Default Constructor.
     */
    Perturbed_Atompair()
      : i(0),
        j(0),
	A_interaction(0),
	B_interaction(0)
    {};

    /**
     * Constructor.
     * @param IEB atom sequence number.
     * @param JEB atom sequence number.
     * @param IETA interaction mode in state A.
     * @param IETB interaction mode in state B.
     */
    Perturbed_Atompair(int ieb, int jeb, int ieta, int ietb)
      : i(ieb),
        j(jeb),
        A_interaction(ieta),
        B_interaction(ietb)
      {};
    
    int i;
    int j;
    int A_interaction;
    int B_interaction;
  };
  
} // simulation

#endif

  
    
      
    
