/**
 * @file perturbed_atompair.h
 * a perturbed atom pair.
 */

#ifndef INCLUDED_PERTURBED_ATOMPAIR_H
#define INCLUDED_PERTURBED_ATOMPAIR_H

namespace simulation
{
  /**
   * @class Perturbed_Atompair
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
     * @param ieb atom sequence number.
     * @param jeb atom sequence number.
     * @param ieta interaction mode in state A.
     * @param ietb interaction mode in state B.
     */
    Perturbed_Atompair(int ieb, int jeb, int ieta, int ietb)
      : i(ieb),
        j(jeb),
        A_interaction(ieta),
        B_interaction(ietb)
      {};
    
    /**
     * atom i sequence number.
     */
    int i;
    /**
     * atom j sequence number.
     */
    int j;
    /**
     * interaction in state A:
     * - 0: excluded
     * - 1: 1,4 pair
     * - 2: normal interaction.
     */
    int A_interaction;
    /**
     * interaction in state B:
     * - 0: excluded
     * - 1: 1,4 pair
     * - 2: normal interaction.
     */
    int B_interaction;
  };
  
} // simulation

#endif

  
    
      
    
