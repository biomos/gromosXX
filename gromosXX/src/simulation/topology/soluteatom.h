/**
 * @file soluteatom.h
 * the soluteatom class.
 */

#ifndef INCLUDED_SOLUTEATOM_H
#define INCLUDED_SOLUTEATOM_H

namespace simulation
{
  /**
   * @class soluteatom
   * holds solute atom information.
   */
  class soluteatom
  {
  public:
    /**
     * @struct soluteatom_struct
     * holds the information.
     */
    struct soluteatom_struct
    {
      std::string name;
      int residue_nr;
      int iac;
    };
    /**
     * accessor.
     */
    soluteatom_struct & operator()(size_t i);
    
    /**
     * add a soluteatom at the end of the list.
     */
    void add(std::string name, int residue_nr, int iac);

  private:
    std::vector<soluteatom_struct> m_information;
    
  };
  
} // simulation

// inline methods
#include "soluteatom.tcc"

#endif
