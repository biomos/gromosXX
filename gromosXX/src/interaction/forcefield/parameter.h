/**
 * @file parameter.h
 * structures that hold the forcefield parameter.
 */

#ifndef INCLUDED_PARAMETER_H
#define INCLUDED_PARAMETER_H

namespace interaction
{
  /**
   * @struct bond_type_struct
   * bond types.
   */
  struct bond_type_struct
  {
    double K;
    double r0;
  };
}

#endif
