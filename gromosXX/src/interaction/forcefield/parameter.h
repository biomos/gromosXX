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

  /**
   * @struct angle_type_struct
   * angle types.
   */
  struct angle_type_struct
  {
    double K;
    double cos0;
  };
  
  /**
   * @struct improper_dihedral_type_struct
   * improper dihedral types.
   */
  struct improper_dihedral_type_struct
  {
    double K;
    double q0;
  };
  
  /**
   * @struct dihedral_type_struct
   * dihedral types.
   */
  struct dihedral_type_struct
  {
    double K;
    double pd;
    int m;
  };
  
  /**
   * @struct lj_parameter_struct
   * Lennard Jones interaction parameter.
   */
  struct lj_parameter_struct
  {
    double c6;
    double c12;
    double cs6;
    double cs12;
  };

}

#endif
