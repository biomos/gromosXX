/**
 * @file interaction_types.h
 * structures that hold the forcefield parameter.
 */

#ifndef INCLUDED_INTERACTION_TYPES_H
#define INCLUDED_INTERACTION_TYPES_H

namespace interaction
{
  /**
   * @struct bond_type_struct
   * bond types.
   */
  struct bond_type_struct
  {
    bond_type_struct(double K, double r0) : K(K), r0(r0) {}
    double K;
    double r0;
  };

  /**
   * @struct angle_type_struct
   * angle types.
   */
  struct angle_type_struct
  {
    angle_type_struct(double K, double cos0) : K(K), cos0(cos0) {}
    double K;
    double cos0;
  };
  
  /**
   * @struct improper_dihedral_type_struct
   * improper dihedral types.
   */
  struct improper_dihedral_type_struct
  {
    improper_dihedral_type_struct(double K, double q0) : K(K), q0(q0) {}
    double K;
    double q0;
  };
  
  /**
   * @struct dihedral_type_struct
   * dihedral types.
   */
  struct dihedral_type_struct
  {
    dihedral_type_struct(double K, double pd, int m) : K(K), pd(pd), m(m) {}
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
    lj_parameter_struct(double c6, double c12, double cs6, double cs12)
      : c6(c6), c12(c12), cs6(cs6), cs12(cs12) {}
    lj_parameter_struct()
      : c6(0), c12(0), cs6(0), cs12(0) {}
    double c6;
    double c12;
    double cs6;
    double cs12;
  };

}

#endif
