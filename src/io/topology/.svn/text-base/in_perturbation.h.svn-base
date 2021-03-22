/**
 * @file in_perturbation.h
 * read in a perturbation topology file (03 format)
 */

/**
 * @page pttopo perturbation topology format
 * @date 24-10-2008
 *
 * - @ref  title
 * - @ref  scaled
 * - @ref  pertatomparam
 * - @ref  pertatompair
 * - @ref  pertbondstretch
 * - @ref  pertbondstretchh
 * - @ref  pertbondangle
 * - @ref  pertbondangleh
 * - @ref  pertimproperdih
 * - @ref  pertimproperdihh
 * - @ref  pertproperdih
 * - @ref  pertproperdihh
 * - @ref  pertpolparam
 * - @ref  mpertatom


@section scaled SCALEDINTERACTIONS block
Specific (nonbonded) interactions can be scaled dependent on lambda.
@verbatim
SCALEDINTERACTIONS
# number of scaled interactions
  1
# scale nonbonded interaction between atoms of energy group 1 and
# atoms of energy group 2:
# from 0.5 (at lambda = 0) to 0.1 (at lambda = 1)
  1     2      0.5     0.1
END
@endverbatim

@section pertatomparam PETATOMPARAM block
Perturbed atom information block
@verbatim
PERTATOMPARAM
# number of perturbed atoms
   2
#  
#   NR RES NAME IAC(A) MASS(A)  CHARGE(A) IAC(B) MASS(B) CHARGE(B)   ALJ  ACRF
     3   1  CZ1     11  12.011       0.15     12  16.011      0.25   1.0  1.0
    17   1  CB2     14  15.035       0.00     13  10.035      0.10   0.0  0.0
END
@endverbatim

@section pertatompair PERTATOMPAIR block
Perturbed NONBPL atom pair block
@verbatim
PERTATOMPAIR
# number of perturbed atom pairs
   1
#  interaction:
#    0 : excluded
#    1 : normal interaction
#    2 : 1,4 interaction
#
#  NR(I) NR(J) INTERACTION(A) INTERACTION(B)
    2     5    2              1
END
@endverbatim

@section pertbondstretch PERTBONDSTRETCH block 
Perturbed bonds NOT involving H-atoms
@verbatim 
PERTBONDSTRETCH
# number of perturbed bonds
	3
#   atom(i) atom(j) bond_type(A) bond_type(B)
       4       6           15	    26
       6      12           15   	25
       3       8           15   	25
END
@endverbatim

@section pertbondstretchh PERTBONDSTRETCHH block
Perturbed bonds involving H-atoms
@verbatim
PERTBONDSTRETCHH
# number of perturbed hydrogen bonds
        3
#   atom(i) atom(j) bond_type(A) bond_type(B)
       4       6           15           26
       6      12           15           25
       3       8           15           25
END
@endverbatim

@section pertbondangle PERTBONDANGLE block
Perturbed bond angles NOT involving H-atoms
@verbatim
PERTBONDANGLE
# number of perturbed bond angles
    3
#    atom(i) atom(j) atom(k) type(A) type(B)
        2       3       8      26       8
        4       6      12      26       7
        3       8      10      26       7
END
@endverbatim

@section pertbondangleh PERTBONDANGLEH block
Perturbed bond angles involving H-atoms
@verbatim
PERTBONDANGLEH
# number of perturbed hydrogen bond angles
    3
#    atom(i) atom(j) atom(k) type(A) type(B)
        2       3       8      26       8
        4       6      12      26       7
        3       8      10      26       7
END
@endverbatim

@section pertimproperdih PERTIMPROPERDIH block
Perturbed improper (harmonic) dihedrals NOT involving H-atoms
@verbatim
PERTIMPROPERDIH
# number of perturbed improper dihedrals
    2
#    atom(i) atom(j) atom(k) atom(l)  type(A) type(B)
       12      13      10       6        1       2
       18      19      13      16        1       2
END
@endverbatim

@section pertimproperdihh PERTIMPROPERDIHH block 
Perturbed improper (harmonic) dihedrals involving H-atoms
@verbatim
PERTIMPROPERDIHH
# number of perturbed hydrogen improper dihedrals
    2
#    atom(i) atom(j) atom(k) atom(l)  type(A) type(B)
       12      13      10       6        1       2
       18      19      13      16        1       2
END
@endverbatim

@section pertproperdih PERTPROPERDIH block
Perturbed (trigonometric) dihedrals NOT involving H-atoms
@verbatim
PERTPROPERDIH
# number of perturbed dihedrals
    2
#    atom(i) atom(j) atom(k) atom(l)  type(A) type(B)
	6      12      13      14        1      17
       12      13      18      19        4      17
END
@endverbatim

@section pertproperdihh PERTPROPERDIHH block
Perturbed (trigonometric) dihedrals involving H-atoms
@verbatim
PERTPROPERDIHH
# number of perturbed hydrogen dihedrals
    2
#    atom(i) atom(j) atom(k) atom(l)  type(A) type(B)
        6      12      13      14        1      17
       12      13      18      19        4      17
END
@endverbatim

@section pertpolparam PERTPOLPARAM block
Perturbed atomic polarisabilities
@verbatim
PERTPOLPARAM
# number of perturbed polarisable atoms 
# (atoms must also appear in PERTATOMPARAM)
  1
#   NR RES NAME   ALPHA(A) E_0(A)   ALPHA(B) E_0(B)
     1  1   OW    0.00093  80.0     0.00093  80.0
END
@endverbatim

@section mpertatom  MPERTATOM block
Multiple perturbed atom information block
@verbatim
MPERTATOM
# NJLA:   number of perturbed atoms
# NPTB:   number of listed perturbation (i.e. number of
#         perturbation states)
# NJLA    NPTB
  6      2
# JLA:    atom sequence number
# ATNAME: atom name
# IAC[i]: integer atom code in state i
# CG[i]:  charge in state i
# ALJ: soft core van der Waals parameter
# ACRF: soft core electrostatic parameter
# JLA  ATNAME   IAC[1] CG[1]   IAC[2] CG[2]    ALJ   ACRF
# name to identify perturbation
                stateA         stateB
# first water
  1    OW       4     -0.82    19     0         1.0   1.0
  2    HW1      18     0.41    19     0         1.0   1.0
  3    HW2      18     0.41    19     0         1.0   1.0
# second water
  4    OW       19     0       4     -0.82      1.0   1.0
  5    HW1      19     0       18     0.41      1.0   1.0
  6    HW2      19     0       18     0.41      1.0   1.0
END
@endverbatim
 */

#ifndef INCLUDED_IN_PERTURBATION_H
#define INCLUDED_IN_PERTURBATION_H

#include "../instream.h"

namespace io
{
  /**
   * @class In_Perturbation
   * reads in a perturbation topology file (03 version)
   * and parses it into Topology
   * @sa topology::Topology
   */
  class In_Perturbation : public GInStream
  {
  public:
    /**
     * Constructor.
     */
    In_Perturbation(std::istream &is);
    /**
     * parse the topology.
     */
    void read(topology::Topology &topo, simulation::Parameter &param, 
              std::ostream & os = std::cout);
    
  };
  
} // io

#endif
