/*
 * This file is part of GROMOS.
 * 
 * Copyright (c) 2011, 2012, 2016, 2018, 2021, 2023 Biomos b.v.
 * See <https://www.gromos.net> for details.
 * 
 * GROMOS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

/* 
 * File:   RdcFuncs.h
 * Author: jra, lnw
 *
 * Created on November 24, 2009, 4:09 PM, improved mid-2015
 */

#ifndef RDCFUNCS_H
#define RDCFUNCS_H

#include <iomanip>
#include <vector>
#include <string>
#include <ostream>

#include "../gromos/Exception.h"
#include "Noe.h"

namespace utils {

// struct for storing data to describe one specific RDC
class rdcparam {
  rdcparam();

  public:
  // atoms defining inter-nuclear vector
  unsigned int mol;
  unsigned int i;
  unsigned int j;
  unsigned int k;
  unsigned int l;
  // parameters
  double w;    // weight factor
  double exp;  // experimental rdc
  double gi;   // gyromagnetic ratio of atom i
  double gj;   // gyromagnetic ratio of atom j
  int type;    // (1) backbone NH,
               // (2) CA:C (peptide backbone),
               // (3) C:N (peptide backbone) (NOTE: C from residue i, N from i+1),
               // (4) CA:HA (using pseudo-atom for HA) (NOTE: cannot be GLY), 
               // (5) N:H (ASN side-chain),
               // (6) N:H (GLN side-chain),
               // (7) H:H (ASN side-chain),
               // (8) H:H (GLN side-chain)
  double rij;  // internuclear distance for atoms i and j (if not calcrij)
  double rik;  // internuclear distance for atoms i and k (if not calcrij)
  double dmax; // maximum possible rdc for atoms ij (and ik) (if assuming rij is constant)


  rdcparam(const unsigned _mol, const unsigned _i, const unsigned _j, const unsigned _k, const unsigned _l,
           double _w, double _exp, double _gi, double _gj, int _type, double _rij, double _rik);

  rdcparam(const rdcparam &rdcp) : mol(rdcp.mol), i(rdcp.i), j(rdcp.j),
                                   k(rdcp.k), l(rdcp.l), w(rdcp.w), exp(rdcp.exp), gi(rdcp.gi), gj(rdcp.gj), type(rdcp.type),
                                   rij(rdcp.rij), rik(rdcp.rik), dmax(rdcp.dmax) {}

  rdcparam &operator=(const rdcparam &rdcp) {
    mol = rdcp.mol;
    i = rdcp.i;
    j = rdcp.j;
    k = rdcp.k;
    l = rdcp.l;
    w = rdcp.w;
    exp = rdcp.exp;
    gi = rdcp.gi;
    gj = rdcp.gj;
    type = rdcp.type;
    rij = rdcp.rij;
    rik = rdcp.rik;
    dmax = rdcp.dmax;
    return *this;
  }

  void recalculate_bond_lengths(const gcore::System &sys);

  friend std::ostream& operator<<(std::ostream& s, const rdcparam& rdc){
    s << "{" << rdc.i << ", " << rdc.j << ", " << rdc.k << ", " << rdc.l << ", " 
      << std::setprecision(5) << rdc.w << ", " << std::scientific << rdc.exp;
    s.unsetf(std::ios_base::floatfield);
    s << ", " << rdc.gi << ", " << rdc.gj << ", " << rdc.type <<  ", " << rdc.rij <<  ", " << rdc.rik <<  ", "
      << std::scientific << rdc.dmax << "}";
    return s;
  }
};

typedef std::vector<rdcparam> rdcdata_t;


// function to read in RDC data
rdcdata_t read_rdc(const std::vector<std::string> &buffer, const gcore::System &sys, bool fit);

// function to read in RDC grousp
std::vector<std::vector<unsigned int> > read_groups(const std::vector<std::string> &buffer, const unsigned n_rdc);

// compute the coefficients of the matrix describing bond vector fluctuations for fitting
void calc_coef_fit(const gcore::System &sys, const rdcdata_t &fit_data, double coef_mat[]);

// compute the coefficients of the matrix describing bond vector fluctuations for back-calculation
void calc_coef_bc(const gcore::System &sys, const rdcdata_t &fit_data, double coef_mat_j[], double coef_mat_k[]);

// get inter spin vector, normalised (\vec r / |\vec r|)
gmath::Vec get_inter_spin_vector_normalised(const gcore::System &sys, const rdcdata_t &fit_data, int index);

//  returns {xx, yy, xy, xz, yz}
std::vector<double> lls_fit(const gcore::System &sys, const rdcdata_t &fit_data, double coef_mat[]);
std::vector<double> svd_fit(const gcore::System &sys, const rdcdata_t &fit_data, double coef_mat[]);

void diagonalize_tensor(const std::vector<double> &t5, std::vector<double> &eval, std::vector<std::vector<double> >  &evec);


// calculate Q according to Cornilescu
double calc_Q(const std::vector<double> &calc, const std::vector<double> &expt);
// calculate the R value
double calc_R(const std::vector<double> &calc, const std::vector<double> &expt);
// calculate the RMSD
double calc_RMSD(const std::vector<double> &calc, const std::vector<double> &expt);

}

#endif /* _RDCFUNCS_H */

