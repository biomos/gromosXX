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

// utils_JvalFuncs.cc
#include "JvalFuncs.h"

#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <cassert>
#include <sstream>
#include <vector>
#include <iostream>
#include <string>
#include <vector>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_vector_double.h>

#include "Value.h"
#include "PropertyContainer.h"
#include "../args/Arguments.h"
#include "../gcore/System.h"
#include "../gcore/Molecule.h"
#include "../gmath/Physics.h"
#include "../gromos/Exception.h"

using namespace std;
using namespace gcore;
using namespace gmath;

using gcore::System;
using args::Arguments;
using utils::JvalFuncs;
using utils::PropertyContainer;

// Constructor
JvalFuncs::JvalFuncs(System &sys, Arguments &args) {}

// Destructor
JvalFuncs::~JvalFuncs() {}

// function to read Jval data
void JvalFuncs::read_jval(vector<string> buffer, const System &sys, 
        vector<JvalData::jvalues> &jval, PropertyContainer &props, double &delta) {

  // local (temporary) storage for Karplus parameters
  vector<double> deltas;
  vector<double> As;
  vector<double> Bs;
  vector<double> Cs;

  // read into jvalues class
  for (unsigned int jj = 1; jj < buffer.size() - 1; jj++) {

    // first get atom numbers
    istringstream is(buffer[jj]);
    int i, j, k, l;
    is >> i >> j >> k >> l;

    // adjust atom numbers to gromos numbering
    i--;
    j--;
    k--;
    l--;

    // define local set of Jval parameters
    JvalData::jvalues jv;

    // get molecule number, offset atoms numbers
    int m = 0, offset = 0;
    for (; i >= sys.mol(m).numAtoms() + offset; m++)
      offset += sys.mol(m).numAtoms();
    jv.mol = m;
    jv.i = i - offset;
    jv.j = j - offset;
    jv.k = k - offset;
    jv.l = l - offset;

    // store the angle definitions as properties
    ostringstream os;
    os << "t%" << jv.mol + 1 << ":"
            << jv.i + 1 << "," << jv.j + 1 << ","
            << jv.k + 1 << "," << jv.l + 1;
    props.addSpecifier(os.str());


    // now get the weight and experimental jvalue
    is >> jv.w >> jv.exp;

    // convert jvalue into gromos units: pointless
    //double jval_cf = pico;
    //jv.exp *= jval_cf;

    // the karplus components are not needed but we check that they are the same
    // as svd fitting only makes sense if all data are of the same type
    double A, B, C;
    is >> delta >> A >> B >> C; // discard H
    if (jj > 1 && (delta != deltas[jj-2] || A != As[jj-2] || B != Bs[jj-2] || C != Cs[jj-2] ))
      throw gromos::Exception("svd_fit", "All Jvalues must be of same type\n" + buffer[jj]);
    deltas.push_back(delta); // note that the function will return the last value of delta read
    As.push_back(A);
    Bs.push_back(B);
    Cs.push_back(C);

    if (is.fail())
      throw gromos::Exception("svd_fit", "Bad line in jval-file\n" + buffer[jj]);

    //std::cout << jv.i << " " << jv.j << " " << jv.k << " " << jv.l << " " << jv.exp << " "
    //<< jv.w << " " << A << " " << B << " " << C << " " << delta << std::endl;

    jval.push_back(jv);

  }

}

// function to read weights for individual frames from file
void JvalFuncs::read_weights(vector<string> buffer, vector<JvalWeights::weights> &weight_data) {

  // read into weights vector
  for (unsigned int d = 1; d < buffer.size() - 1; d++) {

    // define local weight struct
    JvalWeights::weights w;
    // get values from file
    istringstream is(buffer[d]);
    is >> w.frame >> w.weight;
    // add to weight_data
    weight_data.push_back(w);

  }
}

// compute the coefficients of the matrix describing bond vector fluctuations
void JvalFuncs::calc_coef(const System &sys, PropertyContainer &fit_props,
        gsl_matrix *coef_mat, int njval, double w, double delta) {

  // compute angles
  fit_props.calc();

  for (int d = 0; d < njval; d++) {

    // calculate the coefficients
    double cosphi = cos((fit_props[d]->getValue().scalar() + delta) * physConst.get_degree2radian());
    double cos2phi = cosphi * cosphi;

    gsl_matrix_set(coef_mat, d, 0,
            gsl_matrix_get(coef_mat, d, 0) + w * cos2phi);
    gsl_matrix_set(coef_mat, d, 1,
            gsl_matrix_get(coef_mat, d, 1) + w * cosphi);
    // last entry is always 1.0
    gsl_matrix_set(coef_mat, d, 2,
            gsl_matrix_get(coef_mat, d, 2) + w * 1.0);

  }
}

// fill a vector with the Jvals
void JvalFuncs::fill_jvalvec(const vector<JvalData::jvalues> &J, gsl_vector *v, int njval)
{
	for (int i=0; i<njval; i++) {
		gsl_vector_set(v,i,J[i].exp);
	}
}

// calculate the Q value (goodness of fit)
double JvalFuncs::calc_Q(gsl_vector *calc, gsl_vector *expt)
{
	int ndat;
	double sdev2,scalc2,tmpf,expt_i,calc_i;
	sdev2 = 0;
	scalc2 = 0;
	ndat = calc->size;
	for (int i=0; i<ndat; i++) {
		calc_i = gsl_vector_get(calc,i);
		expt_i = gsl_vector_get(expt,i);
		tmpf = calc_i-expt_i;
		sdev2 += tmpf*tmpf;
		scalc2 += calc_i*calc_i;
	}
	return sqrt(sdev2/scalc2);
}


