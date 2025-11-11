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

/**
 * @file matrix_overlap.cc
 * calculates overlap between two matrices
 */

/**
 * @page programs Program Documentation
 *
 * @anchor matrix_overlap
 * @section matrix_overlap calculates the overlap between two matrices
 * @author @ref bh
 * @date 21-10-2010
 *
 * This program makes use of the GSL library to calculate the overlap
 * between two matrices. Considering the following equation for the
 * difference between matrices M1 and M2
 *
 * @f[ D = \sqrt{tr((\sqrt{M1}-\sqrt{M2})^2)} @f]
 *
 * where tr is the trace and the square root operator corresponds to the matrix square root and is
 * calculated according to the following steps. Consider a matrix A that can be 
 * diagonalized by
 *
 * @f[ T = V^{-1} A V @f]
 *
 * The square root of the elements of the diagonal matrix is taken. This procedure
 * corresponds to the application of a normal scalar square root operator
 * to all the elements of the diagonal matrix. Third, the square root of the
 * diagonal matrix is used to calculate the square root of the matrix as
 *
 * @f[ A^{1/2} = V T^{1/2} V^{-1} @f]
 *
 * The (normalized) overlap is given by 1 minus the difference D of the matrices
 * divided by the normalization factor as shown below  
 *
 * @f[ O = 1 - \frac{D}{\sqrt{tr(M1)+tr(M2)}} @f]
 *
 *
 *
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@m1</td><td>&lt;matrix 1 &gt; </td></tr>
 * <tr><td> \@m2</td><td>&lt;matrix 2 &gt; </td></tr>
 * <tr><td> \@dim</td><td>&lt;dimension &gt; </td></tr>
 * 
 * </table>
 *
 *
 * Example:
 * @verbatim
  matrix_overlap 
    @m1             matrix1
    @m2             matrix2
    @dim            dimension
 
 @endverbatim
 *
 * <hr>
 */

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>

#include <gsl/gsl_cblas.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_vector_double.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/gromos/Exception.h"

using namespace std;
using namespace gcore;
using namespace bound;
using namespace args;
using namespace gmath;

int readin(string filename, stringstream &input);

int main(int argc, char** argv) {
  Argument_List knowns;
  knowns << "m1" << "m2" << "dim";
       

  string usage = "# " + string(argv[0]);
  usage += "\n\t@m1           <matrix1>\n";
  usage += "\t@m2             <matrix2>\n";
  usage += "\t@dim             <integer: dimension>\n";

  try {
    Arguments args(argc, argv, knowns, usage);

    int dim = args.getValue<int>("dim");
    if(dim <= 0)
      throw gromos::Exception("overlap_edyn", "dim must be bigger than zero");


    std::string m1;
    if(args.count("m1") > 0)
      m1 = args["m1"]/*.c_str()*/;

    std::string m2;
    if(args.count("m2") > 0)
      m2 = args["m2"]/*.c_str()*/;

    double mtx1[dim][dim];
    //ifstream infile(m1.c_str());
    stringstream infile;
    readin(m1.c_str(), infile);
    //infile.open(m1);

    double d;
    for(int i = 0; i < dim; i++) {
      for(int j = 0; j < dim; j++) {
        infile >> d;
        //infile >> mtx1[i][j];
        mtx1[i][j] = d;
      }
    }
    //infile.close();

    double mtx2[dim][dim];
    //ifstream infile2(m2.c_str());
    stringstream infile2;
    readin(m2.c_str(), infile2);
    //infile.open(m2);

    for(int i = 0; i < dim; i++) {
      for(int j = 0; j < dim; j++) {
        infile2 >> mtx2[i][j];
      }
    }
    //infile2.close();
/*
    cout << endl << "MATRIX 1:" << endl;
    for(int i=0; i<dim; i++){
      for(int j=0; j<dim; j++){
        cout << " " << mtx1[i][j]<< "  ";
      }
      cout<<endl;
    }
  */
    //Calculating trace of the cov matrices
    double trace_cov1=0;
    for (int i=0;i<dim;++i){
      trace_cov1+=mtx1[i][i];
    }
    if ( trace_cov1< 0){
      throw gromos::Exception("overlap_edyn",
	    " trace of the covariance matrix 1 is negative. "
	    "This will make problems with the rest of the calculations.\n");
    }

    double trace_cov2=0;
    for (int i=0;i<dim;++i){
      trace_cov2+=mtx2[i][i];
    }
    if ( trace_cov2< 0){
      throw gromos::Exception("overlap_edyn",
	    " trace of the covariance matrix 2 is negative. "
	    "This will make problems with the rest of the calculations.\n");
    }

    // TEST
    cout << "TRACE M1:  " << trace_cov1 <<endl;
    cout << "TRACE M2:  " << trace_cov2 <<endl;

    //STARTING TO PLAY WITH GSL:
    int i, j;
    gsl_matrix * g_m1 = gsl_matrix_alloc(dim, dim);
    gsl_matrix * g_m2 = gsl_matrix_alloc(dim, dim);
    gsl_matrix_set_zero(g_m1);
    gsl_matrix_set_zero(g_m2);
    for(i = 0; i < dim; i++) {
      for(j = 0; j < dim; j++) {
        gsl_matrix_set(g_m1, i, j, mtx1[i][j]);
        gsl_matrix_set(g_m2, i, j, mtx2[i][j]);
      }
    }

    // TEST: Printing M1 and M2:
    /*
    cout << "Matrix 1:" << endl;
    for(i = 0; i < dim; i++){
      for(j = 0; j < dim; j++) {
        double bla = gsl_matrix_get(g_m1, i, j);
        cout << bla << "  ";
      }
      cout << endl;
    }

    cout << "Matrix 2:"<<endl;
    for(i = 0; i < dim; i++){
      for(j = 0; j < dim; j++) {
        double bla = gsl_matrix_get(g_m2, i, j);
        cout << bla << "  ";
      }
      cout << endl;
    }
     */



    //M1

    gsl_vector *eval1 = gsl_vector_alloc (dim);
    gsl_matrix *evec1 = gsl_matrix_alloc (dim, dim);
    gsl_eigen_symmv_workspace * w1 = gsl_eigen_symmv_alloc (dim);
    gsl_eigen_symmv (g_m1, eval1, evec1, w1);
    gsl_eigen_symmv_free(w1);

    //M2

    gsl_vector *eval2 = gsl_vector_alloc (dim);
    gsl_matrix *evec2 = gsl_matrix_alloc (dim, dim);
    gsl_eigen_symmv_workspace * w2 = gsl_eigen_symmv_alloc (dim);
    gsl_eigen_symmv (g_m2, eval2, evec2, w2);
    gsl_eigen_symmv_free(w2);

    //PRINTING EIGENVALS AND EIGENVECS
    /*
    cout << "Eigenvalues and eigenvectors of Matrix 1:"<<endl;
    for(i = 0; i < dim; i++){
      double val = gsl_vector_get(eval1, i);
      cout << val << "  ";
      for(j = 0; j < dim; j++) {
        double bla = gsl_matrix_get(evec1, i, j);
        cout << bla << "  ";
      }
      cout << endl;
    }

    cout << "Eigenvalues and eigenvectors of Matrix 2:"<<endl;
    for(i = 0; i < dim; i++){
      double val = gsl_vector_get(eval2, i);
      cout << val << "  ";
      for(j = 0; j < dim; j++) {
        double bla = gsl_matrix_get(evec2, i, j);
        cout << bla << "  ";
      }
      cout << endl;
    }
     */

    // Construct diag matrix
    gsl_matrix *diag1 = gsl_matrix_alloc (dim, dim);
    gsl_matrix *diag2 = gsl_matrix_alloc (dim, dim);
    gsl_matrix_set_zero (diag1);
    gsl_matrix_set_zero (diag2);

    for(i = 0; i < dim; i++){
      gsl_matrix_set(diag1, i, i, sqrt(abs(gsl_vector_get(eval1, i))));
      gsl_matrix_set(diag2, i, i, sqrt(abs(gsl_vector_get(eval2, i))));
    }

    /*
    for(i = 0; i < dim; i++){
      for(j = 0; j < dim; j++) {
        double bla = gsl_matrix_get(diag1, i, j);
        cout << bla << "  ";
      }
      cout << endl;
    }
    cout << endl << "print eval1";
    for(i = 0; i < dim; i++){
      double bla = gsl_vector_get(eval1, i);
      cout << bla << endl;
    }*/

    // Now that I have sqrt(D) I can calculate the sqrt(M)

    gsl_matrix *root1 = gsl_matrix_alloc (dim, dim);
    gsl_matrix *root2 = gsl_matrix_alloc (dim, dim);
    gsl_matrix_set_zero (root1);
    gsl_matrix_set_zero (root2);

    gsl_matrix *temp1 = gsl_matrix_alloc (dim, dim);
    gsl_matrix *temp2 = gsl_matrix_alloc (dim, dim);
    gsl_matrix_set_zero (temp1);
    gsl_matrix_set_zero (temp2);

    // Calculating the inverse of the (orthogonal) eigenvectors matrix
    gsl_matrix *inv1 = gsl_matrix_alloc (dim, dim);
    gsl_matrix *inv2 = gsl_matrix_alloc (dim, dim);
    gsl_matrix_set_zero (inv1);
    gsl_matrix_set_zero (inv2);
    gsl_matrix_transpose_memcpy(inv1, evec1);
    gsl_matrix_transpose_memcpy(inv2, evec2);

    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
                  1.0, evec1, diag1,
                  0.0, temp1);
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
                  1.0, evec2, diag2,
                  0.0, temp2);

    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
                  1.0, temp1, inv1,
                  0.0, root1);
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
                  1.0, temp2, inv2,
                  0.0, root2);

    gsl_matrix_free (temp1);
    gsl_matrix_free (inv1);
    gsl_matrix_free (evec1);
    gsl_matrix_free (temp2);
    gsl_matrix_free (inv2);
    gsl_matrix_free (evec2);

    gsl_matrix_sub (root1, root2);
    gsl_matrix *square = gsl_matrix_alloc (dim, dim);
    gsl_matrix_set_zero (square);

    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
                  1.0, root1, root1,
                  0.0, square);

    double trace=0;
    for(int i=0; i<dim; i++){
      trace += gsl_matrix_get(square, i, i);
    }

    double diff = sqrt(trace);

    //cout << "TRACE  " << trace;
    cout << "ABSOLUTE DIFFERENCE: "<< diff<<endl;

    cout << "Normalizing..." << endl;

    double root_trace = sqrt(trace_cov1+trace_cov2);
    double norm_overlap = 1-(diff/root_trace);
    cout << "NORMALIZED OVERLAP: "<< norm_overlap<<endl;

    gsl_matrix_free(root1);
    gsl_matrix_free(root2);
    gsl_matrix_free(square);
    gsl_matrix_free(diag1);
    gsl_matrix_free(diag2);
    gsl_matrix_free(g_m1);
    gsl_matrix_free(g_m2);
    gsl_vector_free(eval1);
    gsl_vector_free(eval2);






    /*for(i = 0; i < dim; i++)
      for(j = 0; j < dim; j++) {
        double bla = gsl_matrix_get(g_m1, i, j);
        cout << bla;
      }
      cout << endl;
    }*/



    
    
    





  }  catch(const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}


int readin(string filename, stringstream &input){

	// Define stream (Convert into C-String)
	ifstream fin (filename.c_str());

	// Check if there is an error
	if (!fin.is_open()){
		cerr << "Could not open file!" << endl;
		return 1;
	}

	// Create a buffer to read in
	std::string s;

	// Read in the whole file
	while(true){
		getline(fin, s);
		// End reading in, if EOF is reached
		if(fin.eof())
			break;

		// Get rid of the comments starting with #
		s = s.substr(0,s.find('#'));

		// If something is left, add s to stream
		if (s!="")
			input << s << endl;

	}

	return 0;
}
