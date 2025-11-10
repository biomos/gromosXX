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

// gmath_Matrix.cc
#include "Matrix.h"

#include <cassert>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_vector_double.h>

#include "Vec.h"

using namespace std;

namespace gmath{
  Matrix::Matrix(int rows, int columns, double value){
    d_rows=rows;
    d_columns=columns;
    d_val=new double* [d_rows];
    for (int i=0;i<d_rows;++i){
      d_val[i]=new double [d_columns];
      for(int j=0;j<d_columns;++j){
	d_val[i][j]=value;
      }
    }
  }

  Matrix::Matrix(const Matrix &mat){
    d_rows=mat.d_rows;
    d_columns = mat.d_columns;
    d_val=new double* [d_rows];
    for (int i=0;i<d_rows;++i){
      d_val[i]=new double [d_columns];
      for(int j=0;j<d_columns;++j){
	d_val[i][j]=mat.d_val[i][j];
      }
    }
  }

  Matrix::Matrix(const Vec &v, const Vec &w){
    d_rows=3;
    d_columns=3;
    d_val=new double* [d_rows];
    for (int i=0;i<d_rows;++i){
      d_val[i]=new double [d_columns];
      for(int j=0;j<d_columns;++j){
	d_val[i][j]=v[i]*w[j];
      }
    }
  }

  Matrix::Matrix(const Vec &u, const Vec &v, const Vec &w){
    d_rows=3;
    d_columns=3;
    d_val=new double* [d_rows];
    for (int i=0;i<d_rows;++i){
      d_val[i]=new double [d_columns];
      d_val[i][0]=u[i];
      d_val[i][1]=v[i];
      d_val[i][2]=w[i];
    }
  }

  Matrix Matrix::transpose()const{
    Matrix m(d_columns, d_rows);
    for(int i=0; i<d_rows; ++i)
      for(int j=0; j<d_columns; ++j)
	m(j,i) = (*this)(i,j);
    return m;
  }

  Matrix &Matrix::operator=(const Matrix &mat){
    if (this != &mat){
      this->~Matrix();
      new(this) Matrix(mat);
    }
    return *this;
  }

  Matrix::~Matrix(){
    for(int i=0;i<d_rows;++i)
      delete[] d_val[i];
    delete[] d_val;
  }

  Matrix Matrix::luDecomp()const{
    assert(d_rows==d_columns);
    Matrix mat(*this);

    gsl_matrix * gsl_mat = gsl_matrix_alloc (mat.rows(), mat.columns());
    gsl_matrix_set_zero (gsl_mat);
    for (int i=0; i < mat.rows(); ++i){
      for (int j=0; j < mat.columns(); j++){
	gsl_matrix_set (gsl_mat, i , j , mat(i,j));
      }
    }  

    int s;
    gsl_permutation * p = gsl_permutation_alloc (mat.rows());
    gsl_linalg_LU_decomp (gsl_mat, p, &s);

    Matrix ret(mat.rows(), mat.columns());
    for (int i=0; i < mat.rows(); ++i){
      for (int j=0; j < mat.columns(); j++){
	ret(i,j)=gsl_matrix_get(gsl_mat, i, j);
      }
    }  

    return (ret);
  }

  Matrix Matrix::invert()const{
    assert(d_rows==d_columns);
    Matrix mat(*this);

    gsl_matrix * gsl_mat = gsl_matrix_alloc (mat.rows(), mat.columns());
    gsl_matrix_set_zero (gsl_mat);
    for (int i=0; i < mat.rows(); ++i){
      for (int j=0; j < mat.columns(); j++){
	gsl_matrix_set (gsl_mat, i , j , mat(i,j));
      }
    }

    int s;
    gsl_permutation * p = gsl_permutation_alloc (mat.rows());
    gsl_linalg_LU_decomp (gsl_mat, p, &s);

    gsl_matrix * inverse = gsl_matrix_alloc (mat.rows(), mat.columns());
    gsl_linalg_LU_invert(gsl_mat, p, inverse);


    Matrix ret(mat.rows(), mat.columns());
    for (int i=0; i < ret.rows(); ++i){
      for (int j=0; j < ret.columns(); j++){
	ret(i,j)=gsl_matrix_get(inverse, i, j);
      }
    }

    return (ret);
  }


  double Matrix::det()const{
    assert(d_rows==d_columns);
    Matrix tmp(*this);
 
    gsl_matrix * gsl_mat = gsl_matrix_alloc (tmp.rows(), tmp.columns());
    gsl_matrix_set_zero (gsl_mat);
    for (int i=0; i < tmp.rows(); ++i){
      for (int j=0; j < tmp.columns(); j++){
	gsl_matrix_set (gsl_mat, i , j , tmp(i,j));
      }
    }  
    gsl_permutation * p = gsl_permutation_alloc (tmp.rows());
    int s; 
    gsl_linalg_LU_decomp (gsl_mat, p, &s);


    double d = gsl_linalg_LU_det(gsl_mat, s);
    return d;
  }

  double Matrix::fastdet3X3Matrix()const{
    assert(d_rows==d_columns);
    Matrix tmp(*this); 
    //The usual elementary text book method:
    //D =   a11*a22*a33 + a21*a32*a13 + a31*a12*a23 
    //    - a11*a23*a32 - a21*a12*a33 - a31*a22*a13
    //is 12 multiplication, 2 addition and 3 subtractions.
    //this is an easy improvement:
    //D =  a11*(a22*a33-a23*a32)
    //   + a21*(a32*a13-a12*a33)
    //   + a31*(a12*a23-a22*a13)
    //is 9 multiplications, 3 subtractions and 2 additions.

    return  (tmp(0,0) * (tmp(1,1) * tmp(2,2) - tmp(1,2) * tmp(2,1))
	     + tmp(1,0) * (tmp(2,1) * tmp(0,2) - tmp(0,1) * tmp(2,2))
	     + tmp(2,0) * (tmp(0,1) * tmp(1,2) - tmp(1,1) * tmp(0,2)));

  }

  double fastdet3X3Matrix(gsl_matrix *gsl_mat) {

    return (gsl_matrix_get(gsl_mat, 0, 0) * (gsl_matrix_get(gsl_mat, 1, 1) * gsl_matrix_get(gsl_mat, 2, 2) 
					     - gsl_matrix_get(gsl_mat, 1, 2) * gsl_matrix_get(gsl_mat, 2, 1))
	    + gsl_matrix_get(gsl_mat, 1, 0) * (gsl_matrix_get(gsl_mat, 2, 1) * gsl_matrix_get(gsl_mat, 0, 2)
					       - gsl_matrix_get(gsl_mat, 0, 1) * gsl_matrix_get(gsl_mat, 2, 2))
	    + gsl_matrix_get(gsl_mat, 2, 0) * (gsl_matrix_get(gsl_mat, 0, 1) * gsl_matrix_get(gsl_mat, 1, 2) 
					       - gsl_matrix_get(gsl_mat, 1, 1) * gsl_matrix_get(gsl_mat, 0, 2)));
  }


  Matrix Matrix::diagonaliseSymmetric(double *eigenValues, bool sort){
    assert(d_rows==d_columns);
    Matrix mat(*this);
  
    gsl_matrix * gsl_mat = gsl_matrix_alloc (mat.rows(), mat.columns());
    gsl_matrix_set_zero (gsl_mat);
    for (int i=0; i < mat.rows(); ++i){
      for (int j=0; j < mat.columns(); j++){
	gsl_matrix_set (gsl_mat, i , j , mat(i,j));
      }
    }  

    gsl_vector *eval = gsl_vector_alloc (mat.rows());
    gsl_matrix *evec = gsl_matrix_alloc (mat.rows(), mat.columns());

    gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (mat.rows());

    gsl_eigen_symmv (gsl_mat, eval, evec, w);

    gsl_eigen_symmv_free(w);

    if (sort) gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_VAL_DESC);

    Matrix ret(mat.rows(), mat.columns());
    for (int i=0; i < mat.rows(); ++i){
      eigenValues[i] = gsl_vector_get(eval, i);
      for (int j=0; j < mat.columns(); j++){
	ret(i,j)=gsl_matrix_get(evec, i, j);
      }
    }  

    return ret;
  }

  Vec operator*(const Matrix &m, const Vec &v){
    assert(m.rows()==3&&m.columns()==3);
    Vec temp;
    for (int i=0;i<3;++i)
      for(int j=0;j<3;++j)
	temp[i]+=m(i,j)*v[j];
    return temp;
  }
  
  Matrix operator*(const Matrix &m1, const Matrix &m2){
    assert(m1.rows()==m2.columns()&&m1.columns()==m2.rows());
    Matrix temp(m1.rows(),m1.columns(),0);
    for(int i=0;i<m1.rows();++i)
      for(int j=0;j<m1.columns();++j)
	for(int k=0;k<m1.columns();++k)
	  temp(i,j)+=m1(i,k)*m2(k,j);
    return temp;
  }

}
