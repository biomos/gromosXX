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

// fit_RotationalFit.cc
//includes explicit calls to gsl now

#include "FastRotationalFit.h"

#include <cassert>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <vector>
#include <sstream>

#include <gsl/gsl_cblas.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_vector_double.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>

#include "PositionUtils.h"
#include "../gcore/System.h"
#include "../gcore/Molecule.h"
#include "../gmath/Matrix.h"
#include "../gmath/Vec.h"
#include "../utils/AtomSpecifier.h"


using gmath::Matrix;
using gmath::Vec;
using namespace fit;
using namespace std;



double FastRotationalFit::rmsd(Matrix const &rot, 
			       vector<Vec> const &ref, 
			       vector<Vec> const &sys)const {
  double rmsd2 = 0;
  double temp;
  if (d_rmsd_spec.size()) {
    for (size_t i = 0; i < ref.size(); ++i) {
      for (size_t j = 0; j < 3; ++j) {
        if (d_rmsd_spec[i]) {
          temp = 0;
          for (int b = 0; b < 3; ++b)
            temp += rot(j, b) * sys[i][b];

          rmsd2 += (ref[i][j] - temp) * (ref[i][j] - temp);
        }
      }
    }

    rmsd2 /= d_rmsd_num_atoms;
  } else {
    for (size_t i = 0; i < ref.size(); ++i) {
      for (size_t j = 0; j < 3; ++j) {
        temp = 0;
        for (int b = 0; b < 3; ++b)
          temp += rot(j, b) * sys[i][b];

        rmsd2 += (ref[i][j] - temp) * (ref[i][j] - temp);
      }
    }

    rmsd2 /= ref.size();
  }

  return sqrt(rmsd2);
}

int FastRotationalFit::fit(utils::AtomSpecifier & ref_spec,
			   utils::AtomSpecifier & sys_spec,
			   gcore::System & sys)const
{
  // sanity checks
  if (ref_spec.size() != sys_spec.size()){
    throw Exception("fitting: same number of atoms from reference and fittee needed");
  }
  
  if ((d_fit_spec.size() > 0) && (d_fit_spec.size() != unsigned(ref_spec.size()))){
    throw Exception("fitting: atoms to be used in fitting don't match reference specifier");
  }
  
  std::vector<Vec> v_ref(ref_spec.size()), 
    v_sys(sys_spec.size());

  for(unsigned int i=0; i<ref_spec.size(); ++i){
    v_ref[i] = ref_spec.pos(i);
    v_sys[i] = sys_spec.pos(i);
  }
  
  Matrix rot(3,3,0.0);
  int error;
  if (d_kabsch_fit)
    error = kabsch_fit(rot, v_ref, v_sys);
  else
    error = fit(rot, v_ref, v_sys);
  
  if (error){
    std::ostringstream os;
    os << "fitting: could not calculate rotation matrix (error " << error << ")";
    throw Exception(os.str());
  }
  
  PositionUtils::rotate(&sys, rot);
  return 0;
}

int FastRotationalFit::fit(vector<Vec> const &ref, 
			   vector<Vec> &sys)const{

  Matrix r(3,3,0);
  int error;

  if (d_kabsch_fit)
    error = kabsch_fit(r,ref, sys);
  else
    error = fit(r,ref, sys);

  if(error)
    return error;
  size_t num = ref.size();
  
  for(size_t n=0; n < num; ++n){
    sys[n] = r * sys[n];
  }
  return 0;
}

  
int FastRotationalFit::fit(Matrix &rot,
			   vector<Vec> const &ref, 
			   vector<Vec> const &sys)const{
  
  size_t num = ref.size();
  
  Matrix U(3,3,0);
  if(d_fit_spec.size()){
    
    for(int i=0;i<3;++i)
      for(int j=0;j<3;++j)
	for(size_t n=0;n<num;++n)
	  if(d_fit_spec[n])
	    U(i,j)+= sys[n][i]* ref[n][j];
    U *= 1.0/d_fit_num_atoms;
  }
  else{
    for(int i=0;i<3;++i)
      for(int j=0;j<3;++j)
	for(size_t n=0;n<num;++n)
	  U(i,j)+= sys[n][i]* ref[n][j];
    U *= 1.0/num;
  }
  
  double det=U.fastdet3X3Matrix();
  
  int signU = ( det>0 ? 1 : -1);


  gsl_matrix * omega = gsl_matrix_alloc (6, 6);
  gsl_matrix_set_zero (omega);

  for(int i=0;i<3;++i){
    for(int j=0;j<3;++j){
     gsl_matrix_set (omega, i, j+3,U(i,j));
     gsl_matrix_set (omega, i+3, j,U(j,i));
    }
  }

  
  double *eigenvals = new double [6];
  
  gsl_vector *eval = gsl_vector_alloc (6);
  gsl_matrix *evec = gsl_matrix_alloc (6,6);

  gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (6);

  gsl_eigen_symmv (omega, eval, evec, w);

  gsl_eigen_symmv_free(w);

  gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_VAL_DESC);

   Matrix Omega(6,6,0);
   for (int i=0; i < 6; ++i){
    eigenvals[i] = gsl_vector_get(eval, i);
   for (int j=0; j < 6; ++j){
     Omega(i,j)=gsl_matrix_get(evec, i, j);
   }
   }

  gsl_matrix_free (omega);
  gsl_matrix_free (evec);
  gsl_vector_free (eval);

  //  if(det<0){
  //  delete[] eigenvals;
  //  return -1;
  //}
  if(det<0 && fabs(eigenvals[1]-eigenvals[2])<1.0e-6){
    std::cout << "eigenval[1] = " << eigenvals[1] << "\n"
	      << "eigenval[2] = " << eigenvals[2] << "\n"
	      << "difference  = " << fabs(eigenvals[1]-eigenvals[2]) << "\n"
	      << "determinant = " << det << "\n"
	      << std::endl;
    delete[] eigenvals;
    return -2;
  }

  // Extract vectors from Omega.  
  Omega *= sqrt(2.0);
  Vec k1(Omega(0,0), Omega(1,0), Omega(2,0));
  Vec k2(Omega(0,1), Omega(1,1), Omega(2,1));
  Vec k3(Omega(0,2), Omega(1,2), Omega(2,2));
  Vec h1(Omega(3,0), Omega(4,0), Omega(5,0));
  Vec h2(Omega(3,1), Omega(4,1), Omega(5,1));
  Vec h3(Omega(3,2), Omega(4,2), Omega(5,2));

  double spat = h1.dot(h2.cross(h3));
  
  // turn 3rd vectors
  if(spat<0){
    h3=-h3;
    k3=-k3;
  }

  rot = Matrix(h1,k1) + Matrix(h2,k2) + signU*Matrix(h3,k3);

  delete[] eigenvals;

  return 0;
}


/* gsl does not provide it */
static inline void gsl_vector_cross(
  const gsl_vector *a,
  const gsl_vector *b,
  gsl_vector *c
) {
  double a0=gsl_vector_get(a,0);
  double a1=gsl_vector_get(a,1);
  double a2=gsl_vector_get(a,2);
  double b0=gsl_vector_get(b,0);
  double b1=gsl_vector_get(b,1);
  double b2=gsl_vector_get(b,2);
  gsl_vector_set(c,0,a1*b2-b1*a2);
  gsl_vector_set(c,1,a2*b0-b2*a0);
  gsl_vector_set(c,2,a0*b1-b0*a1);
}

#define NORM_EPS 0.00000001

int FastRotationalFit::kabsch_fit(Matrix &rot,
				  vector<Vec> const &ref, 
				  vector<Vec> const &sys)const{

  const unsigned int size = ref.size();
  if (size != sys.size())
    return -1;
  
  unsigned int i,j,k;
  // double n = 1.0 / size;

  int U_ok=0;

  // gsl_vector *cx=gsl_vector_alloc(3);     /* centroid of X */
  // gsl_vector *cy=gsl_vector_alloc(3);     /* centroid of Y */
  gsl_matrix *U=gsl_matrix_alloc(3,3);    /* rotation matrix */
  gsl_matrix *R=gsl_matrix_alloc(3,3);    /* Kabsch's R */
  gsl_matrix *RTR=gsl_matrix_alloc(3,3);  /* R_trans * R (and Kabsch's bk) */
  gsl_eigen_symmv_workspace *espace=gsl_eigen_symmv_alloc(3);
  gsl_matrix *evec=gsl_matrix_alloc(3,3); /* eigenvectors (and Kabsch's ak) */
  gsl_vector *eval=gsl_vector_alloc(3);   /* vector of eigenvalues */

  // should already be at cog
  /*
  // compute centroid of X
  gsl_vector_set_zero(cx);
  for(i=size;i>0;) {
    gsl_vector_const_view row=gsl_matrix_const_row(X,--i);
    gsl_vector_add(cx,&row.vector);
  } 
  gsl_vector_scale(cx,n);

  // compute centroid of Y
  gsl_vector_set_zero(cy);
  for(i=size;i>0;) {
    gsl_vector_const_view row=gsl_matrix_const_row(Y,--i);
    gsl_vector_add(cy,&row.vector);
  } 
  gsl_vector_scale(cy,n);

  // move X to origin
  for(i=size;i>0;) {
    gsl_vector_view row=gsl_matrix_row(X,--i);
    gsl_vector_sub(&row.vector,cx);
  }
  // move Y to origin
  for(i=size;i>0;) {
    gsl_vector_view row=gsl_matrix_row(Y,--i);
    gsl_vector_sub(&row.vector,cy);
  }
  */

  if(size==1) {
    /* just one point, so U is trival */
    gsl_matrix_set_identity(U);
  }
  else {
    /* compute R */
    gsl_matrix_set_zero(R);
    for(k=size;k>0;) {
      --k;
      for(i=3;i>0;) {
        --i;
        for(j=3;j>0;) {
          --j;
          gsl_matrix_set(R,i,j,
			 gsl_matrix_get(R,i,j)+
			 // gsl_matrix_get(Y,k,i)*gsl_matrix_get(X,k,j)
			 ref[k][i] * sys[k][j]
			 );
        }
      }
    }

    /* compute RTR = R_trans * R */
    gsl_matrix_set_zero(RTR);
    gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,R,R,0.0,RTR);

    /* compute orthonormal eigenvectors */
    gsl_eigen_symmv(RTR,eval,evec,espace);  /* RTR will be modified! */
    gsl_eigen_symmv_sort(eval,evec,GSL_EIGEN_SORT_VAL_DESC);
    if(gsl_vector_get(eval,1)>NORM_EPS) {
      /* compute ak's (as columns of evec) and bk's (as columns of RTR) */
      double norm_b0,norm_b1,norm_b2;
      gsl_vector_const_view a0=gsl_matrix_const_column(evec,0);
      gsl_vector_const_view a1=gsl_matrix_const_column(evec,1);
      gsl_vector_view a2=gsl_matrix_column(evec,2);
      gsl_vector_view b0=gsl_matrix_column(RTR,0);
      gsl_vector_view b1=gsl_matrix_column(RTR,1);
      gsl_vector_view b2=gsl_matrix_column(RTR,2);
      gsl_vector_cross(&a0.vector,&a1.vector,&a2.vector); /* a2 = a0 x a1 */
      gsl_blas_dgemv(CblasNoTrans,1.0,R,&a0.vector,0.0,&b0.vector);
      norm_b0=gsl_blas_dnrm2(&b0.vector);
      gsl_blas_dgemv(CblasNoTrans,1.0,R,&a1.vector,0.0,&b1.vector);
      norm_b1=gsl_blas_dnrm2(&b1.vector);
      if(norm_b0>NORM_EPS&&norm_b1>NORM_EPS) {
        gsl_vector_scale(&b0.vector,1.0/norm_b0);         /* b0 = ||R * a0|| */
        gsl_vector_scale(&b1.vector,1.0/norm_b1);         /* b1 = ||R * a1|| */
        gsl_vector_cross(&b0.vector,&b1.vector,&b2.vector);  /* b2 = b0 x b1 */

        norm_b2=gsl_blas_dnrm2(&b2.vector);
        if(norm_b2>NORM_EPS) {
          /* we reach this point only if all bk different from 0 */
          /* compute U = B * A_trans (use RTR as B and evec as A) */
          gsl_matrix_set_zero(U); /* to avoid nan */
          gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.0,RTR,evec,0.0,U);
        }
        else {
          U_ok=1;
          gsl_matrix_set_identity(U);
        }
      }
      else {
        U_ok=1;
        gsl_matrix_set_identity(U);
      }
    }
    else {
      U_ok=1;
      gsl_matrix_set_identity(U);
    }
  }

  /*
  // cx and cy are zero...
  // compute t = cy - U * cx
  gsl_vector_memcpy(t,cy);
  gsl_blas_dgemv(CblasNoTrans,-1.0,U,cx,1.0,t);
  */

  /*
  // no scaling factor...
  if(s) {
    // let us compute the optimal scaling as well
    // s = <Y,UX> / <UX,UX>
    *s=1.0;
    if(U_ok&&size>1) {
      double dom=0.0;
      double nom=0.0;
      double dom_i,nom_i;
      gsl_vector *Uxi=gsl_vector_alloc(3);
      for(i=size;i>0;) {
        gsl_vector_const_view row_x=gsl_matrix_const_row(X,--i);
        gsl_vector_const_view row_y=gsl_matrix_const_row(Y,i);
        gsl_vector_set_zero(Uxi);
        gsl_blas_dgemv(CblasNoTrans,1.0,U,&row_x.vector,1.0,Uxi);
        gsl_blas_ddot(&row_y.vector,Uxi,&nom_i);
        nom+=nom_i;
        gsl_blas_ddot(Uxi,Uxi,&dom_i);
        dom+=dom_i;
      }
      *s=nom/dom;
      gsl_vector_free(Uxi);
    }
  }
  */

  gsl_vector_free(eval);
  gsl_matrix_free(evec);
  gsl_eigen_symmv_free(espace);
  gsl_matrix_free(RTR);
  gsl_matrix_free(R);
  // gsl_vector_free(cy);
  // gsl_vector_free(cx);

  if (U_ok == 0){
    for(int i=0; i<3; ++i)
      for(int j=0; j<3; ++j)
	rot(i,j) = gsl_matrix_get(U, i, j);
  }
  else rot = Matrix(3,3,0.0);

  gsl_matrix_free(U);

  return U_ok;

}
