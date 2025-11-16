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

// gmath_Stat.tpp
#ifndef INCLUDED_GMATH_STAT_CC
#define INCLUDED_GMATH_STAT_CC

#include "Stat.h"

#include <cmath>
#include <iostream>
#include <ostream>
#include <vector>

#include "Distribution.h"
#include "Vec.h"
#include "../gromos/Exception.h"


namespace gmath
{
  inline double sqrt(double d)
  {
    return ::sqrt(d);
  }
  
  inline gmath::Vec sqrt(gmath::Vec const & v)
  {
    return gmath::Vec(::sqrt(v[0]),::sqrt(v[1]), ::sqrt(v[2]));
  }

  inline gmath::Vec operator/(gmath::Vec const &v1, gmath::Vec const &v2)
  {
    return gmath::Vec(v1[0]/v2[0], v1[1]/v2[1], v1[2]/v2[2]);
  }
  
  template<typename T>
  Stat<T>::Stat()
    : d_counter(0),
      d_ave(T()),
      d_msd(T()),
      d_ee(T()),
      d_avedone(false),
      d_lnexpavedone(false),
      d_msddone(false),
      d_eedone(false),
      d_distdone(false),
      d_dist(0, 1, 1)
  { 
  }
  
  template<typename T>
  void Stat<T>::addval(T val)
  {
    d_vals.push_back(val);
    d_counter++;
    d_avedone=false;
    d_lnexpavedone=false;
    d_msddone=false;
    d_eedone = false;
    if (d_distdone) d_dist.add(val);
  }
  
  template<typename T>
  T Stat<T>::ave()const
  {
    if(!d_avedone){
      d_ave = this->subave(0,d_counter);
      d_avedone=1;
    }
    return d_ave;
  }
  
  template<typename T>
  T Stat<T>::subave(int b, int e)const
  {
    T ave = 0;
    
    //calculate the average
    for(int i=b;i<e;i++){
      ave += d_vals[i];
    }
    return ave/(e-b);
  }

  template<typename T>
  T Stat<T>::lnexpave()const
  {
    if(!d_lnexpavedone){
      d_lnexpave = this->lnexpsubave(0,d_counter);
      d_lnexpavedone=1;
    }
    return d_lnexpave;
  }

  template<typename T>
  T Stat<T>::lnexpsubave(int b, int e)const
  {
    T  lnexpave = d_vals[b];

    //calculate the average
    for(int i=b+1;i<e;i++) {
      //log<exp(d_vals[i])>
      // for numerical reasons follow B.A. Berg, Comput. Phys. Comm. 153 (2003) 397
      lnexpave = std::max(lnexpave, d_vals[i])
              + log(1.0 + exp(std::min(lnexpave, d_vals[i]) - std::max(lnexpave, d_vals[i])));
    }
    return lnexpave - log(double(e - b));
  }

  template<typename T>
  T Stat<T>::lnXexpave(Stat<T> X, Stat<T> Y, int &sign) {
    if (X.d_counter != Y.d_counter)
      throw gromos::Exception("Stat", "Can't calculate the ln|<Xexp[Y]>|. Unequal number of elements.");
    // calculate the average ln|<Xexp(Y)>|
    // for numerical reasons follow B.A. Berg, Comput. Phys. Comm. 153 (2003) 397
    // set ave to starting value; subtract it again in the end
    double eps = Y.d_vals[0];
    T ave = eps;
    int sign_ave = 1;
    // loop over all data
    for (int i = 0; i < X.d_counter; i++) {
      // check for the sign of X
      int signX = 1;
      if (X.d_vals[i] < 0) signX = -1;
      // |X|*exp[Y]
      T term = log(fabs(X.d_vals[i])) + Y.d_vals[i];
      // case I: both signs are positive: lnC = ln (A+B)
      // and case II: X and ave are negative -A -B = - (A + B)
      //if( (signX > 0 && sign_ave > 0) || (signX < 0 && sign_ave < 0) ){
      if (signX * sign_ave > 0) {
        // sum up; case I: sign stays positive; case II: sign stays negative
        ave = std::max(ave, term) + log(1.0 + exp(std::min(ave, term) - std::max(ave, term)));
      } // case III: X is negative and ave is positive ln|C| = ln|ave  -term|
        // and case IV: X positive and ave negative ln|C| = ln |term -ave|
        // else if( (signX < 0 && sign_ave > 0) || (signX > 0 && sign_ave < 0) ){
      else if (signX * sign_ave < 0) {
        // determine sign of ave (case IV, so far negative)
        if (term > ave) sign_ave = +1 * signX;
        // for case III the sign is the opposite, that is why we multiply by signX
        // which is +1 for case IV (no change) and -1 for case III (sign change)
        // sum up                         !
        ave = std::max(ave, term) + log(1.0 - exp(std::min(ave, term) - std::max(ave, term)));
      } else
        throw gromos::Exception("Stat", "Error when evaluating ln|<Xexp[Y]>|.");
    } // loop over all data
    // subtract eps again, determine sign

    if (sign_ave < 0) {
      // -ave -eps = - (ave + eps); sign stays untouched
      ave = std::max(ave, eps) + log(1.0 + exp(std::min(ave, eps) - std::max(ave, eps)));
    } else {
      // +ave -eps; determine sign; ave so far positive
      if (ave < eps) sign_ave = -1;
      ave = std::max(ave, eps) + log(1.0 - exp(std::min(ave, eps) - std::max(ave, eps)));
    }
    sign = sign_ave;

    return ave - log(double(X.d_counter));

  }

  template<typename T>
  T Stat<T>::covariance(Stat<T> X, Stat<T> Y) {
    if (X.d_counter != Y.d_counter)
      throw gromos::Exception("Stat", "Can't calculate the covariance. Unequal number of elements.");
    T sumX = 0, sumY = 0, ssum = 0;
    for(int i = 0; i < X.d_counter; i++){
      sumX += X.d_vals[i];
      sumY += Y.d_vals[i];
      ssum += X.d_vals[i] * Y.d_vals[i];
    }
    sumX /= X.d_counter;
    sumY /= Y.d_counter;
    ssum /= X.d_counter;
    
    T cov = ssum - sumX * sumY;
    return cov;
  }

  template<typename T>
  T Stat<T>::lnexpcovariance(Stat<T> X, Stat<T> Y, int & sign) {
    //ln{Cov(exp(X),exp(Y))}=ln{<exp(X)exp(Y)>-<exp(X)><exp(Y)>}
    // for numerical reasons follow B.A. Berg, Comput. Phys. Comm. 153 (2003) 397
    if (X.d_counter != Y.d_counter)
      throw gromos::Exception("Stat", "Can't calculate the covariance. Unequal number of elements.");
    T cov = X.d_vals[0] + Y.d_vals[0];
    //calculate the average ln<exp(X)exp(Y)>
    for (int i = 1; i < X.d_counter; i++) {
      T d_vals_prod = X.d_vals[i] + Y.d_vals[i];
      cov = std::max(cov, d_vals_prod)
              + log(1.0 + exp(std::min(cov, d_vals_prod) - std::max(cov, d_vals_prod)));
    }
    // "divide" by the number of data points
    cov = cov - log(double(X.d_counter));

    // calculate ln{<exp(X)exp(Y)>-<exp(X)><exp(Y)>}
    // store ln{<exp(X)><exp(Y)>}
    T ave_prod = X.lnexpave() + Y.lnexpave();
    sign = +1;
    if (cov < ave_prod)
      sign = -1;
    cov = std::max(cov, ave_prod) + log(1.0 - exp(std::min(cov, ave_prod) - std::max(cov, ave_prod)));
    
    return cov;
  }

  template<typename T>
  T Stat<T>::msd()const
  {
    if(!d_msddone){
      T sum=0, ssum=0;
      for(int i=0; i<d_counter; i++){
	sum+=d_vals[i];
	ssum+=d_vals[i]*d_vals[i];
      }
      sum/=d_counter;
      ssum/=d_counter;
      d_msd = ssum - sum*sum;
    }
    return d_msd;
  }

  template<typename T>
  T Stat<T>::rmsd()const
  {
    return sqrt(this->msd());
  }

  template<typename T>
  T Stat<T>::min()const
  {
    T m=d_vals[0];
    for(int i=1; i<d_counter; ++i)
      if(d_vals[i] < m) m = d_vals[i];
    return m;
  }
  
  template<typename T>
  T Stat<T>::max()const
  {
    T m = d_vals[0];
    for(int i=1; i<d_counter; ++i)
      if(d_vals[i] > m) m=d_vals[i];
    return m;
  }
  
  template<typename T>
  T Stat<T>::ee()const
  {
    if(!d_eedone){
      // first prepare the blocks
      double blksz=50;
      int old=2;
      while(4*blksz<d_counter){
	d_blocksize.push_back(int(blksz));
	old=int(blksz);
	while(old==int(blksz)) blksz = blksz*1.07177;
      }
      
      int Nblocks=d_blocksize.size();
      T rmsd2, ave=0;
      T runave=this->ave();
      T runrmsd=this->rmsd();
      std::vector<T> fit(Nblocks), x(Nblocks);
      
      for(int j=0; j<Nblocks; j++){
	int Nblcki=d_counter/d_blocksize[j];
	
	// The rmsd of the property we are interested in, weighted by the
	// average energy of the blocks
	rmsd2=0;
	for(int i=0; i<Nblcki; i++){
	  ave = this->subave(i*d_blocksize[j],(i+1)*d_blocksize[j]);
	  rmsd2+=(ave-runave)*(ave-runave);
	}
	rmsd2/=Nblcki;
	fit[j]=(d_blocksize[j]*rmsd2) / runrmsd / runrmsd;
	x[j]=1.0/d_blocksize[j];
	
      }
      T sx=0, sf=0,sfx=0,sxx=0;
      for(int i=0; i<Nblocks;i++){
	sx+=x[i];
	sf+=fit[i];
	sfx+=x[i]*fit[i];
	sxx+=x[i]*x[i];
      }
      
      T a, b;
      a=(sf*sx/Nblocks-sfx)/(sx*sx/Nblocks-sxx);
      b = (sf - a*sx)/Nblocks;

      d_ee=sqrt(b/d_counter)*runrmsd;
    }
    
    return d_ee;
  }

  template<typename T>
  T Stat<T>::stat_ineff(Stat<T> X, Stat<T> Y) {
    if (X.d_counter != Y.d_counter)
      throw gromos::Exception("Stat", "Can't calculate the statistical inefficiency. Unequal number of elements.");
    // Get the length of the timeseries
    int N = X.d_counter;
    // Initialize statistical inefficiency estimate with uncorrelated value
    T statisticalInefficiency = 1.0;
    // Compute means and variance
    T mu_X = X.ave();
    T mu_Y = Y.ave();
    T sigma2_XY = gmath::Stat<T>::covariance(X, Y);
    if (sigma2_XY == 0) {
      std::cerr << "Variance zero! Unable to compute the statistical inefficiency" << std::endl;
    }
    /*
       Accumulate the integrated correlation time by computing the normalized correlation time at
       increasing values of t.  Stop accumulating if the correlation function goes negative, since
       this is unlikely to occur unless the correlation function has decayed to the point where it
       is dominated by noise and indistinguishable from zero.
     */
    int t = 1;
    int increment = 1;
    do {
      //  Compute unnormalized correlation function for time t.
      // C = sum( A_n(1:(N-t))*B_n((1+t):N) + B_n(1:(N-t))*A_n((1+t):N) ) / (2.0 * (N-t))
      double C = 0.0;
      for (int n = 1; n <= (N - t); n++) {
        C += X.d_vals[n - 1] * Y.d_vals[n + t - 1] + Y.d_vals[n - 1] * X.d_vals[n + t - 1];
      }
      C = C / double(2.0 * (N - t));
      // Compute normalized fluctuation correlation functions from unnormalized correlation functions.
      C = (C - mu_X * mu_Y) / sigma2_XY;

      // Terminate if the correlation function has crossed zero.
      if (C <= 0) break;

      // Accumulate contribution to the statistical inefficiency.
      statisticalInefficiency = statisticalInefficiency + 2.0 * C * (1.0 - double(t) / double(N)) * increment;
      // Increment t and the amount by which we increment t.
      t = t + increment;
      increment = increment + 1;
    } while (t < N - 1);
    return statisticalInefficiency;

  }

  template<typename T>
  T Stat<T>::lnexp_stat_ineff(Stat<T> X, Stat<T> Y) {
    if (X.d_counter != Y.d_counter)
      throw gromos::Exception("Stat", "Can't calculate the statistical inefficiency. Unequal number of elements.");
    // Get the length of the timeseries
    int N = X.d_counter;
    // Initialize statistical inefficiency estimate with uncorrelated value
    T statisticalInefficiency = 0.0;
    // Compute means and variance
    T mu_X = X.lnexpave();
    T mu_Y = Y.lnexpave();
    int sign = 1;
    T sigma2_XY = gmath::Stat<T>::lnexpcovariance(X, Y, sign);
    if (exp(sigma2_XY) == 0) {
      std::cerr << "Variance zero! Unable to compute the statistical inefficiency" << std::endl;
    }
    /*
       Accumulate the integrated correlation time by computing the normalized correlation time at
       increasing values of t.  Stop accumulating if the correlation function goes negative, since
       this is unlikely to occur unless the correlation function has decayed to the point where it
       is dominated by noise and indistinguishable from zero.
     */
    int t = 1;
    int increment = 1;
    do {
      //  Compute unnormalized correlation function for time t.
      // C = sum( A_n(1:(N-t))*B_n((1+t):N) + B_n(1:(N-t))*A_n((1+t):N) ) / (2.0 * (N-t))
      // double C = 0.0;
      T prod1 = X.d_vals[0] + Y.d_vals[t];
      T prod2 = Y.d_vals[0] + X.d_vals[t];
      double C = std::max(prod1,prod2) + log(1.0 + exp(std::min(prod1,prod2) - std::max(prod1,prod2)));
      for (int n = 2; n <= (N - t); n++) {
        //C += X.d_vals[n - 1] * Y.d_vals[n + t - 1] + Y.d_vals[n - 1] * X.d_vals[n + t - 1];
        // compute x_n*y_{n+t} and y_n*x_{n+t}
        prod1 = X.d_vals[n - 1] + Y.d_vals[n + t - 1];
        prod2 = Y.d_vals[n - 1] + X.d_vals[n + t - 1];
        // now add x_n*y_{n+t} to y_n*x_{n+t}
        T term = std::max(prod1,prod2) + log(1.0 + exp(std::min(prod1,prod2) - std::max(prod1,prod2)));
        // now add everything to C
        C = std::max(C,term) + log(1.0 + exp(std::min(C,term) - std::max(C,term)));
      }
      //C = C / double(2.0 * (N - t));
      C -= log(double(2.0 * (N - t)));
      
      // Compute normalized fluctuation correlation functions from unnormalized correlation functions.
      //C = (C - mu_X * mu_Y) / sigma2_XY;
      // 1. comput mu_X * mu_Y
      T aveprod = mu_X + mu_Y;
      // 2. subtract; sign is determined by sigma2_XY and by "C < aveprod"
      if (C < aveprod) sign *= -1;
      C = std::max(C,aveprod) + log(1.0 - exp(std::min(C,aveprod) - std::max(C,aveprod)));
      // 3. divide by sigma2_XY
      C -= sigma2_XY;

      // Terminate if the correlation function has crossed zero.
      if (sign <= 0) break;

      // determine the term to add to the si
      double term = log(2.0) + C + log(1.0 - double(t) / double(N)) + log(double(increment));
      // Accumulate contribution to the statistical inefficiency.
      // because of the break statement above we can be sure that C is always positive
      if (sign > 0){
        // add the term
        // Accumulate contribution to the statistical inefficiency.
        statisticalInefficiency = std::max(statisticalInefficiency,term)
                + log(1.0 + exp(std::min(statisticalInefficiency,term) - std::max(statisticalInefficiency,term)));
      }
      else
        throw gromos::Exception("Stat","time correlation function negative. Forgotten break statement?");
      
      // Increment t and the amount by which we increment t.
      t = t + increment;
      increment = increment + 1;
    } while (t < N - 1);
    return statisticalInefficiency;
  }

  template<typename T>
  std::vector<T> const & Stat<T>::data() const{
    return d_vals;
  }

  template<typename T>
  gmath::Distribution const & Stat<T>::distribution()const
  {
    if(d_distdone)
      return d_dist;
    else
      throw Distribution::Exception("call dist_init first");
  }
  
  template<typename T>
  gmath::Distribution const & Stat<T>::dist_init(double lower, double upper, int nsteps, bool periodic)
  {
    d_dist = gmath::Distribution(lower, upper, nsteps);
    //put all values in it
    for(int i=0; i<d_counter; i++) {
      double d = d_vals[i];
      if(periodic) {
        double period = upper - lower;
        while(d >= upper) {
          d -= period;
        }
        while(d < lower) {
          d += period;
        }
      }
      d_dist.add(d);
    }
    d_distdone=1;
    
    //return the distribution
    return d_dist;
  }

  template<typename T>
  void Stat<T>::subtract_average()
  {
    double ave=this->ave();
    for(int i=0; i<d_counter; i++)
      d_vals[i]-=ave;
    d_avedone=0;
  }
}

#endif
