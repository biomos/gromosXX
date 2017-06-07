/**
 * @file colvar.cc
 * methods used by Colvars
 */
 
#include <limits>
#include "../../stdheader.h"
#include "../../interaction/special/colvar/colvar.h"

 
inline
double interaction::fastpow(double base, int exp)
{
    if(exp<0){
      exp=-exp;
      base=1.0/base;
    }
    double result = 1.0;
    while (exp)
    {
        if (exp & 1)
            result *= base;
        exp >>= 1;
        base *= base;
    }

    return result;
}

// copied from Plumed-2.2.0 switching function "rational"
double interaction::switchingfunction(double rdist,double&dfunc,int nn,int mm) {
      // Very small non-zero number
      const double epsilon(std::numeric_limits<double>::epsilon());
      
      double result;
      if(2*nn==mm){
// if 2*N==M, then (1.0-rdist^N)/(1.0-rdist^M) = 1.0/(1.0+rdist^N)
        double rNdist=fastpow(rdist,nn-1);
        double iden=1.0/(1+rNdist*rdist);
        dfunc = -nn*rNdist*iden*iden;
        result = iden;
      } else {
        if(rdist>(1.-100.0*epsilon) && rdist<(1+100.0*epsilon)){
           result=nn/mm;
           dfunc=0.5*nn*(nn-mm)/mm;
        }else{
           double rNdist=fastpow(rdist,nn-1);
           double rMdist=fastpow(rdist,mm-1);
           double num = 1.-rNdist*rdist;
           double iden = 1./(1.-rMdist*rdist);
           double func = num*iden;
           result = func;
           dfunc = ((-nn*rNdist*iden)+(func*(iden*mm)*rMdist));
        }
      }
    return result;
}

