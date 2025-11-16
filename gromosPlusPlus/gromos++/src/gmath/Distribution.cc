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

// gmath_Distribution.cc
#include "Distribution.h"

#include <cmath>
#include <ostream>
#include <vector>
#include <iomanip>

using gmath::Distribution;

using namespace std;

namespace gmath{
  
Distribution::Distribution(double begin, double end, int nsteps):
  d_count(nsteps)
{
  
  if(begin>=end) 
    throw Distribution::Exception("Upper boundary should be higher than lower");
  if(nsteps<1) 
    throw Distribution::Exception("You need at least one step");
  d_step=(end-begin)/(nsteps);
  
  for(int i=0;i<nsteps;i++){
    d_count[i]=0;
    
  }
  d_nsteps=nsteps;
  d_begin=begin;
  d_end=end;
  d_sum=0.0;
  d_num=0;
}

  Distribution::Distribution(Distribution const & d):
    d_count(d.d_nsteps)
  {
    if(d.d_begin>=d.d_end) 
      throw Distribution::Exception("Upper boundary should be higher than lower");
    if(d.d_nsteps<1) 
      throw Distribution::Exception("You need at least one step");
    d_step=(d.d_end-d.d_begin)/(d.d_nsteps);
  
    for(int i=0;i<d.d_nsteps;i++){
      d_count[i]=0;
      
    }
    d_nsteps=d.d_nsteps;
    d_begin=d.d_begin;
    d_end=d.d_end;
    d_sum=0.0;
    d_num=0;
  }

   
void Distribution::write(std::ostream &os)const
{
  
  for(int i=0;i<d_nsteps;i++)
    os << setw(8) << d_begin+(i+0.5)*d_step << "\t" 
       << setw(5) << d_count[i] << endl;
}

void Distribution::write_normalized(std::ostream &os)const
{
  int nval = nVal();
  if (nval == 0) nval = 1;
  
  for(int i=0;i<d_nsteps;i++)
    os << setw(8) << d_begin+(i+0.5)*d_step << "\t" 
       << setw(5) << double(d_count[i]) / (nval * d_step) << endl;
}
  
double Distribution::add(const double value)
{
  if(value>=d_begin&&value<d_end){
     
     unsigned int q=int((value-d_begin)/d_step);
    if(q<d_count.size()){
     this->d_count[q]++;
     this->d_sum+=value;
     this->d_num++;
     return value;
    }
  }
  
  return value+1;
}

int Distribution::getbin(const double value)
{
  if(value>=d_begin&&value<d_end){
     
    unsigned int q=int((value-d_begin)/d_step);
    return q;
  }
  return -1;
}

bool Distribution::inrange(const double value) {
  if(value>=d_begin&&value<d_end) return true;
  else return false;
}

double Distribution::rmsd()const
{
  double sumdiff=0;
  double avr=this->ave();
  for(int i=0;i<d_nsteps;i++){
    double diff=avr - (d_begin+(i+0.5)*d_step);
    sumdiff+=d_count[i]*diff*diff;
  }
  return sqrt(sumdiff/d_num);
}

double Distribution::maxValAt() const
{
  int x_max = 0;
  int y_max = d_count[0];
  for(int i = 1 ; i < d_nsteps; ++i) {
    if(d_count[i] > y_max) {
      x_max = i;
      y_max = d_count[i];
    }
  }
  return d_begin+(x_max+0.5)*d_step;
}
 
void Distribution::clear()
{
  for(int i=0;i<d_nsteps;i++){
    d_count[i]=0;
  }
  d_num=0;
  d_sum=0.0;
}

}

