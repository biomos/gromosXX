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
// gmath_StatDisk.tpp
#ifndef INCLUDED_GMATH_STATDISK_CC
#define INCLUDED_GMATH_STATDISK_CC

#include "StatDisk.h"

#include <cmath>
#include <cstdio>
#include <cassert>
#include <fstream>
#include <ios>
#include <stdio.h>
#include <string>
#include <vector>

#include "../gromos/Exception.h"
#include "Vec.h"

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
  StatDisk<T>::StatDisk(std::string file, size_t autoflush)
    : d_counter(0),
      d_ave(T()),
      d_msd(T()),
      d_ee(T()),
      d_avedone(false),
      d_msddone(false),
      d_eedone(false),
      d_file(file),
      d_autoflush(autoflush)
  { 
    open(file);
  }
  
    template<typename T>
  StatDisk<T>::StatDisk(size_t autoflush)
    : d_counter(0),
      d_ave(T()),
      d_msd(T()),
      d_ee(T()),
      d_avedone(false),
      d_msddone(false),
      d_eedone(false),
      d_file(""),
      d_autoflush(autoflush)
  { 
  }

  template<typename T>
  StatDisk<T>::~StatDisk()
  {
    d_out.close();
    remove(d_file.c_str());
  }

  template<typename T>
  void StatDisk<T>::addval(T val)
  {
    d_vals.push_back(val);
    d_counter++;
    d_avedone=false;
    d_msddone=false;
    d_eedone=false;
    
    if (d_vals.size() == d_autoflush)
      flush();
  }
  
  template<typename T>
  T StatDisk<T>::ave()
  {
    if(!d_avedone){
      d_ave = this->subave(0,d_counter);
      d_avedone=1;
    }
    return d_ave;
  }
  
  template<typename T>
  T StatDisk<T>::subave(int b, int e)
  {
    //calculate the average
    T ave = 0;
    flush();
    
    assert(b >= 0 && b <= e);
    assert(e <= d_counter);
    
    std::ifstream in(d_file.c_str(), std::ios::binary);
    if (!in.is_open())
      throw gromos::Exception("StatDisk", "Cannot open file.");
    
    for(int i=0; i<b; ++i) {
      T tmp;
      in.read(reinterpret_cast<char*>(&tmp), sizeof(T));
    }
    
    for(int i=b;i<e;i++){
      T tmp;
      in.read(reinterpret_cast<char*>(&tmp), sizeof(T));
      ave += tmp;
    }
    in.close();
    return ave/(e-b);
  }
  
  template<typename T>
  T StatDisk<T>::msd()
  {
    if(!d_msddone) {
      flush();
      std::ifstream in(d_file.c_str(), std::ios::binary);
      if (!in.is_open())
        throw gromos::Exception("StatDisk", "Cannot open file.");
      T sum=0, ssum=0;
      for(int i=0; i<d_counter; i++){
        T tmp;
        in.read(reinterpret_cast<char*>(&tmp), sizeof(T));
	sum+=tmp;
	ssum+=tmp*tmp;
      }
      in.close();
      
      sum/=d_counter;
      ssum/=d_counter;
      d_msd = ssum - sum*sum;
    }
    return d_msd;
  }

  template<typename T>
  T StatDisk<T>::rmsd()
  {
    return sqrt(this->msd());
  }

  template<typename T>
  T StatDisk<T>::min()
  {
    flush();
    std::ifstream in(d_file.c_str(), std::ios::binary);
    if (!in.is_open())
      throw gromos::Exception("StatDisk", "Cannot open file.");
    T m;
    in.read(reinterpret_cast<char*>(&m), sizeof(T));
    for(int i=1; i<d_counter; ++i) {
      T tmp;
      in.read(reinterpret_cast<char*>(&tmp), sizeof(T));
      if(tmp < m) m=tmp;
    }
    in.close();
    return m;
  }
  
  template<typename T>
  T StatDisk<T>::max()
  {
    flush();
    std::ifstream in(d_file.c_str(), std::ios::binary);
    if (!in.is_open())
      throw gromos::Exception("StatDisk", "Cannot open file.");
    T m;
    in.read(reinterpret_cast<char*>(&m), sizeof(T));
    for(int i=1; i<d_counter; ++i) {
      T tmp;
      in.read(reinterpret_cast<char*>(&tmp), sizeof(T));
      if(tmp > m) m=tmp;
    }
    in.close();
    return m;
  }
  
  template<typename T>
  T StatDisk<T>::ee()
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
  void StatDisk<T>::flush()
  {
    size_t num = d_vals.size();
    for(size_t i = 0; i < num; ++i)
      d_out.write(reinterpret_cast<char *>(&d_vals[i]), sizeof(T));
    d_vals.clear();
    d_out.flush();  
  }
  
  template<typename T>
  void StatDisk<T>::open(std::string file)
  {
    d_file = file;
    d_out.open(d_file.c_str(), std::ios::binary);
    if (!d_out.is_open())
      throw gromos::Exception("StatDisk", "Cannot write scratch file: " + d_file);
  }
}

#endif
