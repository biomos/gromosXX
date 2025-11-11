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

// gcore_LJException.h

#ifndef INCLUDED_GCORE_LJEXCEPTION
#define INCLUDED_GCORE_LJEXCEPTION

#include <set>

namespace gcore{

  /**
   * Class LJException
   * Purpose: contains a LJ exception
   *
   * Description:
   * Contains the atoms and type making up a LJException. The atoms are sorted
   * to have the lowest atom number first.
   *
   * @class LJException
   * @author A. Eichenberger
   * @ingroup gcore
   */
class LJException{
  int d_a[2];
  int d_type;
  std::set<int> d_cond;
  int d_ind;
  // not implemented
  LJException();
 public:
  /**
   * LJException constructor
   * constructs a new LJException defined by atoms a and b. These atoms are stored
   * such that a<b.
   * @param a,b atom numbers of the atoms making the LJException
   */
  LJException(int a, int b);
  /**
   * LJException copy constructor
   * Constructs a new LJException and copies the specied LJException into it
   * @param & LJException to be copied
   */
  LJException(const LJException &);
  /**
   * LJException deconstructor
   */
  ~LJException(){}
  /**
   * Member operator =, copies one LJException into the other
   */
  LJException &operator=(const LJException &b);
  /**
   * Method setType sets the LJException type of the LJException
   * @param i LJException type to be set
   */
  void setType(int i){d_type=i;}
  /**
   * Member operator [], returns the atom number of the i-th atom in the LJException
   * @param i atom index in the LJException (0 or 1)
   * @return Atom number of the i-th atom in the LJException
   */
  int operator[](int i)const{return d_a[i];}
  /**
   * Member operator [], returns the atom number of the i-th atom in the LJException
   * @param i atom index in the LJException (0 or 1)
   * @return Atom number of the i-th atom in the LJException
   */
  int &operator[](int i){return d_a[i];}
  /**
   * Accessor, returns the LJException type of the LJException
   */
  int type()const{return d_type;}
  /**
   * Accessor, returns the set of conditions
   */
  const std::set<int> cond() const { return d_cond; }
  /**
   * Accessor, returns the set of conditions
   */
  std::set<int> & cond() { return d_cond; }
  /**
   * Adds a condition to the LJ Exception
   */
  void addCond(int i) { d_cond.insert(i); }
  /**
   * Returns the number of LJ exception conditions
   */
  int numcond() const { return d_cond.size(); }
  /**
   * Returns the inicator numbor which tells you on what atoms the conditions
   * apply:
   * 0: on both atoms
   * 1: on the first atom
   * 2: on the second atom
   */
  int indicate()const {return d_ind;}
  int & indicate() {return d_ind;}
};
/**
 * @relates LJException
 * Operator < compares two LJExceptions a and b
 *
 * a is defined to be smaller than b if<br><ol>
 * <li> The atom number of the first atom in LJException a is lower than the
 *      atom number of the first atom in LJException b (a[0]<b[0]), or
 * <li> The atom numbers of the first atoms in both LJExceptions are the same
 *      but the second atom has a lower atom number in LJException a (a[0]==b[0]
 *      && a[1]<b[1]), or
 * <li> All atoms are the same, but the LJException type of LJException a has a lower
 *      number than the LJException type of LJException b (a[0]==b[0] && a[1]==b[1] &&
 *      a.type() < b.type()).
 * </ol>
 * &param a, b LJExceptions to be compared
 * &return 1 if a<b; 0 otherwise
 */
int operator<(const LJException &a, const LJException &b);
}
#endif
