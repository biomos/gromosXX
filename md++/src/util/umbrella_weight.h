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
 * @file umbrella_weight.h
 * fancy weighting methods for umbrellas
 */

#ifndef INCLUDED_UMBRELLA_WEIGHT_H
#define	INCLUDED_UMBRELLA_WEIGHT_H

namespace util {
  /**
   * @class Umbrella_Weight
   * represents the weight of a certain configuration
   */
  class Umbrella_Weight {
  public:
    virtual ~Umbrella_Weight() {}
    /**
     * get the weight of this configuration
     */
    virtual double get_weight() const = 0;
    operator double() { return get_weight(); }
    /**
     * increment the weight of this configuration
     */
    virtual void increment_weight() = 0;
    void operator++() { increment_weight(); }
    /**
     * write the weight to a stream
     */
    virtual void write(std::ostream & os) const = 0;
    /**
     * read the weight from a stream
     */
    virtual void read(std::istream & is) = 0;
  };

  /**
   * @class Umbrella_Weight_Factory
   * create instances of umbrella weights
   */
  class Umbrella_Weight_Factory {
  public:
    /**
     * get an instance
     */
    virtual Umbrella_Weight * get_instance() = 0;
    /**
     * clean up
     */
    virtual ~Umbrella_Weight_Factory();
  protected:
    /**
     * save the instances in here to make sure they are deleted
     */
    std::vector<Umbrella_Weight*> instances;
  };

  /**
   * @class Number_Of_Visits_Umbrella_Weight
   * A simple number of visits umbrella weight
   */
  class Number_Of_Visits_Umbrella_Weight : public Umbrella_Weight {
  public:
    virtual ~Number_Of_Visits_Umbrella_Weight() {}
    Number_Of_Visits_Umbrella_Weight() { weight = 1; }
    virtual double get_weight() const { return double(weight); }
    virtual void increment_weight() { ++weight; }
    virtual void write(std::ostream & os) const;
    virtual void read(std::istream & is) { is >> weight; }
  protected:
    unsigned int weight;
  };

  /**
   * @class Number_Of_Visits_Umbrella_Weight_Factory
   * The factory for the number of visits umbrella weight
   */
  class Number_Of_Visits_Umbrella_Weight_Factory : public Umbrella_Weight_Factory {
  public:
    Number_Of_Visits_Umbrella_Weight_Factory() {}
    virtual Umbrella_Weight * get_instance() {
      Umbrella_Weight * inst = new Number_Of_Visits_Umbrella_Weight();
      instances.push_back(inst);
      return inst;
    }
  };
}

std::ostream & operator<<(std::ostream & os, const util::Umbrella_Weight & w);
std::istream & operator>>(std::istream & is, util::Umbrella_Weight & w);

#endif

