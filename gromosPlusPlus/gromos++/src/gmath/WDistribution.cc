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

// gmath_WDistribution.cc
#include "WDistribution.h"

#include <cmath>
#include <ostream>
#include <vector>
#include <iomanip>

#include "Distribution.h"

using gmath::WDistribution;

using namespace std;

namespace gmath {

  WDistribution::WDistribution(double begin, double end, int nsteps) :
  Distribution(begin, end, nsteps), d_count_d_stat(nsteps), d_count_d(nsteps) {
    for (int i = 0; i < nsteps; i++) {
      d_count_d[i] = 0.0;
    }
    d_num_d = 0.0;
  }

  WDistribution::WDistribution(WDistribution const & d) :
  Distribution(d), d_count_d_stat(d.d_nsteps), d_count_d(d.d_nsteps) {
    for (int i = 0; i < d.d_nsteps; i++) {
      d_count_d[i] = 0;
    }
    d_num_d = 0;
  }

  double WDistribution::add(const double value, const double weight) {
    if (value >= d_begin && value < d_end) {

      unsigned int q = int((value - d_begin) / d_step);
      if (q < d_count.size()) {
        /* first store all the values in a Stat object. This allows to use later
         * the function lnexpave of this class for summing up the numbers. This
         * is done to avoid numerical problems when evaluating exponential functions.
         */
        this->d_count_d_stat[q].addval(weight);
        this->d_num_d_stat.addval(weight);
        this->d_num++;
        return value;
      }
    }
    return value + 1;
  }

  void WDistribution::write(std::ostream &os) {
    for (int i = 0; i < d_nsteps; i++) {
      if (this->d_count_d_stat[i].n() > 0)
        // sum up (i.e. average and then multiply by n()) values in bin 
        this->d_count_d[i] = exp(this->d_count_d_stat[i].lnexpave()) * this->d_count_d_stat[i].n();
    }

    for (int i = 0; i < d_nsteps; i++)
      os << setw(8) << d_begin + (i + 0.5) * d_step << "\t"
            << setw(5) << d_count_d[i] << endl;
  }

  void WDistribution::write_normalized(std::ostream &os) {

    for (int i = 0; i < d_nsteps; i++) {
      if (this->d_count_d_stat[i].n() > 0)
        /* sum up (i.e. average and then multiply by n()) values in bin and
         * divide by the total sum of weights. Division by the total sum of weights
         * is done in 2 steps: 1) here: divide by the average weight (in a Distribution
         * this is one); 2) divide by the number of data points (done below when
         * writing out).
         * The average of weights is subtracted (=^ divison) within the loop for
         * numerical reasons.
         */
        this->d_count_d[i] = exp(this->d_count_d_stat[i].lnexpave() - this->d_num_d_stat.lnexpave())
        * this->d_count_d_stat[i].n();
    }

    int nval = this->d_num_d_stat.n();
    if (nval == 0)
      nval = 1;
    // write out normalized
    for (int i = 0; i < d_nsteps; i++) {
      os << setw(8) << d_begin + (i + 0.5) * d_step << "\t"
              << setw(5) << d_count_d[i] / double(nval * d_step) << endl;
    }
  }


}

