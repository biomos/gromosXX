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
#include "Correlation.h"

#include <cstdlib>
#include <iostream>
#include <vector>

#include "../gromos/Exception.h"

using gmath::Correlation;

using namespace std;

int main() {
  try {

    vector<double> a;
    vector<double> b;
    srand(23123);

    for (int j = 0; j < 1000; j++) {
      a.push_back(rand() * 1000.0 / RAND_MAX);
      b.push_back(j);
    }
    Correlation c(a, b);
    Correlation d(a, b);
    c.calc_direct();
    //for(int i=0; i< c.size(); i++){
    //	cout << i << "\t" << c[i] << endl;
    //  }
    d.calc_fft();

    for (int i = 0; i < c.size(); i++) {
      cout << i << "\t" << c[i] << "\t" << d[i] << endl;
    }

    vector<double> w(1000);
    vector<double> s(1000);
    d.spectrum(w, s, 0.002, 0.5);
    for (size_t i = 0; i < w.size(); i++)
      cout << w[i] << "\t" << s[i] << endl;

    return 0;
  }  catch (gromos::Exception e) {
    cerr << e.what() << endl;
    return 1;
  }
}
