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

// gio_OutCif.h

#ifndef INCLUDED_GIO_OUTCIF
#define INCLUDED_GIO_OUTCIF

#include "OutCoordinates.h"

#include <string>

using namespace std;

namespace gcore {
class System;
}

namespace utils {
class AtomSpecifier;
}

namespace gio {
class OutCif_i;
/**
 * Class OutCif
 * is of type OutCoordinates and defines how a single coordinate file
 * is printed out in the mmCIF format
 *
 * @class OutCif
 * @author @ref ega
 * @date 07.11.2023
 * @ingroup gio
 */

class OutCif : public OutCoordinates {
  OutCif_i *d_this;
  string data_name;
  double factor;
  bool renumber;

  // prevent copying and assignment
  OutCif(const OutCif &);
  OutCif &operator=(const OutCif &);

public:
  OutCif(double factor = 10.0, bool renumber = false);
  OutCif(ostream &os, double factor = 10.0, bool renumber = false);
  ~OutCif();
  void select(const string &thing);
  void open(ostream &os);
  void close();
  void writeTitle(const string &title);
  void writeTimestep(const int step, const double time);
  OutCif &operator<<(const gcore::System &sys);
  OutCif &operator<<(const utils::AtomSpecifier &atoms);
};
} // namespace gio

#endif
