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

#ifndef INCLUDED_ARGS_OUTFORMATPARSER
#define INCLUDED_ARGS_OUTFORMATPARSER

#include <string>

namespace gio {
class OutCoordinates;
}

namespace args {
class Arguments;
/**
 * @class OutformatParser
 * @ingroup args
 * @author N. Schmid
 * @date 04.11.2010
 *
 * Used to parse the coordinate output format argument (usually named
 * \@outformat). The following formats are supported:
 *
 * <table>
 * <tr><th>Argument</th><th>Description</th></tr>
 * <tr><td>cnf</td><td>Configuration format containing the POSITION
 * block.</td></tr> <tr><td>trc</td><td>Trajectory format containing the
 * POSITIONRED block.</td></tr> <tr><td>por</td><td>Position restraints
 * specification format.</td></tr> <tr><td>pdb [&lt;factor to convert length
 * unit to Angstrom, 10.0&gt; and/or &lt;"renumber" keyword to start numbering
 * at 1 at each molecule&gt;] </td><td>Protein Data Bank (PDB) format.</td></tr>
 * <tr><td>pqr [&lt;factor to convert length unit to Angstrom, 10.0&gt; and/or
 * &lt;"renumber" keyword to start numbering at 1 at each molecule&gt;]
 * </td><td>Modified Protein Data Bank (PDB) format.</td></tr> <tr><td>vmdam
 * [&lt;factor to convert length unit to Angstrom, 10.0&gt;]</td><td> VMD's
 * Amber Coordinates format.</td></tr>
 * </table>
 */
class OutformatParser {
private:
  OutformatParser(const OutformatParser &op);
  OutformatParser &operator=(const OutformatParser &op);

public:
  /**
   * parse the arguments and return the OutCoordinates for the correct format
   * @param args the arguments
   * @param ext the extension as a string (e.g. '.cnf')
   * @param argname the argument name used
   * @return pointer to the OutCoordinates
   */
  static gio::OutCoordinates *parse(Arguments &args, std::string &ext,
                                    std::string argname = "outformat");
};
} // namespace args

#endif /* INCLUDED_ARGS_OUTFORMATPARSER */
