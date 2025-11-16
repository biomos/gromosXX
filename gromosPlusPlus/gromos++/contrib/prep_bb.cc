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
 * @file prep_bb.cc
 * prepares a building block
 */

#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <map>

/**
 * @page contrib Contrib Program Documentation
 *
 * @anchor prep_bb
 * @section prep_bb prepares a building block
 * @author @ref co ega
 * @date 16. 3. 2005
 *
 * Program prep_bb takes information from a specification or PDB file and
creates
 * a building block. If required, it makes suggestions for the choice of
parameters.
 *
 * The specification file has to be given in a certain format:
 * @verbatim
TITLE
ETOH
END
ATOMS
# specification of the atoms
# name    element
  C1      C
  C2      C
  O2      O
  H2      H
 END
 BONDS
 # specification of bonds
 # i    j    type (1=single, 2=double, 3=triple, 4=delocalized/aromatic)
   1    2    1
   2    3    1
   3    4    1
 END
 @endverbatim
 * The file consists of three blocks: a TITLE block which holds the name of the
 * building block, an ATOMS block which specifies the names of the atoms and a
 * BONDS block which specifies the bonds. Atomic elements and bond types are
 * only required for graph based parameter suggestions (see below) and can be
omitted.
 *
 * If a PDB file is given bonds are automatically detected: all atoms pairs
 * with a interatomic distance lower than a given limit (\@bound) are considered
 * to be bonded.
 *
 * The order of the atoms in the resulting building block can be controled by
 * specifying the first atom using the \@reorder argument.
 *
 * If force-field files (and a graph library) are given, the
 * program interacts (\@interact) with the user and provides advice for the
choice of
 * appropriate parameters. Else, the program chooses invalid default parameters
 * and writes a building block skeletton for further manual processing.
 *
 * prep_bb is able the predict parameters using a graph based algorithm in a
 * accurate way. In order to use this features you have to provide a
specification
 * file (the information in a PDB is not sufficient) and give element and bond
 * type information (see above). In addition, a graph library is needed to
translate
 * the existing building blocks in the MTB files to a graph. The graph library
should
 * look like this:
  * @verbatim
TITLE
Graph library for force field
END
FORCEFIELD
53A6
END
ELEMENTMAPPING
# mass  element (capitals only. i.e. Cl should be CL)
1       H
3       C
4       C
5       C
...
END
BONDMAPPING
# bond type   bond order
# bond orders: 1.0 : single
#              1.5 : aromatic bond/delocalized
#              2.0 : double bond, 3.0: triple bond.
1   1.0
2   1.0
...
 END
 @endverbatim
 *
 * The resulting building block is written to BUILDING.out
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@spc</td><td>&lt;specification file&gt; OR</td></tr>
 * <tr><td> \@pdb</td><td>&lt;PDB file&gt; </td></tr>
 * <tr><td> \@bound</td><td>&lt;upper bound for bond length (for PDB)&gt;
</td></tr>
 * <tr><td> \@reorder</td><td>&lt;first atom&gt;</td></tr>
 * <tr><td>[\@build</td><td>&lt;building block file&gt;]</td></tr>
 * <tr><td>[\@param</td><td>&lt;parameter file&gt;]</td></tr>
 * <tr><td>[\@interact</td><td></td></tr>
 * <tr><td>[\@graph_library</td><td>&lt;library to translate BB into
graphs.&gt;]</td></tr>
 * </table>
 *
 * Example:
 * @verbatim
   prep_bb
     @pdb           ligand.pdb
     @build         mtb53a6.dat
     @param         ifp53a6.dat
     @interact
     @graph_library graphlib.53a6
   @endverbatim

 * <hr>
 */

#include <cassert>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "../src/args/Arguments.h"
#include "../src/gcore/Angle.h"
#include "../src/gcore/AngleType.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/BbSolute.h"
#include "../src/gcore/Bond.h"
#include "../src/gcore/BondType.h"
#include "../src/gcore/BuildingBlock.h"
#include "../src/gcore/Dihedral.h"
#include "../src/gcore/DihedralType.h"
#include "../src/gcore/Exclusion.h"
#include "../src/gcore/GromosForceField.h"
#include "../src/gcore/Improper.h"
#include "../src/gcore/ImproperType.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gio/Ginstream.h"
#include "../src/gio/InBuildingBlock.h"
#include "../src/gio/InParameter.h"
#include "../src/gio/OutBuildingBlock.h"
#include "../src/gmath/Vec.h"
#include "../src/gromos/Exception.h"
#include "../src/utils/CommandLine.h"
#include "../src/utils/FfExpert.h"
#include "../src/utils/FfExpertGraph.h"

namespace cmd = utils::CommandLine;

static const char program_name[] = "prep_bb";
constexpr const char *BLUE = "\033[1;34m";
constexpr const char *YELLOW = "\033[1;33m";
constexpr const char *GREEN = "\033[1;32m";
constexpr const char *RED = "\033[1;31m";
constexpr const char *RESET = "\033[0m";

template <typename key_type, typename value_type>
std::map<value_type, key_type>
swap_key_value(const std::map<key_type, value_type> &m) {
  std::map<value_type, key_type> res;
  typename std::map<key_type, value_type>::const_iterator it = m.begin();
  typename std::map<key_type, value_type>::const_iterator to = m.end();
  for (; it != to; ++it) {
    res[it->second] = it->first;
  }
  return res;
}

std::string read_input(std::string file, std::vector<std::string> &atom,
                       std::vector<std::string> &element,
                       std::vector<gcore::Bond> &bond);
std::string read_pdb(std::string file, std::vector<std::string> &atom,
                     std::vector<std::string> &element,
                     std::vector<gcore::Bond> &bond, double bondbound);
/*
string read_mol(string file, vector<string> &atom, vector<string> & element,
vector<Bond> & bond);
*/
std::set<int> neighbours(int a, const std::vector<gcore::Bond> &bond);
void forward_neighbours(int a, std::vector<gcore::Bond> &bond,
                        std::set<int> &cum, int prev);
void add_neighbours(int i, std::vector<gcore::Bond> &bond, std::set<int> &at,
                    std::vector<int> &order);
int count_bound(int i, std::vector<gcore::Bond> &bond, std::set<int> &at,
                std::set<int> &j);
std::set<int> ring_atoms(std::vector<gcore::Bond> &bond, int numatoms);

std::string &operator++(std::string &s) { return s = s + "  "; }

std::string &operator--(std::string &s) {
  return s = s.substr(0, s.length() - 2);
}

const int bignumber = 1000;

class PrepBBCore {
private:
  gcore::BbSolute *bbs{nullptr};

public:
  explicit PrepBBCore(const std::string &name) : bbs{new gcore::BbSolute} {
    bbs->setResName(name.substr(0, name.size() - 1));
  }
  ~PrepBBCore() { delete bbs; }

  void addAtom2BB(const std::vector<int> &order,
                  const std::vector<std::string> &atom_name, bool interact,
                  bool has_graph, const utils::FfExpert &exp,
                  utils::FfExpertGraph *graph,
                  const gcore::GromosForceField &gff,
                  const std::vector<gcore::Bond> &bond,
                  const std::map<int, int> &newnum, const std::set<int> &aro,
                  const std::set<int> &bt_aro);
  bool addBond2BB(bool interact, const utils::FfExpert &exp,
                  const gcore::GromosForceField &gff,
                  std::vector<gcore::Bond> &bondnew);
  bool addAngle2BB(bool interact, const utils::FfExpert &exp,
                   const gcore::GromosForceField &gff,
                   std::vector<gcore::Angle> &newangles);
  bool addImproper2BB(bool interact, const utils::FfExpert &exp,
                      const gcore::GromosForceField &gff,
                      std::vector<gcore::Improper> &newimpropers);
  bool addDihedral2BB(bool interact, const utils::FfExpert &exp,
                      const gcore::GromosForceField &gff,
                      std::vector<gcore::Dihedral> &newdihedrals);
  void writeBB();
};

void PrepBBCore::addAtom2BB(
    const std::vector<int> &order, const std::vector<std::string> &atom_name,
    bool interact, bool has_graph, const utils::FfExpert &exp,
    utils::FfExpertGraph *graph, const gcore::GromosForceField &gff,
    const std::vector<gcore::Bond> &bond, const std::map<int, int> &newnum,
    const std::set<int> &aro, const std::set<int> &bt_aro) {
  double totcharge = 0;
  int numchargegroup = 0;
  int i = 0;
  while (i < static_cast<int>(order.size())) {
    std::string aname = atom_name[order[i]].substr(0, 1);
    int mass = 0;
    int iac = 0;
    double charge = 0.0;
    int chargegroup = 0;

    if (interact) {
      std::cerr << BLUE
                << "\n--------------------------------------------------------"
                << "----------------------\n";
      std::cerr << "Considering atom " << i + 1 << " ("
                << atom_name.at(order[i]) << ")\n"
                << RESET;

      if (has_graph) {
        std::vector<std::vector<utils::Vertex>> hits;

        exp.substructure2iac(order[i], *graph, hits);
        std::vector<std::map<int, double>> stat(hits.size());

        bool found = false;
        for (size_t radius = hits.size(); radius-- > 0;) {
          if (hits[radius].empty()) {
            found = true;
          }
          for (const auto &vertex : hits[radius]) {
            if (stat[radius].count(vertex.iac)) {
              stat[radius][vertex.iac] += 100.0 / hits[radius].size();
            } else {
              stat[radius][vertex.iac] = 100.0 / hits[radius].size();
            }
          }
        }

        if (found) {
          std::cerr << "\n"
                    << "Suggested Integer Atom Codes for atom (" << aname
                    << ") "
                    << "based on substructure matching: \n";
          bool first = true;
          for (size_t radius = hits.size(); radius-- > 0;) {
            if (hits[radius].empty()) {
              continue;
            }

            std::cerr << std::setw(8) << radius << ": ";
            std::cerr.precision(2);
            std::cerr.setf(std::ios::fixed, std::ios::floatfield);

            if (first) {
              std::cerr << RED;
            }

            std::map<double, int> sorted_stat = swap_key_value(stat[radius]);

            std::map<double, int>::reverse_iterator it = sorted_stat.rbegin(),
                                                    to = sorted_stat.rend();
            for (; it != to; ++it) {
              std::cerr << (it->second + 1) << " : " << std::setw(5)
                        << it->first << " %  ";
            }

            if (first) {
              first = false;
              std::cerr << RESET;
            }

            for (size_t hit = 0; hit < hits[radius].size(); ++hit) {
              if (hit % 6 == 0) {
                std::cerr << std::endl << std::setw(10) << " ";
              }

              const auto &h = hits[radius][hit];
              std::cerr << h.iac + 1 << " " << h.residue << "[" << h.name
                        << "]";
              if (hit != hits[radius].size() - 1) {
                std::cerr << ", ";
              }
            }
            std::cerr << std::endl;
          }
        } else {
          std::cerr << "\n"
                    << "No Integer Atom Code found based on topological "
                       "similarity!\n";
        }
      }
      std::cerr.precision(4);

      std::vector<utils::FfExpert::counter> ocList;
      // get IAC
      exp.name2iac(aname, ocList);
      if (!ocList.empty()) {
        size_t ap = utils::sort(ocList, true);

        std::cerr << "\n"
                  << "Suggested Integer Atom Code based on the first letter "
                  << "of the name (" << aname << ") :\n";

        for (size_t iList = 0; iList < ocList.size(); ++iList) {
          if (iList == ap) {
            std::cerr << RED;
          }
          std::cerr << std::setw(10) << ocList[iList].type + 1 << " :"
                    << std::setw(5) << gff.atomTypeName(ocList[iList].type)
                    << " occurs " << std::setw(3) << ocList[iList].occurence
                    << " times\n";

          if (iList == ap) {
            std::cerr << RESET;
          }
        }
        std::cerr << GREEN << std::setw(10) << "undo"
                  << " : In case you want to go back to the previous atom\n"
                  << RESET;
      }

      std::ostringstream prompt;
      prompt << "Give IAC ( 1 - " << gff.numAtomTypeNames() << " ): ";

      std::string line = cmd::getLine(prompt.str());
      if (line == "undo") {
        if (i > 0) {
          std::cerr << YELLOW << "\nUndo: going back to previous atom\n\n"
                    << RESET;
          --i;
        } else {
          std::cerr << YELLOW
                    << "\nAlready at the first atom, cannot go back\n\n"
                    << RESET;
        }
        continue;
      }

      std::istringstream s(line);
      s >> iac;
      --iac;

      if (iac < 0 || iac >= gff.numAtomTypeNames()) {
        std::cerr << RED << "This IAC (" << iac + 1 << ") is not "
                  << "defined in the given parameter file" << RESET;
        continue;
      }

      if (has_graph) {
        // save the iac in the vertex
        graph->vertices()[order[i]].iac = iac;
      }

      // get mass
      exp.iac2mass(iac, ocList);
      if (!ocList.empty()) {
        size_t ap = utils::sort(ocList, true);

        std::cerr << "\n"
                  << "Suggested Mass Code based on the IAC "
                  << "of the atom (" << iac + 1 << ") :\n";

        for (size_t iList = 0; iList < ocList.size(); iList++) {
          if (iList == ap) {
            std::cerr << RED;
          }
          std::cerr << std::setw(10) << ocList[iList].type + 1 << " :"
                    << std::setw(8) << gff.findMass(ocList[iList].type)
                    << " occurs " << std::setw(3) << ocList[iList].occurence
                    << " times\n";

          if (iList == ap) {
            std::cerr << RESET;
          }
        }
        std::cerr << GREEN << std::setw(10) << "undo"
                  << " : In case you want to go back to the previous atom\n"
                  << RESET;
      }

      line = cmd::getLine("Give Masstype: ");
      if (line == "undo") {
        if (i > 0) {
          std::cerr << YELLOW << "\nUndo: going back to previous atom\n\n"
                    << RESET;
          --i;
        } else {
          std::cerr << YELLOW
                    << "\nAlready at the first atom, cannot go back\n\n"
                    << RESET;
        }
        continue;
      }

      s.clear();
      s.str(line);
      s >> mass;
      --mass;

      if (gff.findMass(mass) == 0.0) {
        std::cerr << RED << "This Masstype (" << mass + 1 << ") is "
                  << "not defined in the parameter file\n";
        continue;
      }

      // charge
      if (has_graph) {
        std::vector<std::vector<utils::Vertex>> hits;
        exp.substructure2iac(order[i], *graph, hits);
        std::vector<std::map<double, double>> stat(hits.size());

        bool found = false;
        for (size_t radius = hits.size(); radius-- > 0;) {
          if (hits[radius].empty()) {
            found = true;
          }

          for (size_t hit = 0; hit < hits[radius].size(); ++hit) {
            const auto &h = hits[radius][hit];
            if (stat[radius].find(h.charge) != stat[radius].end()) {
              stat[radius][h.charge] += 100.0 / hits[radius].size();
            } else {
              stat[radius][h.charge] += 100.0 / hits[radius].size();
            }
          }
        }

        if (found) {
          std::cerr << "\n"
                    << "Suggested charge for atom (" << aname << ") "
                    << "based on substructure matching: \n";
          bool first = true;
          for (size_t radius = hits.size(); radius-- > 0;) {
            if (hits[radius].empty()) {
              continue;
            }

            std::cerr << std::setw(8) << radius << ": ";
            std::cerr.setf(std::ios::fixed, std::ios::floatfield);

            if (first) {
              std::cerr << RED;
            }

            std::map<double, double> sorted_stat = swap_key_value(stat[radius]);
            std::map<double, double>::reverse_iterator it =
                                                           sorted_stat.rbegin(),
                                                       to = sorted_stat.rend();
            for (; it != to; ++it) {
              std::cerr.precision(4);
              std::cerr << it->second;
              std::cerr.precision(2);
              std::cerr << " : " << std::setw(5) << it->first << " %  ";
            }
            if (first) {
              first = false;
              std::cerr << RESET;
            }

            std::cerr.precision(3);
            for (size_t hit = 0; hit < hits[radius].size(); ++hit) {
              const auto &h = hits[radius][hit];
              if (hit % 6 == 0) {
                std::cerr << std::endl << std::setw(10) << " ";
              }

              std::cerr << std::setw(5) << h.charge << " " << h.residue << "["
                        << h.name << "]";
              if (hit != hits[radius].size() - 1) {
                std::cerr << ", ";
              }
            }
            std::cerr << std::endl;
          }
        } else {
          std::cerr << "\n"
                    << "No charge found based on topological similarity!\n";
        }
      }

      std::cerr.precision(4);

      // get charge
      exp.iac2charge(iac, ocList);
      if (!ocList.empty()) {
        size_t ap = utils::sort(ocList, false);

        std::cerr << "\n"
                  << "Suggested Charge based on the IAC "
                  << "of the atom (" << iac + 1 << ") :\n";

        for (size_t iList = 0; iList < ocList.size(); ++iList) {
          if (iList == ap) {
            std::cerr << RED;
          }
          std::cerr << std::setw(10) << exp.charge(ocList[iList].type)
                    << " occurs " << std::setw(3) << ocList[iList].occurence
                    << " times\n";

          if (iList == ap) {
            std::cerr << RESET;
          }
        }
        std::cerr << GREEN << std::setw(10) << "undo"
                  << " : In case you want to go back to the previous atom\n"
                  << RESET;

        if (numchargegroup) {
          std::cerr << "Charge group so far (" << numchargegroup << " atom";
          if (numchargegroup > 1) {
            std::cerr << "s";
          }
          std::cerr << ") has net charge of " << totcharge << "\n\n";
        }
      }

      line = cmd::getLine("Give Charge: ");
      if (line == "undo") {
        if (i > 0) {
          std::cerr << YELLOW << "\nUndo: going back to previous atom\n\n"
                    << RESET;
          --i;
        } else {
          std::cerr << YELLOW
                    << "\nAlready at the first atom, cannot go back\n\n"
                    << RESET;
        }
        continue;
      }

      s.clear();
      s.str(line);
      s >> charge;
      totcharge += charge;
      ++numchargegroup;

      std::cerr << "\n"
                << "Suggested Chargegroup code based on the total charge "
                << "of this\n"
                << "charge group (" << totcharge << "; " << numchargegroup
                << " atom";
      if (numchargegroup > 1) {
        std::cerr << "s";
      }
      std::cerr << "): ";
      if ((static_cast<int>(std::rint(totcharge * 1000)) % 1000) == 0) {
        std::cerr << "1\n\n";
      } else {
        std::cerr << "0\n\n";
      }

      std::set<int> zeroone;
      zeroone.insert(0);
      zeroone.insert(1);
      cmd::getValueFromSet<int>(chargegroup, zeroone,
                                "Give Chargegroup code (0/1): ");

      if (chargegroup == 1) {
        numchargegroup = 0;
        totcharge = 0.0;
      }
    }

    // done gathering data. Prepare atom
    gcore::AtomTopology a_top;
    a_top.setName(atom_name[order[i]]);
    a_top.setMass(mass);
    a_top.setIac(iac);
    a_top.setCharge(charge);
    a_top.setChargeGroup(chargegroup);
    bbs->setResNum(i, 0);
    // It is a building block -> 14 exclusions are not strictly necessary
    // but in the case of aromaticity we might need it
    std::set<int> exclusions;
    std::set<int> exclusions14;
    std::set<int> nb1 = neighbours(order[i], bond);
    std::set<int> nb2;
    std::set<int> nb3;
    std::set<int> tmp;

    for (int inb : nb1) {
      exclusions.insert(newnum.at(inb));
      nb2 = neighbours(inb, bond);
      for (int jnb : nb2) {
        exclusions.insert(newnum.at(jnb));
        tmp = neighbours(jnb, bond);
        for (int knb : tmp) {
          nb3.insert(knb);
        }
      }
    }
    for (int inb : nb3) {
      if (!exclusions.count(newnum.at(inb))) {
        exclusions14.insert(newnum.at(inb));
      }
    }

    gcore::Exclusion ex, ex14;
    for (int inb : exclusions) {
      if (inb > static_cast<int>(i)) {
        ex.insert(inb);
      }
    }

    for (int inb : exclusions14) {
      bool i_is_aromatic = aro.count(i) || bt_aro.count(i);
      bool inb_is_aromatic = aro.count(inb) || bt_aro.count(inb);

      if (i_is_aromatic && inb_is_aromatic) {
        if (inb > static_cast<int>(i)) {
          ex.insert(inb);
        }
      } else if (inb > static_cast<int>(i)) {
        ex14.insert(inb);
      }
    }

    a_top.setExclusion(ex);
    a_top.setExclusion14(ex14);
    if (interact) {
      std::cerr << "\n\tDetermined exclusions and 1,4-exclusions.\n";
    }

    bbs->addAtom(a_top);
    ++i;
  }
}

bool PrepBBCore::addBond2BB(bool interact, const utils::FfExpert &exp,
                            const gcore::GromosForceField &gff,
                            std::vector<gcore::Bond> &bondnew) {
  if (interact) {
    std::cerr << BLUE
              << "======= ATOMIC INFORMATION GATHERED ===================="
              << "======================" << RESET << '\n';

    std::cerr << BLUE << "======= IS WRITTEN TO FILE "
              << "============================="
              << "======================\n\n\n";

    std::cerr << "======= BOND PARAMETERS ================================"
              << "======================" << RESET << "\n\n";
  }

  std::vector<utils::FfExpert::counter> ocList;
  int i = 0;
  while (i < static_cast<int>(bondnew.size())) {
    if (interact) {
      gcore::Bond iacbond(bbs->atom(bondnew[i][0]).iac(),
                          bbs->atom(bondnew[i][1]).iac());
      std::cerr << BLUE
                << "\n--------------------------------------------------------"
                << "----------------------\n";
      std::cerr << "Considering bond " << bondnew[i][0] + 1 << " ("
                << bbs->atom(bondnew[i][0]).name() << ")  -  "
                << bondnew[i][1] + 1 << " (" << bbs->atom(bondnew[i][1]).name()
                << ")\n"
                << RESET;

      exp.iac2bond(iacbond, ocList);
      if (!ocList.empty()) {
        size_t ap = utils::sort(ocList);
        std::cerr << "\n"
                  << "Suggested Bondtype based on the IAC "
                  << "of the atoms (" << iacbond[0] + 1 << ", "
                  << iacbond[1] + 1 << ") :\n";
        for (size_t iList = 0; iList < ocList.size(); ++iList) {
          if (ap == iList) {
            std::cerr << RED;
          }
          std::cerr << std::setw(10) << ocList[iList].type + 1
                    << " : Kb = " << std::setw(10) << std::setprecision(1)
                    << gff.bondType(ocList[iList].type).fc()
                    << "; b0 = " << std::setw(7) << std::setprecision(3)
                    << gff.bondType(ocList[iList].type).b0() << " occurs "
                    << std::setw(3) << ocList[iList].occurence << " times\n";
          if (ap == iList) {
            std::cerr << RESET;
          }
        }
        std::cerr << GREEN << std::setw(10) << "undo"
                  << " : In case you want to go back to the previous bond\n"
                  << RESET;
        std::cerr << YELLOW << std::setw(10) << "goback"
                  << " : In case you want to go back to Atom Parameters\n"
                  << RESET;
      }

      // Get user input
      std::ostringstream prompt;
      prompt << "Give Bondtype ( 1 - " << gff.numBondTypes() << " ) : ";

      std::string line = cmd::getLine(prompt.str());
      if (line == "goback") {
        std::cerr << YELLOW << "\nUndo: going back to Atom Parameters\n\n"
                  << RESET;
        return true;
      } else if (line == "undo") {
        if (i > 0) {
          std::cerr << YELLOW << "\nUndo: going back to previous bond\n\n"
                    << RESET;
          --i;
        } else {
          std::cerr << YELLOW
                    << "\nAlready at the first bond, cannot go back\n\n"
                    << RESET;
        }
        continue;
      }

      std::istringstream s(line);
      int bt;
      s >> bt;
      if (s.fail()) {
        std::cerr << RED << "Invalid input. Try again.\n" << RESET;
        continue;
      }

      if (bt < 1 || bt > gff.numBondTypes()) {
        std::cerr << RED << "This Bondtype (" << bt << ") is not "
                  << "defined in the given parameter file" << RESET;
        continue;
      }

      bondnew[i].setType(bt - 1);
    }
    bbs->addBond(bondnew[i]);
    ++i;
  }
  return false;
}

bool PrepBBCore::addAngle2BB(bool interact, const utils::FfExpert &exp,
                             const gcore::GromosForceField &gff,
                             std::vector<gcore::Angle> &newangles) {
  std::vector<utils::FfExpert::counter> ocList;
  if (interact && !newangles.empty()) {
    std::cerr << BLUE
              << "\n\n======= ANGLE PARAMETERS ==============================="
              << "======================\n"
              << RESET;
  }

  int i = 0;
  while (i < static_cast<int>(newangles.size())) {
    auto &angle = newangles[i];

    if (interact) {
      int i0 = angle[0];
      int i1 = angle[1];
      int i2 = angle[2];

      gcore::Angle iacangle(bbs->atom(i0).iac(), bbs->atom(i1).iac(),
                            bbs->atom(i2).iac());

      std::cerr << BLUE
                << "\n--------------------------------------------------------"
                << "----------------------\n";
      std::cerr << "Considering angle " << i0 + 1 << " ("
                << bbs->atom(i0).name() << ")  -  " << i1 + 1 << " ("
                << bbs->atom(i1).name() << ")  -  " << i2 + 1 << " ("
                << bbs->atom(i2).name() << ")\n"
                << RESET;

      exp.iac2angle(iacangle, ocList);
      if (!ocList.empty()) {
        size_t ap = utils::sort(ocList);

        std::cerr << "\nSuggested Angletype based on the IAC "
                  << "of the atoms (" << iacangle[0] + 1 << ", "
                  << iacangle[1] + 1 << ", " << iacangle[2] + 1 << ") :\n";

        for (size_t iList = 0; iList < ocList.size(); ++iList) {
          if (ap == iList) {
            std::cerr << RED;
          }
          std::cerr << std::setw(10) << ocList[iList].type + 1
                    << " : Kb = " << std::setw(10) << std::setprecision(1)
                    << gff.angleType(ocList[iList].type).fc()
                    << "; t0 = " << std::setw(8) << std::setprecision(3)
                    << gff.angleType(ocList[iList].type).t0() << " occurs "
                    << ocList[iList].occurence << " times\n";
          if (ap == iList) {
            std::cerr << RESET;
          }
        }
        std::cerr << GREEN << std::setw(10) << "undo"
                  << " : In case you want to go back to the previous angle\n"
                  << RESET;
        std::cerr << YELLOW << std::setw(10) << "goback"
                  << " : In case you want to go back to Bond Parameters\n"
                  << RESET;
      }

      std::ostringstream prompt;
      prompt << "Give Angletype ( 1 - " << gff.numAngleTypes() << " ) : ";

      std::string line = cmd::getLine(prompt.str());
      if (line == "goback") {
        std::cerr << YELLOW << "\nUndo: going back to Bond Parameters\n\n"
                  << RESET;
        return true;
      } else if (line == "undo") {
        if (i > 0) {
          std::cerr << YELLOW << "\nUndo: going back to previous angle\n\n"
                    << RESET;
          --i;
        } else {
          std::cerr << YELLOW
                    << "\nAlready at the first angle, cannot go back\n\n"
                    << RESET;
        }
        continue;
      }

      std::istringstream s(line);
      int bt;
      s >> bt;
      if (s.fail()) {
        std::cerr << RED << "Invalid input. Try again.\n" << RESET;
        continue;
      }

      if (bt < 1 || bt > gff.numAngleTypes()) {
        std::cerr << RED << "This Angletype (" << bt << ") is not "
                  << "defined in the given parameter file" << RESET;
        continue;
      }

      angle.setType(bt - 1);
    }
    bbs->addAngle(angle);
    ++i;
  }
  return false;
}

bool PrepBBCore::addImproper2BB(bool interact, const utils::FfExpert &exp,
                                const gcore::GromosForceField &gff,
                                std::vector<gcore::Improper> &newimpropers) {
  std::vector<utils::FfExpert::counter> ocList;
  if (interact && !newimpropers.empty()) {
    std::cerr << BLUE
              << "\n\n======= IMPROPER PARAMETERS ============================"
              << "======================\n"
              << RESET;
  }

  int i = 0;
  while (i < static_cast<int>(newimpropers.size())) {
    auto &improper = newimpropers[i];
    if (interact) {
      int i0 = improper[0];
      int i1 = improper[0];
      int i2 = improper[0];
      int i3 = improper[0];

      gcore::Improper iacimproper(bbs->atom(i0).iac(), bbs->atom(i1).iac(),
                                  bbs->atom(i2).iac(), bbs->atom(i3).iac());

      std::cerr << BLUE
                << "\n--------------------------------------------------------"
                << "----------------------\n";
      std::cerr << "Considering improper " << i0 + 1 << " ("
                << bbs->atom(i0).name() << ")  -  " << i1 + 1 << " ("
                << bbs->atom(i1).name() << ")  -  " << i2 + 1 << " ("
                << bbs->atom(i2).name() << ")  -  " << i3 + 1 << " ("
                << bbs->atom(i3).name() << ")\n"
                << RESET;

      exp.iac2improper(iacimproper, ocList);
      if (ocList.size()) {
        size_t ap = utils::sort(ocList);

        std::cerr << "\n"
                  << "Suggested Impropertype based on the IAC "
                  << "of the atoms (" << iacimproper[0] + 1 << ", "
                  << iacimproper[1] + 1 << ", " << iacimproper[2] + 1 << ", "
                  << iacimproper[3] + 1 << ") :\n";

        for (size_t iList = 0; iList < ocList.size(); ++iList) {
          if (iList == ap) {
            std::cerr << RED;
          }

          std::cerr << "\t\t" << std::setw(10) << ocList[iList].type + 1
                    << " : Kb = " << std::setw(10) << std::setprecision(1)
                    << gff.improperType(ocList[iList].type).fc()
                    << "; q0 = " << std::setw(10) << std::setprecision(3)
                    << gff.improperType(ocList[iList].type).q0() << " occurs "
                    << ocList[iList].occurence << " times\n";

          if (iList == ap) {
            std::cerr << RESET;
          }
        }
        std::cerr << GREEN << std::setw(10) << "undo"
                  << " : In case you want to go back to the previous improper\n"
                  << RESET;
        std::cerr << YELLOW << std::setw(10) << "goback"
                  << " : In case you want to go back to Angle Parameters\n"
                  << RESET;
      }

      std::ostringstream prompt;
      prompt << "Give Impropertype ( 1 - " << gff.numImproperTypes() << " ) : ";

      std::string line = cmd::getLine(prompt.str());
      if (line == "goback") {
        std::cerr << YELLOW << "\nUndo: going back to Angle Parameters\n\n"
                  << RESET;
        return true;
      } else if (line == "undo") {
        if (i > 0) {
          std::cerr << YELLOW << "\nUndo: going back to previous improper\n\n"
                    << RESET;
          --i;
        } else {
          std::cerr << YELLOW
                    << "\nAlready at the first improper, cannot go back\n\n"
                    << RESET;
        }
        continue;
      }

      std::istringstream s(line);
      int bt;
      s >> bt;
      if (s.fail()) {
        std::cerr << RED << "Invalid input. Try again.\n" << RESET;
        continue;
      }

      if (bt < 1 || bt > gff.numImproperTypes()) {
        std::cerr << RED << "This Impropertype (" << bt
                  << ") is not defined in the given "
                  << "parameter file\n"
                  << RESET;
        continue;
      }

      if (gff.improperType(bt - 1).q0() != 0) {
        std::cerr << RED << "You might want to modify the order of "
                  << " the atoms for this improper dihdral\n"
                  << "to make sure you get the correct stereochemistry\n"
                  << RESET;
      }
      improper.setType(bt - 1);
    }
    bbs->addImproper(improper);
    ++i;
  }
  return false;
}

bool PrepBBCore::addDihedral2BB(bool interact, const utils::FfExpert &exp,
                                const gcore::GromosForceField &gff,
                                std::vector<gcore::Dihedral> &newdihedrals) {
  std::vector<utils::FfExpert::counter> ocList;
  if (interact && !newdihedrals.empty()) {
    std::cerr << BLUE
              << "\n\n======= DIHEDRAL PARAMETERS ============================"
              << "======================\n\n"
              << RESET;
  }

  int i = 0;
  while (i < static_cast<int>(newdihedrals.size())) {
    auto &dihedral = newdihedrals[i];
    if (interact) {
      int i0 = dihedral[0];
      int i1 = dihedral[1];
      int i2 = dihedral[2];
      int i3 = dihedral[3];
      gcore::Dihedral iacdihedral(bbs->atom(i0).iac(), bbs->atom(i1).iac(),
                                  bbs->atom(i2).iac(), bbs->atom(i3).iac());
      std::cerr << BLUE
                << "\n--------------------------------------------------------"
                << "----------------------\n";
      std::cerr << "Considering dihedral " << i0 + 1 << " ("
                << bbs->atom(i0).name() << ")  -  " << i1 + 1 << " ("
                << bbs->atom(i1).name() << ")  -  " << i2 + 1 << " ("
                << bbs->atom(i2).name() << ")  -  " << i3 + 1 << " ("
                << bbs->atom(i3).name() << ")\n"
                << RESET;

      exp.iac2dihedral(iacdihedral, ocList);
      if (!ocList.empty()) {
        size_t ap = utils::sort(ocList);

        std::cerr << "\n"
                  << "Suggested Dihedraltype based on the IAC "
                  << "of the atoms (" << iacdihedral[0] + 1 << ", "
                  << iacdihedral[1] + 1 << ", " << iacdihedral[2] + 1 << ", "
                  << iacdihedral[3] + 1 << ") :\n";
        for (size_t i = 0; i < ocList.size(); ++i) {
          if (ap == i) {
            std::cerr << RED;
          }
          std::cerr << std::setw(10) << ocList[i].type + 1
                    << " : Kd = " << std::setw(6) << std::setprecision(1)
                    << gff.dihedralType(ocList[i].type).fc()
                    << "; pd = " << std::setw(4) << std::setprecision(1)
                    << gff.dihedralType(ocList[i].type).pd()
                    << "; m = " << std::setw(2) << std::setprecision(1)
                    << gff.dihedralType(ocList[i].type).np() << " occurs "
                    << ocList[i].occurence << " times\n";
          if (ap == i) {
            std::cerr << RESET;
          }
        }
        std::cerr << GREEN << std::setw(10) << "undo"
                  << " : In case you want to go back to the previous dihedral\n"
                  << RESET;
        std::cerr << YELLOW << std::setw(10) << "goback"
                  << " : In case you want to go back to Angle Parameters\n"
                  << RESET;
      }

      std::ostringstream prompt;
      prompt << "Give Dihedraltype ( 1 - " << gff.numDihedralTypes() << " ) : ";

      std::string line = cmd::getLine(prompt.str());
      if (line == "goback") {
        std::cerr << YELLOW << "\nUndo: going back to Angle Parameters\n\n"
                  << RESET;
        return true;
      } else if (line == "undo") {
        if (i > 0) {
          std::cerr << YELLOW << "\nUndo: going back to previous dihedral\n\n"
                    << RESET;
          --i;
        } else {
          std::cerr << YELLOW
                    << "\nAlready at the first dihedral, cannot go back\n\n"
                    << RESET;
        }
        continue;
      }

      std::istringstream s(line);
      int bt;
      s >> bt;
      if (s.fail()) {
        std::cerr << RED << "Invalid input. Try again.\n" << RESET;
        continue;
      }

      if (bt < 1 || bt > gff.numDihedralTypes()) {
        std::cerr << RED << "This Dihedraltype (" << bt << ") is not defined "
                  << "in the given parameter file\n"
                  << RESET;
        continue;
      }
      dihedral.setType(bt - 1);
    }
    bbs->addDihedral(dihedral);
    ++i;
  }
  return false;
}

class PrepBBLogic {
private:
  bool interact{false};
  bool has_graph{false};
  bool goBack{false};
  utils::FfExpert exp;
  utils::FfExpertGraph *graph{nullptr};
  gcore::GromosForceField *gff{nullptr};
  PrepBBCore *bbc{nullptr};
  std::set<int> aro;
  std::set<int> bt_aro;

  std::vector<std::string> atom_name;
  std::vector<int> order;
  std::map<int, int> newnum;

  std::vector<gcore::Bond> bond;
  std::vector<gcore::Bond> bondnew;
  std::vector<gcore::Angle> newangles;
  std::vector<gcore::Improper> newimpropers;
  std::vector<gcore::Dihedral> newdihedrals;

public:
  explicit PrepBBLogic(const args::Arguments &args);
  ~PrepBBLogic();

  void CreateGraph(const std::vector<std::string> &elements);
  void addAtom();
  void addBond();
  void addAngle();
  void addImproper();
  void addDihedral();
  void setBondedTerms();
  void clearBondedTerms();
  void printBB() { bbc->writeBB(); }

  bool shouldGoBack() const { return goBack; }
  void setGoBack() { goBack = true; }
  void clearGoBack() { goBack = false; }
};

PrepBBLogic::PrepBBLogic(const args::Arguments &args) {
  if (args.count("interact") >= 0) {
    interact = true;
  }

  bool suggest = false;
  if (args.count("build") > 0) {
    suggest = true;
    gcore::BuildingBlock mtb;

    auto iter = args.lower_bound("build");
    auto to = args.upper_bound("build");
    for (; iter != to; ++iter) {
      gio::InBuildingBlock ibb(iter->second);
      mtb.addBuildingBlock(ibb.building());
    }

    if (args.count("graph_library") && !args.count("pdb")) {
      utils::FfExpertGraphMapper mapper(args["graph_library"]);
      exp.learn(mtb, &mapper);
      has_graph = true;
    } else {
      exp.learn(mtb);
    }

    if (!interact) {
      throw gromos::Exception(program_name,
                              "specification of building block "
                              "and or parameter file for suggested "
                              "parameters only takes effect if interact "
                              "is specified");
    }
  }

  if (args.count("param") > 0) {
    gio::InParameter ipp(args["param"]);
    gff = new gcore::GromosForceField(ipp.forceField());
  } else if (suggest) {
    throw gromos::Exception(program_name,
                            "If you specify a building block "
                            "for suggestions, I also need a parameter file");
  }

  // read input
  std::vector<std::string> elements;
  std::string name;
  if (args.count("spec") > 0) {
    name = read_input(args["spec"], atom_name, elements, bond);
  } else if (args.count("pdb") > 0) {
    // let's disable graph suggestions for now. PDB doesn't know about
    // bond order
    has_graph = false;
    double bondbound = 0.0;
    if (args.count("bound") > 0) {
      bondbound = std::stod(args["bound"]);
    }
    name = read_pdb(args["pdb"], atom_name, elements, bond, bondbound);
  }

  if (!atom_name.empty() && bond.empty()) {
    throw gromos::Exception(program_name, "No bonds defining determined");
  }

  std::set<int> at;
  for (size_t i = 0; i < atom_name.size(); ++i) {
    at.insert(i);
  }

  // If necessary, get first atom and reorder the atoms
  if (args.count("reorder") > 0) {
    int first_atom = std::stoi(args["reorder"]) - 1;

    order.push_back(first_atom);
    at.erase(first_atom);
    add_neighbours(first_atom, bond, at, order);
  } else {
    for (size_t i = 0; i < atom_name.size(); ++i) {
      order.push_back(i);
    }
  }

  if (order.size() != atom_name.size()) {
    throw gromos::Exception(program_name, "ordering of atoms failed");
  }

  // do some bookkeeping
  for (size_t i = 0; i < order.size(); ++i) {
    newnum[order[i]] = i;
  }
  for (const auto &bd : bond) {
    gcore::Bond b(newnum[bd[0]], newnum[bd[1]]);
    bondnew.push_back(b);
  }

  if (interact && args.count("reorder") > 0) {
    std::cerr << "\n--------------------------------------------------------"
              << "----------------------\n";
    std::cerr << "Atoms have been renumbered starting with atom "
              << order[0] + 1 << " (" << atom_name[order[0]] << ")\n";
    std::cerr << "New atom list:\n";
    for (size_t i = 0; i < order.size(); ++i) {
      std::cerr << "\t" << std::setw(8) << i + 1 << std::setw(8)
                << atom_name[order[i]] << " (was " << order[i] + 1 << ")\n";
    }
    std::cerr << "\n";
    std::cerr << "Bonds have been adopted accordingly\n";
    std::cerr << "\n========================================================"
              << "======================\n";
  }

  // determine atoms that are in rings
  std::set<int> ra = ring_atoms(bondnew, order.size());

  if (interact && !ra.empty()) {
    std::cerr << "\n--------------------------------------------------------"
              << "----------------------\n";
    std::cerr << "Bonds indicate ring system involving atoms ";

    for (int atom : ra) {
      std::cerr << atom + 1 << " ";
    }
    std::cerr << "\n";
    if (cmd::getYesNo("Is this an aromatic ringsystem? (y/n) ")) {
      for (int atom : ra) {
        std::set<int> nb = neighbours(atom, bondnew);
        if (nb.size() < 4) {
          aro.insert(atom);
          for (int neighbour : nb) {
            if (!ra.count(neighbour)) {
              bt_aro.insert(neighbour);
            }
          }
        }
      }
    }

    if (cmd::getYesNo("Do you want to specify different atoms as being part of "
                      "an aromatic ring? (y/n) ")) {
      int number = 0;
      cmd::getValue<int>(number, "Give number of aromatic atoms: ");

      std::cerr << "Give " << number << " atom numbers involving an "
                << "aromatic ring: ";
      for (int j = 0; j < number; ++j) {
        int n;
        std::ostringstream prompt;
        prompt << " - " << (j + 1) << " atom: ";
        cmd::getValue<int>(n, prompt.str());
        aro.insert(n - 1);
      }
      for (int atom : aro) {
        for (int neighbour : neighbours(atom, bondnew)) {
          if (!aro.count(neighbour)) {
            bt_aro.insert(neighbour);
          }
        }
      }
    }
    std::cerr << BLUE
              << "========================================================"
              << "======================" << RESET;
  }

  bbc = new PrepBBCore(name);
  CreateGraph(elements);
}

PrepBBLogic::~PrepBBLogic() {
  if (interact) {
    std::cerr << BLUE
              << "\n\n======= PREPBB HAS GATHERED ALL INFORMATION AND WRITTEN"
              << " A BUILDING BLOCK=====" << RESET;
  }

  delete graph;
  delete bbc;
  delete gff;
}

void PrepBBLogic::CreateGraph(const std::vector<std::string> &elements) {
  if (has_graph) {
    graph = new utils::FfExpertGraph(atom_name, elements, bond, order);
  }
}

void PrepBBLogic::addAtom() {
  bbc->addAtom2BB(order, atom_name, interact, has_graph, exp, graph, *gff, bond,
                  newnum, aro, bt_aro);
}

void PrepBBLogic::addBond() {
  bool goBack = bbc->addBond2BB(interact, exp, *gff, bondnew);
  if (goBack) {
    setGoBack();
  }
}

void PrepBBLogic::addAngle() {
  bool goBack = bbc->addAngle2BB(interact, exp, *gff, newangles);
  if (goBack) {
    setGoBack();
  }
}

void PrepBBLogic::addImproper() {
  bool goBack = bbc->addImproper2BB(interact, exp, *gff, newimpropers);
  if (goBack) {
    setGoBack();
  }
}

void PrepBBLogic::addDihedral() {
  bool goBack = bbc->addDihedral2BB(interact, exp, *gff, newdihedrals);
  if (goBack) {
    setGoBack();
  }
}

void PrepBBLogic::clearBondedTerms() {
  newangles.clear();
  newimpropers.clear();
  newdihedrals.clear();
}

void PrepBBLogic::setBondedTerms() {
  for (size_t i = 0; i < order.size(); i++) {
    std::set<int> nb = neighbours(i, bondnew);
    std::vector<int> vnb(nb.begin(), nb.end());

    const auto add_angle_combinations = [&](const std::vector<int> &v) {
      for (size_t a = 0; a < v.size(); ++a) {
        for (size_t b = a + 1; b < v.size(); ++b) {
          newangles.emplace_back(v[a], i, v[b]);
        }
      }
    };

    switch (vnb.size()) {
    case 1:
      break;
    case 2:
      newangles.push_back(gcore::Angle(vnb[0], i, vnb[1]));
      break;
    case 3:
      add_angle_combinations(vnb);
      newimpropers.push_back(gcore::Improper(i, vnb[0], vnb[1], vnb[2]));
      break;
    case 4:
    case 5:
      add_angle_combinations(vnb);
      break;
    default:
      std::ostringstream os;
      os << "Don't know how to create angles for 0 or 6 bonds to atom "
         << i + 1;

      throw gromos::Exception(program_name, os.str());
    }
  }

  for (const auto &bond : bondnew) {
    int i0 = bond[0];
    int i1 = bond[1];

    std::set<int> nb1 = neighbours(i0, bondnew);
    nb1.erase(i1);
    std::set<int> nb2 = neighbours(i1, bondnew);
    nb2.erase(i0);

    if (!nb1.empty() && !nb2.empty()) {
      if (aro.count(i0) && aro.count(i1)) {
        // we are in an aromatic ring
        int a = -1, b = -1;
        for (int n : nb1) {
          if (aro.count(n)) {
            a = n;
          }
        }
        for (int n : nb2) {
          if (aro.count(n)) {
            b = n;
          }
        }
        if (a == -1 || b == -1) {
          throw gromos::Exception(program_name,
                                  "Error trying to determine "
                                  "improper dihedrals for aromatic ring");
        }
        newimpropers.push_back(gcore::Improper(a, i0, i1, b));
      } else {
        newdihedrals.push_back(
            gcore::Dihedral(*nb1.begin(), i0, i1, *nb2.begin()));
      }
    }
  }
}

void PrepBBCore::writeBB() {
  std::ofstream fout("BUILDING.out");
  gio::OutBuildingBlock obb(fout);
  obb.writeSingle(*bbs, gio::OutBuildingBlock::BBTypeSolute);
  fout.close();
  std::cerr << "\nBuilding block was written to BUILDING.out" << std::endl;
}

class State {
public:
  virtual ~State() = default;
  virtual void handle(class Context *ctx) = 0;
};

class AtomState : public State {
public:
  void handle(Context *ctx) override;
};

class BondState : public State {
public:
  void handle(Context *ctx) override;
};

class AngleState : public State {
public:
  void handle(Context *ctx) override;
};

class ImproperState : public State {
public:
  void handle(Context *ctx) override;
};

class DihedralState : public State {
public:
  void handle(Context *ctx) override;
};

class Context {
private:
  State *atomState{nullptr};
  State *bondState{nullptr};
  State *angleState{nullptr};
  State *improperState{nullptr};
  State *dihedralState{nullptr};

  State *state;
  PrepBBLogic *pbb{nullptr};
  bool exitRequested{false};

public:
  explicit Context(PrepBBLogic *prep)
      : pbb{prep}, atomState{new AtomState}, bondState{new BondState},
        angleState{new AngleState}, improperState{new ImproperState},
        dihedralState{new DihedralState} {
    state = atomState;
  }
  ~Context() {
    delete atomState;
    delete bondState;
    delete angleState;
    delete improperState;
    delete dihedralState;
  }

  void run() {
    while (!exitRequested) {
      state->handle(this);
    }
  }

  void setState(State *s) { state = s; }
  PrepBBLogic *getPrepBB() const { return pbb; }
  State *getAtomState() const { return atomState; }
  State *getBondState() const { return bondState; }
  State *getAngleState() const { return angleState; }
  State *getImproperState() const { return improperState; }
  State *getDihedralState() const { return dihedralState; }
  void requestExit() { exitRequested = true; }
};

void AtomState::handle(Context *ctx) {
  ctx->getPrepBB()->addAtom();
  ctx->setState(ctx->getBondState());
}

void BondState::handle(Context *ctx) {
  ctx->getPrepBB()->addBond();

  if (ctx->getPrepBB()->shouldGoBack()) {
    ctx->getPrepBB()->clearGoBack();
    ctx->setState(ctx->getAtomState());
  } else {
    ctx->getPrepBB()->setBondedTerms();
    ctx->setState(ctx->getAngleState());
  }
}

void AngleState::handle(Context *ctx) {
  ctx->getPrepBB()->addAngle();
  if (ctx->getPrepBB()->shouldGoBack()) {
    ctx->getPrepBB()->clearGoBack();
    ctx->getPrepBB()->clearBondedTerms();
    ctx->setState(ctx->getBondState());
  } else {
    ctx->setState(ctx->getImproperState());
  }
}

void ImproperState::handle(Context *ctx) {
  ctx->getPrepBB()->addImproper();
  if (ctx->getPrepBB()->shouldGoBack()) {
    ctx->getPrepBB()->clearGoBack();
    ctx->setState(ctx->getAngleState());
  } else {
    ctx->setState(ctx->getDihedralState());
  }
}

void DihedralState::handle(Context *ctx) {
  ctx->getPrepBB()->addDihedral();
  if (ctx->getPrepBB()->shouldGoBack()) {
    ctx->getPrepBB()->clearGoBack();
    ctx->setState(ctx->getAngleState());
  } else {
    ctx->requestExit();
  }
}

int main(int argc, char *argv[]) {
  args::Argument_List knowns;
  knowns << "spec" << "pdb" << "bound" << "reorder" << "build" << "param"
         << "interact" << "graph_library";

  std::string usage = std::string("# ") + program_name;
  usage += "\n\t@spec         <specifications> OR\n";
  usage += "\t@pdb            <pdb-file>\n";
  usage += "\t@bound          <upper bound for bondlength (for pdb)>\n";
  usage += "\t@reorder        <first atom>\n";
  usage += "\t[@build         <building block file>]\n";
  usage += "\t[@param         <parameter file>]\n";
  usage += "\t[@interact]\n";
  usage += "\t[@graph_library <file>]\n";

  try {
    args::Arguments args(argc, argv, knowns, usage);

    PrepBBLogic pbb(args);

    Context machine(&pbb);

    machine.run();

    /* pbb.addAtom();
    pbb.addBond();
    pbb.setBondedTerms();
    pbb.addAngle();
    pbb.addImproper();
    pbb.addDihedral(); */

    pbb.printBB();
  } catch (const gromos::Exception &e) {
    std::cerr << e.what() << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

std::string read_input(std::string file, std::vector<std::string> &atom,
                       std::vector<std::string> &element,
                       std::vector<gcore::Bond> &bond) {
  std::vector<std::string> buffer;
  gio::Ginstream gin(file);
  gin.getblock(buffer);

  if (buffer[buffer.size() - 1].find("END") != 0) {
    throw gromos::Exception(program_name, "Input file " + gin.name() +
                                              " is corrupted. No END in " +
                                              buffer[0] + " block. Got\n" +
                                              buffer[buffer.size() - 1]);
  }

  if (buffer[0] != "ATOMS") {
    throw gromos::Exception(program_name,
                            "ATOMS block expected in file " + file);
  }

  for (size_t i = 1; i < buffer.size() - 1; i++) {
    std::string name, elem;
    std::istringstream is(buffer[i]);
    is >> name;
    atom.push_back(name);
    if (is >> elem) {
      element.push_back(elem);
    }
  }

  gin.getblock(buffer);
  if (buffer[0] != "BONDS") {
    throw gromos::Exception(program_name,
                            "BONDS block expected in file " + file);
  }

  int a, b;
  for (size_t i = 1; i < buffer.size() - 1; i++) {
    std::istringstream is(buffer[i]);
    is >> a >> b;
    if (is.fail()) {
      throw gromos::Exception(program_name, "bad bond in BOND block");
    }
    bond.push_back(gcore::Bond(a - 1, b - 1));
    int type;
    if (is >> type) {
      bond.back().setType(type);
    }
  }

  return gin.title();
}

std::string read_pdb(std::string file, std::vector<std::string> &atom,
                     std::vector<std::string> &element,
                     std::vector<gcore::Bond> &bond, double bondbound) {
  std::ifstream fin(file.c_str());
  std::string inPdbLine;
  std::string resName;
  double bondbound2 = bondbound * bondbound;

  std::vector<gmath::Vec> pos;

  while (!fin.eof()) {
    getline(fin, inPdbLine);
    if (inPdbLine.substr(0, 4) == "ATOM" ||
        inPdbLine.substr(0, 6) == "HETATM") {
      std::istringstream is(inPdbLine.substr(12, 5));
      std::string name;
      is >> name;
      atom.push_back(name);
      resName = inPdbLine.substr(17, 4);
      gmath::Vec v;
      v[0] = atof(inPdbLine.substr(30, 8).c_str());
      v[1] = atof(inPdbLine.substr(38, 8).c_str());
      v[2] = atof(inPdbLine.substr(46, 8).c_str());
      pos.push_back(v);

      // get the element and remove white space
      std::string elementStr = inPdbLine.substr(77, 2).c_str();
      std::string::size_type i;
      while ((i = elementStr.find(" ")) != std::string::npos)
        elementStr.erase(i);

      element.push_back(elementStr);
    }
    if (inPdbLine.substr(0, 6) == "CONECT") {
      std::istringstream is(inPdbLine.substr(6, 80));
      int i, j;
      is >> i >> j;
      bond.push_back(gcore::Bond(i - 1, j - 1));
    }
  }
  if (bond.size() == 0) {
    for (unsigned int i = 0; i < atom.size(); i++) {
      for (unsigned int j = i + 1; j < atom.size(); j++) {
        if ((pos[i] - pos[j]).abs2() < bondbound2) {
          bond.push_back(gcore::Bond(i, j));
        }
      }
    }
  }

  return resName;
}

/*
string read_mol(string file, vector<string> &atom, vector<string> & element,
vector<Bond> & bond)
{
  ifstream fin(file.c_str());
  vector<string> buffer;
  for(unsigned int i = 0; !fin.eof(); ++i) {
    string line;
    getline(fin, line);
    if (i > 3) // skip the header
      buffer.push_back(line);
  }
  fin.close();

  vector<string>::const_iterator it = buffer.begin(),
          to = buffer.end();

  istringstream _lineStream(*it);
  unsigned int num_atoms, num_bonds;
  _lineStream >> num_atoms >> num_bonds;
  if (_lineStream.fail())
    throw gromos::Exception(program_name, "Could not read number of atoms/bonds
from molfile.");

  for(unsigned int i = 0; i < num_atoms; ++i) {
    _lineStream.str(*(++it));
    double pos;
    std::string e;
    _lineStream >> pos >> pos >> pos >> e; // discard the position information
    if (_lineStream.fail())
      throw gromos::Exception(program_name, "bad atom in molfile.");

    element.push_back(e);
    atom.push_back(e + (i + 1));
  }

  for(unsigned int i = 0; i < num_bonds; ++i) {
    _lineStream.str(*(++it));
    unsigned int i, j, t;
    _lineStream >> i >> j >> t;
    if (_lineStream.fail())
      throw gromos::Exception(program_name, "bad bond in molfile.");
    i--; j--;
    if (i < 0 || i >= atom.size() || j < 0 || j >= atom.size())
      throw gromos::Exception(program_name, "bad atoms in bond in molfile.");

    Bond bond(i,j);
    bond.setType(t);
    element.push_back(bond);
  }
  return "CMPD";
}
*/

std::set<int> neighbours(int a, const std::vector<gcore::Bond> &bond) {
  std::set<int> tmp;
  for (unsigned int i = 0; i < bond.size(); i++) {
    if (bond[i][0] == a) {
      tmp.insert(bond[i][1]);
    }
    if (bond[i][1] == a) {
      tmp.insert(bond[i][0]);
    }
  }
  return tmp;
}

void forward_neighbours(int a, std::vector<gcore::Bond> &bond,
                        std::set<int> &cum, int prev) {
  std::set<int> tmp = neighbours(a, bond);
  tmp.erase(prev);
  for (std::set<int>::iterator it = tmp.begin(), to = tmp.end(); it != to;
       ++it) {
    // cout << "\tadding " << *it << " (from " << a << ")" << endl;
    if (cum.count(*it)) {
      // cout << " is this a ring! atom? " << *it << endl;
      return;
    } else {
      cum.insert(*it);
      forward_neighbours(*it, bond, cum, a);
    }
  }
}

std::set<int> ring_atoms(std::vector<gcore::Bond> &bond, int numatoms) {
  std::set<int> ra;
  for (int i = 0; i < numatoms; i++) {
    // cout << "analysing ring-possibility for " << i << endl;

    std::set<int> tmp;
    forward_neighbours(i, bond, tmp, -1);
    if (tmp.count(i))
      ra.insert(i);
  }
  return ra;
}

void add_neighbours(int i, std::vector<gcore::Bond> &bond, std::set<int> &at,
                    std::vector<int> &order) {
  // cout << "called add_neighbours with i: " << i << endl;

  std::set<int> nb = neighbours(i, bond);
  // cout << "\t" << i << " has " << nb.size() << " neighbours" << endl;

  std::map<int, int> cb;
  std::set<int>::const_iterator it = nb.begin(), to = nb.end();
  for (; it != to; ++it) {
    // cout << "\t\tcounting the bonds for every neighbour" << endl;
    std::set<int> counted;
    counted.insert(i);

    if (at.count(*it))
      cb[*it] = count_bound(*it, bond, at, counted);
    else
      cb[*it] = bignumber;

    // cout << "\t\t" << *it << "\t" << cb[*it] << endl;
  }
  for (std::map<int, int>::iterator mit = cb.begin(), mto = cb.end();
       mit != mto; ++mit) {
    if (mit->second == bignumber)
      nb.erase(mit->first);
  }
  // cout << "\tcorrected neighbour size " << nb.size() << endl;

  while (nb.size()) {
    int min = bignumber;
    int nextatom = -1;

    for (it = nb.begin(), to = nb.end(); it != to; ++it) {
      if (cb[*it] < min) {
        min = cb[*it];
        nextatom = *it;
      }
    }
    if (nextatom != -1) {
      if (at.count(nextatom)) {
        nb.erase(nextatom);
        at.erase(nextatom);
        // cout << "added atom " << nextatom << endl;
        order.push_back(nextatom);
        add_neighbours(nextatom, bond, at, order);
      } else {
        nb.erase(nextatom);
        // cout << "ring! atom was gone already" << i << " - " << nextatom <<
        // endl;
      }
    }
  }
}

int count_bound(int i, std::vector<gcore::Bond> &bond, std::set<int> &at,
                std::set<int> &j) {
  std::set<int> nb = neighbours(i, bond);
  int counter = 0;
  std::set<int>::const_iterator it = nb.begin(), to = nb.end();
  for (; it != to; ++it) {
    if (at.count(*it) && !j.count(*it)) {
      counter++;
      j.insert(i);
      counter += count_bound(*it, bond, at, j);
    }
  }
  return counter;
}
