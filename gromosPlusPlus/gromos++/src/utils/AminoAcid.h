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

#ifndef INCLUDED_AMINOACID
#define INCLUDED_AMINOACID

#include <string>
#include <map>
#include <vector>

namespace utils {

  struct gromosAminoAcid {
    std::string acid;
    std::string base;
    std::map<std::string, std::vector<std::string> > Hdonors;
    std::map<std::string, std::vector<std::string> > Hacceptors;
    double pKa;
    double pKb;
    double pKc;
  };

  class gromosAminoAcidLibrary {

  private:
    std::string version;
    std::map<std::string, gromosAminoAcid> lib;
    //std::map<std::string, std::vector<std::string> > Hdonors;
    //std::map<std::string, std::vector<std::string> > Hacceptors;

  public:
    void load(std::string &fname);
    void loadHardcoded45A4(void);
    void loadHardcoded53A6(void);
    void writeLibrary(std::ostream &os, std::string title = "");
    std::string pdb2acid(std::string PDBname);
    std::string pdb2base(std::string PDBname);
    double pKa(std::string PDBname);
    double pKb(std::string PDBname);
    double pKc(std::string PDBname);
    std::vector<std::string> rHdonors(std::string PDBname, std::string GROMOSname);
    std::vector<std::string> rHacceptors(std::string PDBname, std::string GROMOSname);
    std::map<std::string, gromosAminoAcid> getAminoAcids();
  };

}
#endif
