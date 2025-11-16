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

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   AmberTopology.h
 * Author: bschroed
 *
 * Created on March 7, 2018, 3:19 PM
 */

#ifndef AMBERTOPOLOGY_H
#define AMBERTOPOLOGY_H
//for .h - here
#include <map>
#include <cassert>
#include <iostream>
#include <sstream>
#include <vector>
#include <set>
#include <cstdio>
#include <cmath>

#include "Ginstream.h"

#include "../gcore/GromosForceField.h"
#include "../gcore/LinearTopology.h"
#include "../gcore/GromosForceField.h"
#include "../gcore/LJException.h"
#include "../gcore/LJType.h"
#include "../gcore/Exclusion.h"
#include "../gcore/AtomTopology.h"
#include "../gcore/MassType.h"
#include "../gcore/BondType.h"
#include "../gcore/Bond.h"
#include "../gcore/Angle.h"
#include "../gcore/AngleType.h"
#include "../gcore/Improper.h"
#include "../gcore/ImproperType.h"
#include "../gcore/Dihedral.h"
#include "../gcore/DihedralType.h"
#include "../gcore/CrossDihedral.h"
#include "../gcore/LinearTopology.h"
#include "../gmath/Physics.h"
#include "../utils/StringOps.h"
#include "../gromos/Exception.h"


using namespace std;
using namespace gcore;

namespace gio{
    
    class AmberTopology  : public Ginstream {
            AmberTopology &operator=(const AmberTopology&);
            
      public:
        GromosForceField d_gff;
        map<string, vector<string>> d_blocks;

          /**
         * Give Constructor the path<string> to amber top file!
         */
        AmberTopology();
        AmberTopology(const AmberTopology&);
        AmberTopology(string s);  // : d_blocks()
        virtual ~AmberTopology();

        /**
         * the init function reads in the whole file into the map of blocks and
         * reads in the topology version
         */
        void init();
        /**
         * parseForceField takes all blocks that end up in the forcefield
         * and stores the information in... d_gff
         */
        void parseFile(LinearTopology &lt, double ljscaling, bool atomic_chargegroups);

        const GromosForceField & forceField() const {
          return d_gff;
        }

        private:
            static string readAmber(string inputfile);

    };
}

#endif /* AMBERTOPOLOGY_H */

