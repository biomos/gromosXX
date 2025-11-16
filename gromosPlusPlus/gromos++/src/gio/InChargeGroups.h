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
 * File:   InChargeGroups.h
 * Author: bschroed
 *
 * Created on March 8, 2018, 5:57 PM
 */

#ifndef INCHARGEGROUPS_H
#define INCHARGEGROUPS_H

#include <vector>
#include <map>
#include <set>

#include "../gromos/Exception.h"

#include "Ginstream.h"
#include "../gcore/System.h"
#include "../gcore/Molecule.h"
#include "../gcore/MoleculeTopology.h"
#include "../gcore/AtomTopology.h"

//#include "../gcore/LinearTopology.h"


using namespace std;
using namespace gcore;

namespace gio{
    class InChargeGroups : public Ginstream  {
         InChargeGroups &operator=(const InChargeGroups&);


    public:
        string inFilePath;
        map<string, vector<string>> d_blocks;
        map<string, map< int, pair<int, set<string>> > > resChargeGroups; //residue, map<chargegroup, set<atoms>>
        InChargeGroups();
        InChargeGroups(const InChargeGroups& orig);
        InChargeGroups(string inChargeGroupFilePath);
        
        virtual ~InChargeGroups();
    
        System mapChargeGroupsOnResidues(System sys);    
        
    private:
        set<string> aminoAcidSet;
        map<int, pair<string,  double> > spreadChargeInGroup(map<int, pair<string,  double> >  AtomCharge, double chargeSum, int atomsNum, int resID, string resName, int chargeGroupID, int ChargeGroupCharge = 0  , set<string> exceptions = set<string>());
        MoleculeTopology mapChargeGroupsOnAtoms(map< int, pair<int, set<string>> > AtomCharges, int resID, MoleculeTopology mTopo);
        MoleculeTopology writeToMoleculeTopo (MoleculeTopology mTopo, map<int, pair<string,  double> > chargeGroupAtoms);
        void parseChargeGroupFile();
        void getChargeGroups();

    };
}
#endif /* INCHARGEGROUPS_H */

