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
#include "Hbond.h"

#include <cassert>
#include <cstdlib>
#include <string>
#include <iostream>
#include <sstream>
#include <vector>

#include "../args/Arguments.h"
#include "../gcore/System.h"
#include "CubeSystem.hcc"
#include "Hbond_calc.h"
#include "Hbond_calc_2c.h"
#include "Hbond_calc_3c.h"


using utils::HB;

HB::HB(gcore::System &sys, args::Arguments &args, HBPara2c hbparas2c, HBPara3c hbparas3c, int dummyIAC)
    :   do3c( args.count("threecenter") >= 0 ? true : false ),
        do_native(false),
        sort_occ( args.count("sort") >= 0 ? true : false ),
        doBridges( args.count("solventbridges") >= 0 ? true : false ),
        reduce( args.count("reducesolvent") >= 0 ? true : false ),
        vacuum(false),
        hb2c_calc(hbparas2c, reduce),
        hb3c_calc(hbparas3c, reduce),
        hb_bridges(reduce),
        hbpara2c(hbparas2c),
        hbpara3c(hbparas3c){

    args::Arguments::const_iterator it_arg = args.lower_bound("pbc");
    std::istringstream is(it_arg->second);

    if(is.str() == "v")
        vacuum = true;

    if (args.count("higherthan") > 0){
        it_arg = args.lower_bound("higherthan");
        std::istringstream iss(it_arg->second);
        iss >> higherthan;
    }
    else
        higherthan = 0;
//cout << "gefore hb2c_calc.setval" << endl;
    hb2c_calc.setval(sys, args, dummyIAC);

    if (doBridges)
        hb_bridges.setval(hb2c_calc,sys, args); //this references the donors, acceptors, num_A_xx to hb2c_calc

    if (do3c)
        hb3c_calc.setval(hb2c_calc,sys, args);

}//end HB::HB()

void HB::clear() {

  hb2c_calc.clear();
  if(do3c)
    hb3c_calc.clear();
  if(doBridges)
    hb_bridges.clear();
}//end HB::clear()

void HB::prepare_native(CubeSystem<int>& cubes_donors,CubeSystem<int>& cubes_acceptors, CubeSystem<Key2c>& cubes_bridges){
    if(!vacuum)
        hb2c_calc.calc_hb(cubes_donors, cubes_acceptors); //calculate hbonds from the reference structure
    else
        hb2c_calc.calc_vac();
    hb2c_calc.store_index(); //store all the hbond indexes of the native hbonds

    if (do3c){
        if(!vacuum)
            hb3c_calc.calc_hb(cubes_donors, cubes_acceptors);
        else
            hb3c_calc.calc_vac();
        hb3c_calc.store_index();
    }
    if (doBridges){
        if(!vacuum)
            hb_bridges.calc_hb(hb2c_calc, cubes_bridges);
        else
            hb_bridges.calc_vac(hb2c_calc);
        hb_bridges.store_index();
    }
    do_native = true;
    clear();
}

void HB::calc(CubeSystem<int>& cubes_donors,CubeSystem<int>& cubes_acceptors, CubeSystem<Key2c>& cubes_bridges) {
  if (do_native){
    hb2c_calc.calc_native();
    if(do3c)
      hb3c_calc.calc_native();
    if(doBridges)
        hb_bridges.calc_native();
  }
  else{
      if(!vacuum){
        hb2c_calc.calc_hb(cubes_donors,cubes_acceptors);
        if(doBridges)
            hb_bridges.calc_hb(hb2c_calc, cubes_bridges); //hb_bridges need their own cubesystem
        if(do3c)
            hb3c_calc.calc_hb(cubes_donors,cubes_acceptors);
      }
      else{
        hb2c_calc.calc_vac();
        if(doBridges)
            hb_bridges.calc_vac(hb2c_calc);
        if(do3c)
            hb3c_calc.calc_vac();
      }
  }
}//end HB::calc()

void HB::settime(double times) {
  hb2c_calc.settime(times);
  if (do3c)
    hb3c_calc.settime(times);
  if(doBridges)
    hb_bridges.settime(times);
}//end HB::settime()

void HB::printstatistics() {


  std::cout << "#\n"
          << "# 2-Centered hydrogen bond D-H..A counted if:\n"
          << "#     Distance H..A is at most " << hbpara2c.maxdist << "\n"
          << "#     Angle    D-H..A is at least " << hbpara2c.minangle << "\n";
  std::cout << "#\n#\n";
  if (do3c) {
    std::cout << "#\n"
            << "# 3-Centered hydrogen bond D-H..A1\n"
            << "#                              \\A2 counted if:\n"
            << "#     Distance H..A1 is at most " << hbpara3c.maxdist << "\n"
            << "#     Distance H..A2 is at most " << hbpara3c.maxdist << "\n"
            << "#     Angle    D-H..A1 is at least " << hbpara3c.minangle << "\n"
            << "#     Angle    D-H..A2 is at least " << hbpara3c.minangle << "\n"
            << "#     Sum of angles D-H..A1, D-H..A2, A1..H..A2 is at least "
            << hbpara3c.minanglesum << "\n"
            << "#     Dihedral angle D..A1..A2..H is at most " << hbpara3c.maxdihedral << "\n";
    std::cout << "#\n#\n";
}
  if(doBridges){
    std::cout << "#\n"
          << "# Solute-Solvent-Solute bridges counted if:\n"
          << "#     Same parameters as 2-centered H-bonds\n"
          << "#     The solvent molecule must form a H-bond-bridge (at least 1 atom) between two solute atoms\n"
          << "#     There are 4 possibilities to form such a bridge: (solute..solvent..solute)\n"
          << "#     1.    A1 .. H .. A2\n"
          << "#     2. D1-H1 .. A .. H2-D2\n"
          << "#     3.    A1 .. Hx-D-Hy .. A2\n"
          << "#     4. D1-H1 .. A-H .. A2\n"
          << "#\n#\n";

  }
  hb2c_calc.printstatistics(sort_occ, higherthan);
  if (do3c)
    hb3c_calc.printstatistics(sort_occ, higherthan);
  if (doBridges)
    hb_bridges.printstatistics(sort_occ, higherthan);

}//end HB::printstatistics()

utils::HBPara2c HB::mk_hb2c_paras(const vector<double> &hbparas) {
  utils::HBPara2c hbparas2c;
  hbparas2c.maxdist = hbparas[DISTANCE];
  hbparas2c.maxdist2 = hbparas[DISTANCE] * hbparas[DISTANCE];
  hbparas2c.minangle = hbparas[ANGLE];
  return hbparas2c;
}// end HB::mk_hb2c_paras()

utils::HBPara3c HB::mk_hb3c_paras(const vector<double> &hbparas) {
  utils::HBPara3c hbparas3c;
  hbparas3c.maxdist = hbparas[DISTANCE];
  hbparas3c.maxdist2 = hbparas[DISTANCE] * hbparas[DISTANCE];
  hbparas3c.minangle = hbparas[ANGLE];
  hbparas3c.minanglesum = hbparas[SUM];
  hbparas3c.maxdihedral = abs(hbparas[DIHEDRAL]);
  return hbparas3c;
}// end HB::mk_hb3c_paras()

//merge input and output map
void HB::merge(utils::HB& input, int traj_num){
    hb2c_calc.merge(input.hb2c_calc, traj_num); //merge the output(=this) and input maps
    if(do3c)
        hb3c_calc.merge(input.hb3c_calc, traj_num);
    if(doBridges)
        hb_bridges.merge(input.hb_bridges, traj_num);
}

