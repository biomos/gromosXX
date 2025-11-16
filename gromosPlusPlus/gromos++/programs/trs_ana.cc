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
 * @file trs_ana.cc
 * extracts time series from (energy) trajectory files
 */

/**
 * @page programs Program Documentation
 *
 * @anchor trs_ana
 * @section trs_ana analyse (energy) trajectories
 * @author @ref mp
 * @date 5. 2. 2016
 *
 * Program trs_ana extracts individual values from  gromos trajectory files and
 * can perform simple mathematical operations on them. 
 *
 * The program is based on @ref ene_ana.  
 * It uses the same library format to define blocks which can be read from any 
 * trajectory file that comes in the Gromos block-format. In contrast to 
 * ene_ana it does not require that all the blocks defined in the library 
 * are present in the trajectory file or in the
 * specified order. It can handle trajectories where not all timesteps contain 
 * the same number of blocks, e.g. when different properties were written to the
 * trajectory at different intervals. The time in the output timeseries will
 * always correspond to the time in the previous TIMESTEP block if no time is 
 * given by the user, else the time will be increased by the given timestep at 
 * every occurrence of a TIMESTEP block.
 * If multiple blocks of the same name occur between two TIMESTEPs, only the last 
 * one will be used.
 * 
 * In the library file one can also define properties to be calculated from 
 * the values that are listed in them. For the selected properties, trs_ana 
 * will calculate the time series, averages, root-mean-square fluctuations and 
 * a statistical error estimate. The error estimate is calculated from block 
 * averages of different sizes, as described in Allen and Tildesley: "Computer 
 * Simulation of Liquids", 1987. If a topology is supplied, the trs_ana 
 * uses this to define the total solute mass (MASS) and the total number of 
 * solute molecules (NUMMOL).
 * 
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@trs</td><td>&lt;trajectory files&gt; (and/or) </td></tr>
 * <tr><td> \@prop</td><td>&lt;properties to monitor&gt; </td></tr>
 * <tr><td> \@library</td><td>&lt;library for property names&gt; [print] </td></tr>
 * <tr><td> [\@topo</td><td>&lt;molecular topology file&gt; (for MASS and NUMMOL)] </td></tr>
 * <tr><td> [\@time</td><td>&lt;@ref utils::Time "time and dt"&gt; (overwrites TIME in the trajectory files)] </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  trs_ana
    @topo       ex.top
    @trs   ex.trs
    @prop       densit
    @library    trs.lib
    @time  0 2
   @endverbatim
 *
 * <hr>
 */
#include <cassert>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <string>
#include <vector>
#include <map>
#include <cmath>

#include "../src/args/Arguments.h"
#include "../src/gio/Ginstream.h"
#include "../src/gio/InTopology.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gmath/Stat.h"
#include "../src/gmath/Physics.h"
#include "../src/utils/EnergyTraj.h"
#include "../src/gromos/Exception.h"

using namespace std;
using namespace args;
using namespace gio;
using namespace gcore;

void print(gmath::Stat<double> &p, string s, vector<double> & time);
void set_standards(utils::EnergyTraj &e);
void read_library(string name, utils::EnergyTraj& e);

int main(int argc, char **argv){

  Argument_List knowns; 
  knowns << "topo" << "time" << "trs" << "prop" << "library";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@trs    <special trajectory file or any other time series in gromos block format>\n";
  usage += "\t@prop        <properties to monitor>\n";
  usage += "\t@library     <library for property names> [print]\n";
  usage += "\t[@topo       <molecular topology file> (for MASS and NUMMOL)]\n";
  usage += "\t[@time     <t and dt> (overwrites TIME in the trajectory files)]\n";

  try{
    Arguments args(argc, argv, knowns, usage);

    // get simulation time either from the user or from the files
    bool usertime=false;
    double t0=0, dt=1; 
    {
      Arguments::const_iterator iter=args.lower_bound("time");
      if(iter!=args.upper_bound("time")){
        t0=atof(iter->second.c_str());
        ++iter;
      }
      if(iter!=args.upper_bound("time")){
        dt=atof(iter->second.c_str());
        usertime=true;
        // a bit ugly: as the time is increased by dt before the first printout: reduce t0
        t0 -= dt;
      }
    }

    // check whether we are doing anything
    if(args.count("trs")<=0 )
      throw gromos::Exception("trs_ana", "no data specified:\n"+usage);
    if(args.count("prop") <=0)
      throw gromos::Exception("trs_ana", "no properties to follow:\n"+usage);
    
    // NEW: require a library 
    if(args.count("library") <=0)
      throw gromos::Exception("trs_ana", "no library file specified:\n"+usage);
    
    // read a library file?
    string library="";
    int print_library=0;
    {
      Arguments::const_iterator iter=args.lower_bound("library"), 
    to=args.upper_bound("library");
      if(iter!=to){
    library=iter->second;
    ++iter;
      }
      if(iter!=to && iter->second == "print") print_library=1;
    }

    // define an energy trajectory
    utils::EnergyTraj etrj;

    // read topology for the mass and the number of molecules
    double mass=0;
    if(args.count("topo")>0){
      
      InTopology it(args["topo"]);
      System sys(it.system());
      etrj.addConstant("NUMMOL", sys.numMolecules());
      // calculate total mass of the system
      for(int m=0; m<sys.numMolecules(); m++){
    for(int a=0; a< sys.mol(m).topology().numAtoms(); a++){
      mass += sys.mol(m).topology().atom(a).mass();
    }
      }
    }
    etrj.addConstant("MASS",mass);
    
    // learn about the variable names how they map to the elements
    read_library(library, etrj);
    
    if(print_library) etrj.write_map();
    
    // which properties do we follow
    vector<string> prop;
    int num_prop=0;
    {
      Arguments::const_iterator iter=args.lower_bound("prop"), 
    to=args.upper_bound("prop");
      while(iter!=to){
        if (etrj.find_property(iter->second)) prop.push_back(iter->second);
        else cerr << "# Warning: Property " << iter->second 
                  << " not found in the library\n";
      ++iter;
      }
      num_prop=prop.size();
      if (num_prop == 0)
            throw gromos::Exception("trs_ana","no valid properties specified");
    }

    // prepare for the statistical information
    vector<gmath::Stat<double> > s(num_prop);
    vector<vector<double> > time(num_prop);

    // define input stream
    Ginstream gin;
    
    Arguments::const_iterator it=args.lower_bound("trs"),
      to=args.upper_bound("trs");
    while (it != to) {
      gin.open(it->second.c_str()); 
      ++it; 
      
      // first we need a TIMESTEP
      std::vector<std::string> buffer;
      gin.getblock(buffer);
      if (buffer[0] == "TIMESTEP") {
          if(usertime)
            t0+=dt;
          else {
            etrj.read_block(buffer, "SPECIALTRJ");
            t0=etrj["TIME[2]"];          
          }
      } else {
        throw gromos::Exception("EnergyTraj", 
                      "The first block needs to be a TIMESTEP, found " + buffer[0]);
      }
      
      // until the file ends
      while(true) {
        gin.getblock(buffer);
        
        if (gin.stream().eof()) {
          // collect the last values of a file
          for(int i=0; i<num_prop; i++){
            utils::EnergyIndex ei = etrj.index(prop[i]);
            double value;
            bool value_exists = etrj.value_ifpossible(ei, value, prop[i]);
            if (value_exists) {
               s[i].addval(value);
               time[i].push_back(t0);
            }
          }  
          gin.close(); 
          break;
        } 
        
        // at every TIMESTEP calculate and store the necessary numbers (which 
        // can be calculated from the blocks found since the last TIMESTEP)
        //  in the stat-classes and then clear the data vector
        if (buffer[0] == "TIMESTEP") { 
          for(int i=0; i<num_prop; i++){
            utils::EnergyIndex ei = etrj.index(prop[i]);
            double value;
            bool value_exists = etrj.value_ifpossible(ei, value, prop[i]);
            if (value_exists) {
               s[i].addval(value);
               time[i].push_back(t0);
            }
          }        
       
          etrj.clear_data();
          
          if(usertime)
            t0+=dt;
          else {
            etrj.read_block(buffer, "SPECIALTRJ");
            t0=etrj["TIME[2]"];          
          }
        } else {
          etrj.read_block(buffer, "SPECIALTRJ");
        } 
      }
    }       
    //print out the statistical information
    cout << setw(12) << "# property"
     << " "
     << setw(15) << "average"
     << " "
     << setw(15) << "rmsd"
     << " "
     << setw(15) << "error est."
     << endl;
    for(int i=0; i<num_prop; i++)
     if (s[i].n() != 0) print(s[i], prop[i], time[i]);
     else  cout<< setw(12) << prop[i]<< setw(15) << "novalues" << setw(15)
               << "novalues"<< setw(15) << "novalues" << endl;
  }
  
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}



void print(gmath::Stat<double> &p, string s, vector<double>& time)
{
  if(p.n()!=int(time.size())) 
    throw gromos::Exception("trs_ana","number of time steps read is not equal"
                " to number of data points to print out");
  
  ostringstream os;
  os << s << ".dat";
  ofstream fout(os.str().c_str());
  fout.precision(9); //set precision of numbers going to ofstream
  fout << "#"
       << setw(14) << "time"
       << " "
       << setw(15) << s
       << endl;
  for(int i=0; i< p.n(); i++){
    fout << setw(15) << time[i]
     << " "
     << setw(15) << p.val(i)
     << endl;
  }
  fout.close();
// and print the averages etc to cout
  cout.precision(9); // set precision of number going to cout
  cout << setw(10) << s
       << " " //put an extra space, that will always pe printed, even if setw space does not suffice, to prevent merging of columns
       << setw(15) << p.ave()
       << " "
       << setw(15) << p.rmsd()
       << " "
       << setw(15) << p.ee()
       << endl;
}

void set_standards(utils::EnergyTraj &e)
{  
  e.addConstant("BOLTZ", gmath::physConst.get_boltzmann());
}

void read_library(string name, utils::EnergyTraj& e)
{
  Ginstream gin;
  
  try{
    
    gin.open(name);
  }
  catch (gromos::Exception ex){
      throw gromos::Exception("read_library", "failed to open library file "
                  +name);
  }
  while(true){
    
    vector<string> buffer;
    gin.getblock(buffer);
    if(gin.stream().eof()) break;
    if(buffer[buffer.size()-1].find("END")!=0)
      throw gromos::Exception("trs_ana", "Library file " + gin.name() +
                  " is corrupted. No END in "+buffer[0]+
                  " block. Got\n"
                  + buffer[buffer.size()-1]);
    string sdum;
    
    if(buffer[0]=="SPECIALTRJ"){
      for(unsigned int i=1; i< buffer.size()-1; i++){
        e.addBlock(buffer[i], buffer[0]);
      }
    }
    

    vector<string> data;
    if(buffer[0]=="VARIABLES"){
      
      set_standards(e);
      
      string bufferstring;
      
      gio::concatenate(buffer.begin()+1, buffer.end(), bufferstring);
      
      istringstream iss(bufferstring);

      // i am aware of the fact that END will also be stored in data.
      // This is used in parsing later on
      while(sdum!="END"){
    iss >> sdum;
    data.push_back(sdum);
      }
      
      // now search for the first appearance of "="
      for(unsigned int i=0; i< data.size(); i++){
    if(data[i]=="="){
      
      // search for either the next appearance or the end
      unsigned int to=i+1;
      for(; to < data.size(); to++) if(data[to]=="=") break;
      
      // parse the expression part into an ostringstream
      ostringstream os;
      for(unsigned int j=i+1; j< to-1; j++) os << " " << data[j]; 
      e.addKnown(data[i-1], os.str());
    }
      }
    }
  }
}
