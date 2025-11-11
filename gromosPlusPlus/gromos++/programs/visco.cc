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
 * @file visco.cc
 * calculates the bulk and shear viscosities
 */

/**
 * @page programs Program Documentation
 *
 * @anchor visco
 * @section visco calculates bulk and shear viscosities
 * @author @ref bh
 * @date 19-1-2010
 *
 * Program visco calculates the bulk and shear viscosities from the elements of
 * the pressure tensor that are written to MD energy trajectory files.
 * In order to access this data, visco makes use of the ene_ana library.
 * To obtain more accurate results from the simulation, the Einstein relation is
 * used instead of the direct evaluation of the autocorrelation function (Green-Kubo formulation).
 * Consider @f$P_{\alpha \beta}@f$ as being an element of the pressure tensor.
 * Consider that @f$G_{\alpha \beta}(t)@f$ is the integral of @f$P_{\alpha \beta}dt@f$:
 *
 * @f[ G_{\alpha \beta}(t) = \int \limits_0^t P_{\alpha \beta}(t)dt @f].
 *
 * We define @f$\eta_{\alpha \beta}@f$ as a viscosity term calculated in terms of the integral of the pressure component (@f$G_{\alpha \beta}@f$).
 * It will be proportional to the mean square "displacements" of @f$G_{\alpha \beta}(t)@f$
 * in the limit of infinit time.
 *
 * @f[ \eta_{\alpha \beta} = \frac{V}{2 k_B T} \lim_{t\to\infty} \langle\frac{[G_{\alpha \beta}(t+\tau) - G_{\alpha \beta}(t)]^2}{\tau}\rangle @f],
 *
 * where V is the volume of the periodic box, @f$k_B@f$ is the Boltzmann constant
 * and @f$T@f$ is the absolute temperature of the system.
 * For isotropic systems, the estimation of the bulk viscosities can be obtained from
 * the average of the viscosity terms obtained from the diagonal components of the pressure
 * tensor:
 *
 * @f[ \eta_{bulk} = (\eta_{xx}+\eta_{yy}+\eta_{zz})/3 @f].
 *
 * The shear viscosities of an isotropic system can be estimated by averaging the
 * viscosity terms obtained from the off-diagonal elements of the pressure tensor:
 *
 * @f[ \eta_{shear} = (\eta_{xy}+\eta_{xz}+\eta_{yz})/3 @f].
 *
 * The time series of the mean square "displacements" of @f$G_{\alpha \beta}(t)@f$
 * are printed to separate files (Gxx_msd.dat, Gyy_msd.dat, Gzz_msd.dat, Gxy_msd.dat, Gxz_msd.dat, Gyz_msd.dat).
 * In view of the poor statistics for long times, it is up to the user to decide the interval
 * for which the least square fitting should be calculated.
 * For convenience, program visco also prints the constant @f$\frac{V}{2 k_B T}@f$ and the conversion factors for the commonly used units.
 *
 *
 *
 *  <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@time</td><td>&lt;time and dt&gt; </td></tr>
 * <tr><td> \@en_files</td><td>&lt;energy files files&gt; </td></tr>
 * <tr><td> \@temp</td><td>&lt;temperature&gt; </td></tr>
 * <tr><td> \@library</td><td>&lt;library file (same as for ene_ana)&gt; </td></tr>
 * </table>
 *
 * 
 * Example:
 * @verbatim
  visco
    [@time      0 0.1]
    @en_files  ex.tre
    @temp      298.15
    @library   ene_ana.lib
   @endverbatim
 *
 */

#include <cassert>
#include <cctype>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <cstdlib>
#include <string>
#include <vector>
#include <map>
#include <cmath>

#include "../src/args/Arguments.h"
#include "../src/gio/Ginstream.h"
#include "../src/gmath/Stat.h"
#include "../src/gmath/Physics.h"
#include "../src/utils/EnergyTraj.h"
#include "../src/gromos/Exception.h"

using namespace std;
using namespace args;
using namespace gio;

void set_library(utils::EnergyTraj &e, string type);
void set_standards(utils::EnergyTraj &e, string type);
void read_library(string name, utils::EnergyTraj& e);
void calcvisco(gmath::Stat<double> &p, string s, vector<double> & time);


int main(int argc, char **argv) {

  Argument_List knowns;
  knowns << "time" << "en_files" << "temp" << "library";

  string usage = "# " + string(argv[0]);
  usage += "\n\t[@time    <time and dt>]\n";
  usage += "\t@en_files   <energy files>\n";
  usage += "\t[@temp   <temperature K; default 298>]\n";
  usage += "\t@library    <ene_ana library>\n";

  

  try {

    Arguments args(argc, argv, knowns, usage);

    // get simulation time either from the user or from the files
    vector<double> timearg = args.getValues<double>("time", 2, false,
          Arguments::Default<double>() << 0.0 << 1.0);
    double t0 = timearg[0];
    double dt = timearg[1];
    bool usertime = false;
    if (args.count("time") > 0) usertime = true;

    vector<double> time;

    if(args.count("en_files") <= 0)
      throw gromos::Exception("visco", "no data specified:\n" + usage);

    // read a library file?
    string library="gromos96";
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

    vector<string> prop;
    prop.push_back("boxvol");
    prop.push_back("diagP1");
    prop.push_back("diagP2");
    prop.push_back("diagP3");
    prop.push_back("offP1");
    prop.push_back("offP2");
    prop.push_back("offP3");
    prop.push_back("offP4");
    prop.push_back("offP5");
    prop.push_back("offP6");

    double Kb = gmath::physConst.get_boltzmann();

    int num_prop = 10;

    double temp = args.getValue<double>("temp", true);

    utils::EnergyTraj etrj;
    double mass=0;
    etrj.addConstant("MASS",mass);

    // learn about the variable names how they map to the elements
    read_library(library, etrj);

    if(print_library) etrj.write_map();

    // prepare for the statistical information
    vector<gmath::Stat<double> > s(num_prop);

    Ginstream gin_en;
    bool do_energy_files     =(args.count("en_files")>0);

    Arguments::const_iterator it_en=args.lower_bound("en_files"),
      to_en=args.upper_bound("en_files");

    if(do_energy_files) {
      gin_en.open(it_en->second.c_str());
      ++it_en;
    }

    while(true){

      // read the numbers into the energy trajectory
      if(do_energy_files){
	int end_of_file=etrj.read_frame(gin_en, "ENERTRJ");
	if(end_of_file){
	  if(it_en!=to_en){
	    gin_en.close();
	    gin_en.open(it_en->second.c_str());
	    ++it_en;
	    //try again...
	    etrj.read_frame(gin_en, "ENERTRJ");
	  }
	  else
	    break;
	}
      }
      // calculate and store the necessary number in the stat-classes
      for(int i=0; i<num_prop; i++)
	s[i].addval(etrj[prop[i]]);
      if(usertime)
	t0+=dt;
      else
	t0=etrj["TIME[2]"];
      time.push_back(t0);
      
    }

    // Print to screen
    cout << "### Visco: multiplication factor and conversion units ###" << endl;


    // Volume calculation
    double volume = s[0].ave();
    
    cout << endl << "Average volume:    "<< volume << endl;
    cout         << "Kb . T        :    "<< Kb * temp << endl;
    double factor_gromos = (0.5 * volume) / (Kb * temp);
    double factor_Punit = factor_gromos * 1.6605655E-5;
    cout << endl << "The multiplication of the angular coefficient of the msd by the "
            << "factor below will give the viscosity in standard unit." << endl;
    cout      <<    "V / (2 Kb T) =     " << factor_gromos << endl;
    cout << endl << "If the standard unit corresponds to the international units "
            << "then multiplication by the factor below will give the viscosity in 'Poise' (P)" <<endl;
    cout      <<    "(V / (2 Kb T)) * 1.6605655E-5 = " << factor_Punit << endl;

    for(int i=1; i<4; i++){
      calcvisco(s[i], prop[i], time); //Getting the diagonal Gcomponent
    }
        

    for(int i=4; i<7; i++){
      calcvisco(s[i], prop[i], time); //Getting the off diagonal Gcomponent
    }
    
  }

  catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

void calcvisco(gmath::Stat<double> &p, string s, vector<double> & time)
{
  if(p.n()!=int(time.size()))
    throw gromos::Exception("visco","number of time steps read is not equal"
			    " to number of data points to print out");

  ostringstream os;
  std::string name;
  bool check_name = false;
  // Assign names to the msd files:
  if (s == "diagP1"){
    name = "Gxx";
    check_name = true;
  }
  if (s == "diagP2"){
    name = "Gyy";
    check_name = true;
  }
  if (s == "diagP3"){
    name = "Gzz";
    check_name = true;
  }
  if (s == "offP1"){
    name = "Gxy";
    check_name = true;
  }
  if (s == "offP2"){
    name = "Gxz";
    check_name = true;
  }
  if (s == "offP3"){
    name = "Gyz";
    check_name = true;
  }

  if (check_name == true)
    os << name << "_msd.dat";
  else
    os << s << "_msd.dat";

  ofstream fout(os.str().c_str());
  fout.precision(9); //set precision of numbers going to ofstream
  fout << "#"
       << setw(14) << "time"
       << setw(15) << s
       << endl;

  vector<double> G;
  vector<double> DG;
  G.push_back(0);//vector to receive integral

  //Calculate integral
  for(int it=1; it< p.n(); it++){
    double dt = time[it]-time[it-1];
    double integ = G[it-1] + dt * p.val(it);
    G.push_back(integ);
  }

  for(int it=0; it< p.n(); it++){
    double sum = 0;
    int counter = 0;
    for(int j = 0; j < p.n() - it; j++) {
      counter++;
      const double d = G[j]-G[j + it];
      sum += d*d;
    }
    //const double disp = sum / counter ;
    const double disp = sum / counter;
    DG.push_back(disp);
    fout << setw(15) << time[it]
	 << setw(15) << disp
	 << endl;
  }
  fout.close();


  //FIT: Do not perform! It should be done by the user
  /*
  double sx = 0;
  double sy = 0;
  double sxx = 0;
  double sxy = 0;
  //for(int i=0; i< (p.n()/10); i++){
  for(int i=0; i< 100; i++){
      sx += time[i];
      sy += DG[i];
      sxx += time[i] * time[i];
      sxy += time[i] * DG[i];
  }
  */
  //double a = (sxy - sx * sy / p.n()) / (sxx - sx * sx / p.n());
  //cout << "Angular coefficient: " << a << endl;
  

  //double vis_component = a / (2 * kb);
  //return vis_component;
}



void set_library(utils::EnergyTraj &e, string type)
{
  if(type=="gromos96"){
    e.addBlock("  block TIMESTEP"                              , "ENERTRJ");
    e.addBlock("    subblock TIME 2 1"                         , "ENERTRJ");
    e.addBlock("  block ENERGY"                                , "ENERTRJ");
    e.addBlock("    subblock ENER 22 1"                        , "ENERTRJ");
    e.addBlock("    subblock ENERES 6 1"                       , "ENERTRJ");
    e.addBlock("    size  NUM_ENERGY_GROUPS"                   , "ENERTRJ");
    e.addBlock("    subblock ENERNB matrix_NUM_ENERGY_GROUPS 4", "ENERTRJ");
    e.addBlock("  block VOLUMEPRESSURE"                        , "ENERTRJ");
    e.addBlock("    subblock VOLPRT 20 1"                      , "ENERTRJ");
    e.addBlock("  block TIMESTEP"                              , "FRENERTRJ");
    e.addBlock("    subblock TIME 2 1 "                        , "FRENERTRJ");
    e.addBlock("  block FREEENERGYLAMBDA "                     , "FRENERTRJ");
    e.addBlock("    subblock ENER 9 1 "                        , "FRENERTRJ");
    e.addBlock("    subblock RLAM 1 1 "                        , "FRENERTRJ");
    e.addBlock("    subblock FREN  22 1 "                      , "FRENERTRJ");
  }
  else if(type=="gromosxx"){
    e.addBlock("  block TIMESTEP"                              , "ENERTRJ");
    e.addBlock("    subblock TIME 2 1"                         , "ENERTRJ");
    e.addBlock("  block ENERGY03"                              , "ENERTRJ");
    e.addBlock("    subblock ENER 16 1"                        , "ENERTRJ");
    e.addBlock("    size  NUM_BATHS"                           , "ENERTRJ");
    e.addBlock("    subblock KINENER NUM_BATHS 3"              , "ENERTRJ");
    e.addBlock("    size  NUM_ENERGY_GROUPS"                   , "ENERTRJ");
    e.addBlock("    subblock BONDED NUM_ENERGY_GROUPS 4"       , "ENERTRJ");
    e.addBlock("    subblock NONBONDED matrix_NUM_ENERGY_GROUPS 2", "ENERTRJ");
    e.addBlock("    subblock SPECIAL NUM_ENERGY_GROUPS 7"      , "ENERTRJ");
    e.addBlock("  block VOLUMEPRESSURE03"                      , "ENERTRJ");
    e.addBlock("    subblock MASS 1 1"                         , "ENERTRJ");
    e.addBlock("    size  NUM_BATHS"                           , "ENERTRJ");
    e.addBlock("    subblock TEMPERATURE NUM_BATHS 4"          , "ENERTRJ");
    e.addBlock("    subblock VOLUME 10 1"                      , "ENERTRJ");
    e.addBlock("    subblock PRESSURE 30 1"                    , "ENERTRJ");
    e.addBlock("  block TIMESTEP"                              , "FRENERTRJ");
    e.addBlock("    subblock TIME 2 1"                         , "FRENERTRJ");
    e.addBlock("  block FREEENERDERIVS03"                      , "FRENERTRJ");
    e.addBlock("    subblock RLAM  1 1"                        , "FRENERTRJ");
    e.addBlock("    subblock FREEENER 16 1"                    , "FRENERTRJ");
    e.addBlock("    size  NUM_BATHS"                           , "FRENERTRJ");
    e.addBlock("    subblock FREEKINENER NUM_BATHS 3"          , "FRENERTRJ");
    e.addBlock("    size  NUM_ENERGY_GROUPS"                   , "FRENERTRJ");
    e.addBlock("    subblock FREEBONDED NUM_ENERGY_GROUPS 4"       , "FRENERTRJ");
    e.addBlock("    subblock FREENONBONDED matrix_NUM_ENERGY_GROUPS 2", "FRENERTRJ");
    e.addBlock("    subblock FREESPECIAL NUM_ENERGY_GROUPS 7"      , "FRENERTRJ");

  }
}

void set_standards(utils::EnergyTraj &e, string type)
{
  e.addConstant("BOLTZ", gmath::physConst.get_boltzmann());

  if(type=="gromos96"){
    e.addKnown("time",   "TIME[2]");
    e.addKnown("totene", "ENER[1]");
    e.addKnown("totkin", "ENER[2]");
    e.addKnown("totpot", "ENER[9]");
    e.addKnown("totvdw", "ENER[18]");
    e.addKnown("totcrf", "ENER[19]");
    e.addKnown("boxvol", "VOLPRT[8]");
    e.addKnown("dE_tot", "FREN[1]");
    e.addKnown("dE_kin", "FREN[3]");
    e.addKnown("dE_pot", "FREN[9]");
    e.addKnown("dE_vdw", "FREN[18]");
    e.addKnown("dE_crf", "FREN[19]");
  }
  else if(type=="gromosxx"){
    e.addKnown("time", "TIME[2]");
    e.addKnown("E_tot", "ENER[1]");
    e.addKnown("E_kin", "ENER[2]");
    e.addKnown("E_pot", "ENER[3]");
    e.addKnown("E_vdw", "ENER[8]");
    e.addKnown("E_crf", "ENER[9]");
    e.addKnown("E_cov", "E_pot - E_vdw - E_crf");
    e.addKnown("E_special", "E_tot - E_pot - E_kin");
    e.addKnown("boxvol", "VOLUME[1]");
    e.addKnown("MASS", "MASS[1]");
    e.addKnown("dE_tot", "FREEENER[1]");
    e.addKnown("dE_kin", "FREEENER[2]");
    e.addKnown("dE_pot", "FREEENER[3]");
    e.addKnown("dE_vdw", "FREEENER[8]");
    e.addKnown("dE_crf", "FREEENER[9]");
    e.addKnown("dE_cov", "dE_pot - dE_vdw - dE_crf");
    e.addKnown("dE_special", "dE_tot - dE_pot - dE_kin");
  }
}

void read_library(string name, utils::EnergyTraj& e)
{
  Ginstream gin;

  try{

    gin.open(name);
  }
  catch (gromos::Exception ex){

    std::transform(name.begin(), name.end(), name.begin(), ::tolower);

    if(name=="gromos96"){
      set_library(e,"gromos96");
      set_standards(e, "gromos96");
    }else if(name=="gromosxx"){
      set_library(e,"gromosxx");
      set_standards(e, "gromosxx");
    }
    else
      throw gromos::Exception("read_library", "failed to open library file "
			      +name+
			      "\ngive gromos96 or gromosxx for standards");
    return;
  }
  while(true){

    vector<string> buffer;
    gin.getblock(buffer);
    if(gin.stream().eof()) break;
    if(buffer[buffer.size()-1].find("END")!=0)
      throw gromos::Exception("visco", "Library file " + gin.name() +
			      " is corrupted. No END in "+buffer[0]+
			      " block. Got\n"
			      + buffer[buffer.size()-1]);
    string sdum;

    if(buffer[0]=="ENERTRJ" || buffer[0]=="FRENERTRJ"){
      for(unsigned int i=1; i< buffer.size()-1; i++){
	e.addBlock(buffer[i], buffer[0]);

      }
    }

    vector<string> data;
    if(buffer[0]=="VARIABLES"){

      set_standards(e, "no");

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


