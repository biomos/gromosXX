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
 * @file tser.cc
 * Calculates time series of properties
 */

/**
 * @page programs Program Documentation
 *
 * @anchor tser
 * @section tser Calculates time series of properties
 * @author @ref mc
 * @date 22. 11. 2004
 *
 * Program tser can calculate structural quantities from a trajectory file and 
 * print the time series and / or a distribution of the value associated with 
 * the requested property. The quantity to be calculated is specified through 
 * a @ref PropertySpecifier and can be any of the structural properties,
 * which can be calculated from atomic positions in the trajectory file. Time 
 * series can later be analysed further with e.g. the program @ref tcf.
 * 
 * Note that the keyword periodic (\@dist) can be used to map all values periodically
 * to the intervall between lower and upper (assuming (upper - lower) is a full
 * period length). This is useful e.g. for the calculation of torsional dihedral
 * angles including a distribution from 0 to 360 degrees. If the keyword periodic
 * is missing, the distribution is done omitting values outside the specified 
 * range.
 * 
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; [&lt;gathermethod&gt;] </td></tr>
 * <tr><td> \@time</td><td>&lt;@ref utils::Time "time and dt"&gt; </td></tr>
 * <tr><td> \@prop</td><td>&lt;@ref PropertySpecifier "properties"&gt; </td></tr>
 * <tr><td> [\@nots</td><td>(do not write time series)] </td></tr>
 * <tr><td> [\@dist</td><td>&lt;steps [min max] [periodic] </td></tr>
 * <tr><td> [\@norm</td><td>(normalise distribution)] </td></tr>
 * <tr><td> [\@solv</td><td>(read in solvent)] </td></tr>
 * <tr><td> \@traj</td><td>&lt;trajectory files&gt; </td></tr>
 * <tr><td> [\@skip</td><td>&lt;skip n first frames&gt;] </td></tr>
 * <tr><td> [\@stride</td><td>&lt;take every n-th frame&gt;] </td></tr>
 * </table>
 * 
 * <b>See also</b> @ref PropertySpecifier "property specifier"
 *
 * Example:
 * @verbatim
  tser
    @topo ex.top
    @pbc  r
    @time 0 0.1
    @prop t%1:1,3,5,6
    @traj ex.tr
    @skip 0
    @stride 1
 @endverbatim
 *
 * @bug Mar 22 2005: nearestImage calls in properties were missing
 * <hr>
 */

#include <cstdlib>
#include <string>
#include <vector>
#include <iomanip>
#include <cmath>
#include <iostream>
#include <cassert>
#include <sstream>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/fit/Reference.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/fit/PositionUtils.h"
#include "../src/utils/Property.h"
#include "../src/utils/PropertyContainer.h"
#include "../src/utils/groTime.h"
#include "../src/gmath/Stat.h"
#include "../src/gromos/Exception.h"

using namespace std;
using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;
using namespace utils;

int main(int argc, char **argv){

  Argument_List knowns;
  knowns << "topo" << "pbc" << "time" << "prop" << "traj" << "skip" << "stride"
         << "nots" << "dist" << "norm" << "solv";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo      <molecular topology file>\n";
  usage += "\t@pbc       <boundary type> [<gathermethod>]\n";
  usage += "\t@time      <time and dt>\n";  
  usage += "\t@prop      <properties>\n";
  usage += "\t[@nots     (do not write time series)]\n";
  usage += "\t[@dist     <steps [min max] [periodic]>]\n";
  usage += "\t[@norm     (normalise distribution)]\n";
  usage += "\t[@solv     (read in solvent)]\n";
  usage += "\t@traj      <trajectory files>\n";
  usage += "\t[@skip     <skip n first frames>]\n";
  usage += "\t[@stride   <take every n-th frame>]\n";
 
  try{
    Arguments args(argc, argv, knowns, usage);

    // get the @time argument
    utils::Time time(args);
    
    bool do_tser = true;
    if (args.count("nots") >= 0)
      do_tser = false;

    // if dist_min and dist_max are
    // not supplied, the maximum value
    // (of Stat) lies outside (!) the
    // distribution and is missed.
    // Still, the behaviour seems to be
    // correct.
    // probably it would be best to add
    // a tiny number to max in this case,
    // but then: what is a tiny number?
    bool do_dist = false;
    bool dist_boundaries = false;
    bool periodic = false;
    double dist_min=0, dist_max=0;
    int dist_steps=0; 
    if (args.count("dist") > 0)
    {
      do_dist = true;
      Arguments::const_iterator iter=args.lower_bound("dist");
      if(iter!=args.upper_bound("dist")){
	std::istringstream is(iter->second);
	is >> dist_steps;
	++iter;
	if (dist_steps <= 0){
	  throw Arguments::Exception("distribution: wrong number of steps specified" +
				     iter->second);
	}
      }
      if(iter!=args.upper_bound("dist")){
	dist_boundaries = true;
	std::istringstream is(iter->second);
	is >> dist_min;
	++iter;
      }
      if(iter!=args.upper_bound("dist")){
	std::istringstream is(iter->second);
	is >> dist_max;
	if (dist_max < dist_min){
	  throw Arguments::Exception("distribution: maximum value smaller than minimum value");
	}
        iter++;
      }
      if(iter!=args.upper_bound("dist")){
        if(iter->second == "periodic") {
          periodic = true;
        } else {
          stringstream msg;
          msg << "distribution: keyword " << iter->second << " not known";
          throw Arguments::Exception(msg.str());
        }
      }
    }

    bool normalize = false;
    if (args.count("norm") != -1)
      normalize = true;

    //  read topology
    args.check("topo",1);
    InTopology it(args["topo"]);
    
    System sys(it.system());
    System refSys(it.system());

    // parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(sys, args);
    // parse gather method
    Boundary::MemPtr gathmethod = args::GatherParser::parse(sys,refSys,args);

    // read in a property
    
    // it's nice to actually store the properties read in somewhere
    // let's use a PropertyContainer, as it can handle some of the
    // work later on
    // the PropertyContainer should know about the system, so it can
    // check the input whether ie the atoms exist in the topology
    PropertyContainer props(sys, pbc);
    {
      std::string prop;
      Arguments::const_iterator iter=args.lower_bound("prop");
      Arguments::const_iterator to=args.upper_bound("prop");
      // we read in all properties specified by the user
      if(iter==to)
	throw Arguments::Exception("no property given");
      for(; iter!=to; iter++){
	string spec=iter->second.c_str();
	prop += " " + spec;
	// props.addSpecifier(spec);
      }    
      // and that's how easy it is to add a standard property
      // like distance, angle, torsional angle
      props.addSpecifier(prop);
    }

    int skip = args.getValue<int>("skip", false, 0);
    int stride = args.getValue<int>("stride", false, 1);

    bool solvent = false;
    if (args.count("solv") != -1)
      solvent = true;

    // define input coordinate
    InG96 ic(skip, stride);

    // title
    if (do_tser){
      cout << "#" << endl;
      cout << "#" << setw(9) << "time";
      cout << "\t\t" << props.toTitle();
      cout << endl;
      std::cout.precision(6);
      std::cout.setf(std::ios::fixed, std::ios::floatfield);
    }
    
    // loop over all trajectories
    for(Arguments::const_iterator 
	  iter=args.lower_bound("traj"),
	  to=args.upper_bound("traj");
	iter!=to; ++iter){

      // open file
      ic.open((iter->second).c_str());
      if (solvent) ic.select("ALL");
      
      // loop over single trajectory
      while(!ic.eof()){
	ic >> sys >> time;
	if (ic.stride_eof()) break;
      
	(*pbc.*gathmethod)();
      
	// calculate the props
	// this is now the place, where a property-container is very handy
	// it knows that we want to loop over all properties and calculate
	// their 'value'.
	props.calc();

	// print the properties
	// this is a time series, so let's just print out the properties
	// the << operator is overloaded for the property container as well
	// as for the single properties
	if (do_tser){
	  cout << time << "\t\t";
	  cout << props;
	}
      }
      ic.close();
    }

    
    if (do_tser){
      cout << "# Averages over run: (<average> <rmsd> <error estimate>)\n" ;
      for(unsigned int i=0; i<props.size(); ++i){
	std::cout << "# " << props[i]->toTitle() << "\t" << props[i]->average() << "\n";
      }
    }
    
    // do we write distributions? (only scalar ones!)
    if (do_dist){
      std::cout.precision(6);
      std::cout.setf(std::ios::fixed, std::ios::floatfield);
      
      for(unsigned int i=0; i<props.size(); ++i){
	
	gmath::Stat<double> & stat = props[i]->getScalarStat();
        if (periodic) {
          stat.dist_init(dist_min, dist_max, dist_steps, true);
        }
        else if (dist_boundaries) {
	  stat.dist_init(dist_min, dist_max, dist_steps, false);
        } else {
	  stat.dist_init(dist_steps);
        }

	cout << "\n#" << endl;  
	cout << "# Distribution of      " << props[i]->toTitle() << endl;
	cout << "# values:              " << stat.n() << endl;
	cout << "# average value:       " << stat.ave() << endl;
	cout << "# rmsd:                " << stat.rmsd() << endl;
	cout << "# error estimate:      " << stat.ee() << endl;
        cout << "# maximum value at:    " << stat.distribution().maxValAt() << endl;

	if (normalize)
	  stat.distribution().write_normalized(cout);
	else
	  stat.distribution().write(cout);
      }
    }
  }
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}



