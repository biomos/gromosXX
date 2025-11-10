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
 * @file postcluster.cc
 * Processes the output from cluster
 */

/**
 * @page programs Program Documentation
 *
 * @anchor postcluster
 * @section postcluster Processes the output from cluster
 * @author @ref mc @ref co
 * @date 22-8-2006
 *
 * Program postcluster can do additional analyses on the output of @ref cluster
 * "cluster". Three different kinds of analyses are currently possible on
 * specified clusters:
 *
 * <ol>
 * <li> postcluster can perform a lifetime-analysis. A lifetime limit can be
 *      specified. This is the number of subsequent structures in the time
 *      series need to have switched to a different cluster before a true 
 *      transition to the new conformation is taken into account. This allows
 *      the user to disregard single events from being counted as a double
 *      transition to and from a new conformation. The program will write out
 *      the number of times a certain cluster is observed, its average lifetime.
 *      In addition it prints for every cluster the number of transitions to and
 *      from the other clusters.</li>
 * <li> postcluster can also be used to analyse combined clusterings, in which
 *      the original structures come from different sources (e.g. different
 *      trajectories). This can be used to assess the overlap in the sampled
 *      conformational space between two simulations. In honour of the infamous
 *      red_blue program, of Xavier Daura, this option is called rgb. By 
 *      specifying the number of frames from every individual source, the
 *      program will write out a file that can easily be used to produce a
 *      bar-plot in which the height of the bar indicates the size of the
 *      cluster and individual colors represent the portions of that cluster 
 *      coming from the different sources.</li>
 * <li> postcluster can be used to write out trajectory files and single
 *      structure files containing the central member structures of the 
 *      clusters. The trajectories can subsequently be used in any other 
 *      analysis program to monitor properties over all structures belonging to
 *      one cluster.</li>
 * </ol>
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@cluster_struct</td><td>&lt;structures file from cluster&gt; </td></tr>
 * <tr><td> \@cluster_ts</td><td>&lt;time-series file from cluster&gt; </td></tr>
 * <tr><td> \@clusters</td><td>&lt;@ref IntegerInputParser "structurespecifier"&gt; </td></tr>
 * <tr><td> [\@lifetime</td><td>&lt;lifetime limit&gt;] </td></tr>
 * <tr><td> [\@rgb</td><td>&lt;red&gt; &lt;green&gt; &lt;blue&gt; ...] </td></tr>
 * <tr><td> [\@traj</td><td>&lt;trajectory files&gt;] </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  postcluster
    @topo            ex.top
    @cluster_struct  cluster_structures.dat
    @cluster_ts      cluster_ts.dat
    @clusters        1-3
    @lifetime        5
    @rgb             200 200        
    @traj            ex.tr
 @endverbatim
 *
 * <hr>
 */


#include <cassert>
#include <cstdlib>
#include <map>
#include <string>
#include <vector>
#include <iomanip>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <set>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/fit/Reference.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gio/Ginstream.h"
#include "../src/gio/InTopology.h"
#include "../src/gio/OutG96.h"
#include "../src/gio/OutG96S.h"
#include "../src/bound/Boundary.h"
#include "../src/utils/IntegerInputParser.h"
#include "../src/gromos/Exception.h"

using namespace std;
using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;
using namespace utils;

/*
class IntegerInputParser: public set<int>
{
public:
  void addSpecifier(string const s, int maxnum);
protected:
  void parse(string const s, int maxnum);
};
*/
class cluster_parameter
{
public:
  int num;
  int maxstruct;
  int skip;
  int stride;
  double t0;
  double dt;
  double cutoff;
  bool free;
  int precision;
  int number_cluster;
  cluster_parameter() : num(0), maxstruct(-1), skip(0), stride(1), 
			t0(0.0), dt(1.0), cutoff(0.0),
			free(true), precision(10000), number_cluster(0) {};
};
 
void read_structure_file(string const s, cluster_parameter & cp,
		    vector< vector< int > > & cluster,
		    vector< int > & cm);
void read_ts_file(string const s, vector< int > &ts);
void split_trajectory(Arguments const &args, IntegerInputParser const & cs,
		      cluster_parameter const &cp, 
		      vector<int> const & centralMember, 
		      vector<int> const & timeSeries);
void determine_lifetime(IntegerInputParser const &cs, 
			cluster_parameter const &cp, 
		        vector<int> const &timeSeries, 
			int const lifeTimeLimit);

void colour_structures(cluster_parameter const & cp, 
		       vector<int> const & timeSeries, 
		       vector<int> const & rgb);


int main(int argc, char **argv){

  Argument_List knowns; 
  knowns << "topo" << "cluster_struct" << "cluster_ts" << "clusters"
         << "lifetime" << "traj" << "rgb" << "pbc";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo            <molecular topology file>\n";
  usage += "\t[@pbc            <pbc> [gathering method]]\n";
  usage += "\t@cluster_struct  <structures file from cluster>\n";
  usage += "\t@cluster_ts      <time-series file from cluster>\n";
  usage += "\t@clusters        <IntegerInputParser>\n";
  usage += "\t[@lifetime       <lifetime limit>]\n";
  usage += "\t[@rgb            <red> <green> <blue> ...]\n";
  usage += "\t[@traj           <trajectory files>]\n";
  
 
  try{
    Arguments args(argc, argv, knowns, usage);

    // read the lifeTimeLimit;
    int lifeTimeLimit = args.getValue<int>("lifetime", false, -1);

    // read the rgb values;
    vector<int> rgb;
    {
      Arguments::const_iterator iter=args.lower_bound("rgb"),
	to=args.upper_bound("rgb");
      int sum=0;
      for(; iter!=to ; ++iter){
	sum+=atoi(iter->second.c_str());
	rgb.push_back(sum);
      }
    } 
    
    // read in the cluster_structures file
    vector< vector< int > > cluster;
    vector< int > timeSeries;
    vector< int > centralMember;
    cluster_parameter cp;

    read_structure_file(args["cluster_struct"], cp, cluster, centralMember);
    read_ts_file(args["cluster_ts"], timeSeries);
    // cout << "ts size " << timeSeries.size() << endl;
    
    // structure specifier
    IntegerInputParser cs;
    {
      Arguments::const_iterator iter=args.lower_bound("clusters"),
	to=args.upper_bound("clusters");
      for(; iter!=to ; ++iter){
	cs.addSpecifier(iter->second, cp.number_cluster);
      }
    }
    
    if(rgb.size()){

      colour_structures(cp, timeSeries, rgb);
    }
    
    if(args.count("traj")>0){
  
      split_trajectory(args, cs, cp, centralMember, timeSeries);

    }

    if(lifeTimeLimit != -1){
      
      determine_lifetime(cs, cp, timeSeries, lifeTimeLimit);
	
    }
    

    

    
    
  }
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

/*
void IntegerInputParser::parse(string const s, int maxnum)
{
  if(s=="ALL" || s=="all"){
    for(int i=0; i<maxnum; i++) insert(i+1);
    return;
  }
  std::string::size_type iterator;
  if((iterator=s.find(',')) != std::string::npos){
    parse(s.substr(0,iterator), maxnum);
    parse(s.substr(iterator+1, std::string::npos), maxnum);
  }
  else{
    istringstream is;
    int rangeBegin, rangeEnd;
    if((iterator=s.find('-')) != std::string::npos){
      is.str(s.substr(0,iterator));
      if(!(is >> rangeBegin))
	throw gromos::Exception("IntegerInputParser", 
				"Invalid begin of range "+ s);
      is.clear();
      is.str(s.substr(iterator+1, std::string::npos));
      if(!(is >> rangeEnd))
	throw gromos::Exception("IntegerInputParser", 
				"Invalid end of range " + s);
      for(int i=rangeBegin; i<= rangeEnd; ++i){
	if(i> maxnum)
	  throw gromos::Exception("IntegerInputParser",
				  "Requested clusternumber too high: "+s);
	insert(i);
      }
    }
    else{
      is.str(s);
      if(!(is >> rangeBegin))
	throw gromos::Exception("IntegerInputParser", 
				"Invalid structure specified "+ s);
      if(rangeBegin > maxnum)
	throw gromos::Exception("IntegerInputParser",
				"Requested clusternumber too high: "+s);
      insert(rangeBegin);
    }
  }
}

void IntegerInputParser::addSpecifier(string const s, int maxnum)
{
  parse(s, maxnum);
}
*/
     
void read_structure_file(string const s, cluster_parameter & cp,
			 vector< vector< int > > & cluster,
			 vector< int > & cm)
{
  Ginstream gin(s);
  
  vector<string> buffer;
  gin.getblock(buffer);
  
  if(buffer[0] != "CLUSTERPARAMETERS")
    throw gromos::Exception("postcluster", "Could not read CLUSTERPARAMETERS "
			    "from file " + s );
  string b, e;
  istringstream is(gio::concatenate(buffer.begin()+1, buffer.end(), b));
  if(!(is >> cp.num >> cp.maxstruct >> cp.skip >> cp.stride
       >> cp.t0 >> cp.dt >> cp.cutoff >> cp.free >> cp.precision
       >> cp.number_cluster >> e))
    throw gromos::Exception("postcluster", 
			    "Error while reading CLUSTERPARAMETERS" );
  if(e.substr(0,3)!="END")
    throw gromos::Exception("postcluster",
			    "Error while readind CLUSTERPARAMETERS; no END");
  int p=cp.precision;
  cp.precision=1;
  for(int i=0; i< p; ++i){
    cp.precision *= 10;
  }

  
  gin.getblock(buffer);
  if(buffer[0] != "CLUSTER")
    throw gromos::Exception("postcluster", "Could not read CLUSTER "
			    "from file " + s );

  unsigned int numClusters=cp.number_cluster;
  if(cp.free) numClusters++;
  
  cm.resize(numClusters);
  string sdum;
  unsigned int inum;
  vector< int > csize(numClusters);
  
  for(unsigned int i=1; i< buffer.size()-1; ++i){
    is.clear();
    is.str(buffer[i]);
    
    if(!(is >> inum >> cm[i-1] >> sdum >> csize[i-1]))
      throw gromos::Exception("postcluster", "Error while reading CLUSTER "
			      "block. " + buffer[i]);
    if(inum != i-1)
      throw gromos::Exception("postcluster", "CLUSTER block is corrupted; "
			      "wrong cluster number " + buffer[i]);
    
  }
  if(buffer[buffer.size()-1].substr(0,3)!= "END")
    throw gromos::Exception("postcluster", "CLUSTER block is corrupted; "
			    "no END marker " + buffer[buffer.size()-1]);

  gin.getblock(buffer);
  if(buffer[0] != "MEMBERS")
    throw gromos::Exception("postcluster", "Could not read MEMBERS "
			    "from file " + s );
  

  is.clear();
  is.str(gio::concatenate(buffer.begin()+1, buffer.end(), b));
  cluster.resize(numClusters);

  gin.close();
  
  for(unsigned int i=0; i< numClusters; ++i){
    if(!(is >> inum) || i!=inum)
      throw gromos::Exception("postcluster", "Error while reading MEMBERS; "
			      "couldn't read cluster number");
    cluster[i].resize(csize[i]);
    
    for(int j=0; j < csize[i]; ++j){
      if(!(is >> cluster[i][j]))
	throw gromos::Exception("postcluster", "Error while reading MEMBERS; "
				"couldn't read a member");
    }
  }
  if(!(is >> sdum) || sdum.substr(0,3) != "END")
    throw gromos::Exception("postcluster",
			    "Error while reading MEMBERS; no END " + sdum);
  
}

  
void read_ts_file(string const s, vector< int > &ts)
{
  ifstream fin(s.c_str());
  string line;
  std::string::size_type iterator;
  int idum, iold = 0, cluster;
  double fdum;
  istringstream is;
  
  while(true){
    std::getline(fin, line);
    if(fin.eof()) break;
    
    if((iterator=line.find('#')) != std::string::npos){
      if(iterator==0) continue;
      line = line.substr(0,iterator);
    }
    
    is.clear();
    is.str(line);
    if(!(is >> fdum >> idum >> cluster))
      throw gromos::Exception("postcluster", "Error while reading time "
			      "series file: " + line);
    if(idum!=iold+1)
      throw gromos::Exception("postcluster", "Error while reading time "
			      "series file: structures not sequential");
    
    ts.push_back(cluster);
    iold++;
  }
}


void split_trajectory(Arguments const &args, IntegerInputParser const &cs,
		      cluster_parameter const &cp, 
		      vector<int> const &centralMember, 
		      vector<int> const &timeSeries)
{

  //  read topology
  InTopology it(args["topo"]);
  System sys(it.system());

    // Parse boundary conditions
    // Boundary *pbc = BoundaryParser::boundary(sys, args);
    // GatherParser
    // Boundary::MemPtr gathmethod = args::GatherParser::parse(sys,refSys,args);

  map<int, OutG96 *> file;
  vector<ofstream *> ff(cs.size());

  int count=0;
  for(IntegerInputParser::const_iterator iter = cs.begin(), 
	to=cs.end(); iter != to; ++iter, ++count){
    ostringstream os, ot;
    os << "cluster_" << *iter << ".trj";
    ot << "Structures belonging to cluster " << *iter << "\n"
       << "According to " << args["cluster_struct"] << "\n"
       << "         and " << args["cluster_ts"];
    ff[count] = new ofstream(os.str().c_str());
	
    file[*iter] = new OutG96(*ff[count]);
    file[*iter]->select("ALL");
    file[*iter]->writeTitle(ot.str());
    
    // cout << *iter << endl;
  }
  
  // define input stream
  int frameNumber=0;
  int frameCluster=0;
      
  for(Arguments::const_iterator iter=args.lower_bound("traj"), 
	to=args.upper_bound("traj"); iter != to; ++iter){
    InG96 ic;
    ic.open(iter->second);
    //ic.select("ALL");
	
    while(!ic.eof()){
      ic.select("ALL");
      ic >> sys;
	  
      if(!((frameNumber - cp.skip) % cp.stride)){
	    
	if(cs.count(timeSeries[frameCluster])){
	      
	  (*file[ timeSeries[ frameCluster ] ]) << sys;
	}  

	// check whether it is a central member
        IntegerInputParser::const_iterator iter2 = cs.find(timeSeries[frameCluster]);
	if( (centralMember[timeSeries[frameCluster]]== frameCluster+1)
            && iter2 != cs.end()) {
	  ostringstream os, ot;
	  os << "cluster_" << timeSeries[frameCluster] << ".cms";
	  ot << "Central member structure belonging to cluster " 
	     << timeSeries[frameCluster] << "\n"
             << "Structure " << frameCluster+1
             << " at time " << (frameCluster) * cp.dt + cp.t0 << "\n"
	     << "According to " << args["cluster_struct"] << "\n"
	     << "         and " << args["cluster_ts"];
	  ofstream yaf(os.str().c_str());
	  OutG96S og96(yaf);
	  og96.select("ALL");
	  og96.writeTitle(ot.str());
	  og96 << sys;
	  og96.close();
	  yaf.close();
	}
	    
	    
	frameCluster++;
      }
	  
      frameNumber++;
	  
	  
    }
    
    ic.close();

  }
      
  // close the files
  for(map<int, OutG96 *>::iterator iter=file.begin(),
	to=file.end(); iter!=to; ++iter){
    iter->second->close();
    delete iter->second;
  }
  for(unsigned int i=0; i< ff.size(); ++i){
    ff[i]->close();
    delete ff[i];
  }
}



void determine_lifetime(IntegerInputParser const &cs, 
			cluster_parameter const &cp, 
		        vector<int> const &timeSeries,
			int const lifeTimeLimit)
{
  // the life time means that you
  // 1. you only consider a cluster visited after so many hits
  // 2. you consider less than that misses, as not happening.
  //    if after such a short break, the old one continues for
  //    at least the life time again, you treat that one as continuing.

  int current = -1;
  int numCluster=cp.number_cluster;
  if(cp.free) numCluster++;
  
  vector<int> lifetime(numCluster, 0);
  vector<int> occurence(numCluster, 0);
  vector< vector < unsigned short > > pathways(numCluster);
  for(int i=0; i < numCluster; i++)
    pathways[i].resize(numCluster,0);
  
  bool ok=false;
  //cout << "cp.number_cluster " << cp.number_cluster << endl;
  // check the pathways in the first lifeTimeLimit - 1 structures
  int previous=timeSeries[0];
  
  for(int i=1; i< lifeTimeLimit-1; ++i){
    ++pathways[timeSeries[i]][previous];
    previous=timeSeries[i];
  }
  
  
  for(unsigned int i=lifeTimeLimit-1; i< timeSeries.size(); ++i){
    //cout << lifetime[0] << " " << lifetime[1] << " " << lifetime[2] << " " << lifetime[3] << endl;
     ++pathways[timeSeries[i]][previous];
     previous = timeSeries[i];
     
    if(current==-1){
      ok=true;
	 
      for(int j=1; j < lifeTimeLimit; ++j){
	if(timeSeries[i]!=timeSeries[i-j]){
	  ok=false;
	  break;
	}
      }
      if(ok) {
	current=timeSeries[i];
	lifetime[current]+=lifeTimeLimit;
      }
    }
    else{
      if(timeSeries[i-1] == timeSeries[i]){
	++lifetime[current];
      }
      else{
	// maybe it is just a break
	ok=false;
	for(int k=1; k<lifeTimeLimit; ++k){
	  ok=true;
	  for(int l=0; l <lifeTimeLimit; ++l){
	    if((i+k+l)>=timeSeries.size() ||
	       current != timeSeries[i+k+l]) {
	      ok=false;
	      break;
	    }
	  }
	  if(ok==true){
	    break;
	  }
	}
	if(ok==false){
	    
	  ++occurence[current];
	  if(lifeTimeLimit!=1)
	    current=-1;
	  else{
	    current=timeSeries[i];
	    ++lifetime[current];
	  }
	}
      }
    }
  }
    
  if(current != -1)	  
    ++occurence[current];

  // write the output
  cout << "#  cluster   visited  <lifetime>\n";
  for(IntegerInputParser::const_iterator iter=cs.begin(), to=cs.end();
      iter!=to; ++iter){
    cout << setw(10) << *iter
	 << setw(10) << occurence[*iter];
    if(occurence[*iter])
      cout << setw(12) << lifetime[*iter] * cp.dt / occurence[*iter];
    else
      cout << setw(12) << 0;
    
    cout << endl;
  }
  int sum=0;
  for(unsigned int j=0; j<lifetime.size(); ++j) sum+=lifetime[j];
  
  cout << "#\n# Number of asocial structures " << cp.num - sum -1
       << " (" << (cp.num - sum - 1) * cp.dt << ")\n"
       << "# (i.e. not being part of any visit)\n#\n";

  for(IntegerInputParser::const_iterator iter=cs.begin(), to=cs.end();
      iter!=to; ++iter){
    cout << "\n-- cluster " << *iter 
	 <<" -------------------------------------------------------\n";
    cout << "coming from: \n\t";
    int count=0;
    
    for(int i=0; i< numCluster; i++){
      if(i!=*iter && pathways[*iter][i] != 0){
	if(count && count % 5 ==0) cout << "\n\t";
	
	ostringstream os;
	os << i << " (" << pathways[*iter][i] << ")";
	cout << setw(10) << os.str();
	++count;
      }
    }
    cout << "\ngoing to: \n\t";
    count=0;
    
    for(int i=0; i< numCluster; i++){
      if(i!=*iter && pathways[i][*iter] != 0){
	if(count && count % 5 ==0) cout << "\n\t";
	
	ostringstream os;
	os << i << " (" << pathways[i][*iter] << ")";
	cout << setw(10) << os.str();
	++count;
      }
    }
    cout << "\n";
    
  }
  
}

void colour_structures(cluster_parameter const & cp, 
		       vector<int> const & timeSeries, 
		       vector<int> const & rgb)
{
  
      
      int numCluster=cp.number_cluster;
      if(cp.free) ++numCluster;
      
      if(rgb.back() != cp.num - 1)
	throw gromos::Exception("postcluster", "rgb input is wrong\n"
				"sum of entries is not equal to the total "
				"number of structures\n");
      vector< vector < int > > rgb_size(numCluster);
      for(int i=0; i< numCluster; ++i)
	rgb_size[i].resize(rgb.size(),0);
      
      // loop over the timeSeries
      int colourCounter=0;
      int tsSize=timeSeries.size();
      
      for(int i=0; i< tsSize; i++){
	if(i >= rgb[colourCounter]) ++colourCounter;
	++rgb_size[timeSeries[i]][colourCounter];
      }
      
      ofstream fout("cluster_rgb.dat");
      int start=0;
      if(cp.free) start=1;
      fout << "# cluster" << setw(rgb.size()*10-1) << "cumulative"
	   << setw(rgb.size()*10+2) << "absolute";
      fout << "\n#       ";
      for(unsigned j=0; j < rgb.size(); ++j){
	fout << setw(10) << rgb.size()-j;
      }
      fout << "  ";
      
      for(unsigned j=0; j < rgb.size(); ++j){
	fout << setw(10) << j+1;
      }
      fout << endl;
      
      
      for(int i=start; i<numCluster; ++i){
	int sum=0;
	for(unsigned j=0; j < rgb.size(); ++j){
	  sum += rgb_size[i][j];
	}
	fout << setw(8) << i;
	for(unsigned j=0; j < rgb.size(); ++j){
	  fout << setw(10) << sum;
	  sum -= rgb_size[i][rgb.size()-1-j];
	}
	fout << "  ";
	for(unsigned j=0; j < rgb.size(); ++j){
	  fout << setw(10) << rgb_size[i][j];
	}
	fout << endl;
      }
}
