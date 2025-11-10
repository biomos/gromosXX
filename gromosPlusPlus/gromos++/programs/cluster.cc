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
 * @file cluster.cc
 * Performs a conformational clustering on a RMSD matrix
 */

/**
 * @page programs Program Documentation
 *
 * @anchor cluster
 * @section cluster Performs a conformational clustering on a RMSD matrix
 * @author @ref mc @ref co
 * @date 22-8-06
 *
 * Program cluster performs a conformational clustering based on a similarity
 * matrix, such as calculated by the program @ref rmsdmat "rmsdmat". The
 * clustering algorithm is the one described in [Proteins 1999, 34, 269 - 280].
 * Structures with rmsd values smaller than a user specified cutoff are
 * considered to be structural neighbours. The structure with the highest
 * number of neighbours is considered to be the central member of the cluster
 * of similar structures forming a conformation. After removing all structures
 * belonging to this first cluster, the procedure is repeated to find the
 * second, third etc. most populated clusters.
 * 
 * One specific structure can be forced to be the central member structure of
 * the first cluster, this can also be the reference structure, by specifying
 * structure number 0. The clustering can be performed on a subset of the
 * matrix, by specifying the maximum number of structures to consider. This
 * allows for an assessment of the development of the number of clusters over
 * time.
 *
 * Depending on the settings used for program @ref rmsdmat, the flag human
 * may need to be specified to ensure proper reading in of the matrix.
 *
 * Clusters may be further analysed using program @ref postcluster 
 * "postcluster".
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@rmsdmat</td><td>&lt;rmsd matrix file name&gt; </td></tr>
 * <tr><td> \@cutoff</td><td>&lt;cutoff&gt; </td></tr>
 * <tr><td> \@time</td><td>&lt; @ref utils::Time "time and dt"&gt; </td></tr>
 * <tr><td> \@precision</td><td>&lt;number of digits in the matrix &gt; </td></tr>
 * <tr><td> [\@maxstruct</td><td>&lt;maximum number of structures to consider&gt;] </td></tr>
 * <tr><td> [\@human</td><td>(use a human readable matrix)] </td></tr>
 * <tr><td> [\@force</td><td>&lt;structure&gt; (force clustering on the indicated structure, 0 is the reference)] </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  cluster   
    @rmsdmat    RMSDMAT.bin
    @cutoff     0.12
    @time       0 1
    @maxstruct  100
    @precision  4
    @human
    @force      23
 @endverbatim
 *
 * <hr>
 */

#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <cmath>

#include "../src/args/Arguments.h"
#include "../src/gio/Ginstream.h"
#include "../src/gio/gzstream.h"
#include "../src/gromos/Exception.h"

using namespace std;
using namespace args;

// change this typedef if you have more than 65536 structures
// typedef unsigned short the_type;
// typedef unsigned int the_type;

class cluster_parameter {
public:
  int num;
  int maxstruct;
  int skip;
  int stride;
  double t0;
  double dt;
  double cutoff;
  bool free;
  int force_ref;
  int precision;
  int number_cluster;

  cluster_parameter() : num(0), maxstruct(-1), skip(0), stride(1),
  t0(0.0), dt(1.0), cutoff(0.0),
  free(true), force_ref(0), precision(10000), number_cluster(0) {
  };
};

template<typename the_type>
void read_matrix(string const filename, vector< vector < the_type > > &matrix,
        bool const human, cluster_parameter & cp, int precision);

template<typename the_type>
int cluster_analysis(Arguments & args);

int main(int argc, char **argv) {

  Argument_List knowns;
  knowns << "rmsdmat" << "cutoff" << "human" << "force" << "maxstruct" << "time" << "precision";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@rmsdmat    <rmsd matrix file name>\n";
  usage += "\t@cutoff     <cutoff>\n";
  usage += "\t@time       <t0> <dt>\n";
  usage += "\t@precision   <number of digits in the matrix (default 4)>\n";
  usage += "\t[@maxstruct <maximum number of structures to consider>]\n";
  usage += "\t[@human     (use a human readable matrix)]\n";
  usage += "\t[@force     <structure> (force clustering on the indicated structure, 0 is the reference)]\n";

  try {
    Arguments args(argc, argv, knowns, usage);
    typedef unsigned int uint;
    typedef unsigned short ushort;

    if (args.getValue<int>("precision", true) > 4) {
      return cluster_analysis<uint > (args);
    } else {
      return cluster_analysis<ushort > (args);
    }
  }  catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

template<typename the_type>
int cluster_analysis(Arguments & args) {
  // create the cluster parameters
  cluster_parameter cp;

  // get the cutoff
  cp.cutoff = args.getValue<double>("cutoff", true);
  if (cp.cutoff <= 0.0) {
    throw gromos::Exception("cluster", "cutoff should be > 0.0");
  }
  
  // read the precision
  int ii = args.getValue<int>("precision", true);
  int precision = 1;
  for (int i = 0; i < ii; ++i) {
    precision *= 10;
  }

  // get the time
  vector<double> timearg = args.getValues<double>("time", 2, false,
          Arguments::Default<double>() << 0.0 << 1.0);
  cp.t0 = timearg[0];
  cp.dt = timearg[1];
  double time = cp.t0;

  // try to read the maxstruct
  cp.maxstruct = args.getValue<int>("maxstruct", false, -1);

  // is the matrix human readable
  bool human = false;
  if (args.count("human") >= 0) human = true;

  // is the clustering free
  if (args.count("force") >= 0) cp.free = false;
  if (args.count("force") > 0) {
    std::istringstream is(args["force"]);
    if (!(is >> cp.force_ref))
      throw gromos::Exception("cluster", "could not read structure number for forced clustering");
  }

  // create the data structure
  vector< vector <the_type> > pairs;

  // read matrix
  read_matrix(args["rmsdmat"], pairs, human, cp, precision);

  // now we are almost done
  size_t num = pairs.size();
  vector<int> taken(num, -1);
  vector<the_type> central_member;
  vector<vector <the_type> > cluster;
  int remaining = num;
  int clustercount = 0;

  // set the first cluster
  central_member.push_back(cp.force_ref);
  cluster.push_back(pairs[cp.force_ref]);

  if (cp.free)
    pairs[0].clear();
  taken[0] = 0;

  if (!cp.free) {
    // mark them as taken
    for (size_t j = 0, jto = pairs[cp.force_ref].size(); j != jto; ++j) {
      taken[pairs[cp.force_ref][j]] = 0;
      pairs[pairs[cp.force_ref][j]].clear();
    }
    // Markus: BUG??? changed order
    // take them out 
    remaining -= pairs[cp.force_ref].size();
    pairs[cp.force_ref].clear();

    for (size_t i = 1; i < num; ++i) {
      vector< the_type > temp;
      for (size_t j = 0, jto = pairs[i].size(); j != jto; ++j) {
        if (taken[pairs[i][j]] == -1) {
          temp.push_back(pairs[i][j]);
        }
      }
      pairs[i] = temp;
    }
  }

  //while(remaining){
  while (true) {

    // search for the one with the largest number of neighbours
    size_t maxsize = 0;
    size_t maxindex = 0;
    for (size_t i = 0; i < num; ++i) {
      if (pairs[i].size() > maxsize) {
        maxsize = pairs[i].size();
        maxindex = i;
      }
    }
    if (!maxsize) break;

    // put them in
    clustercount++;
    central_member.push_back(maxindex);
    cluster.push_back(pairs[maxindex]);

    // and take them out
    remaining -= pairs[maxindex].size();

    for (size_t j = 0, jto = pairs[maxindex].size(); j != jto; ++j) {
      taken[pairs[maxindex][j]] = clustercount;
      pairs[pairs[maxindex][j]].clear();
    }
    pairs[maxindex].clear();
    for (size_t i = 0; i < num; ++i) {
      vector< the_type > temp;
      for (size_t j = 0, jto = pairs[i].size(); j != jto; ++j) {
        if (taken[pairs[i][j]] == -1) {
          temp.push_back(pairs[i][j]);
        }
      }
      pairs[i] = temp;
    }
  } // while remaining

  // NOW PRODUCE OUTPUT
  // time series
  {
    ofstream fout("cluster_ts.dat");
    for (size_t i = 1; i < num; ++i) {
      fout << setw(10) << time
              << setw(10) << i
              << setw(10) << taken[i] << endl;
      time += cp.dt;
    }
  }
  // cluster contents
  cp.number_cluster = cluster.size();
  if (cp.free) cp.number_cluster--;

  {
    ofstream fout("cluster_structures.dat");
    fout << "TITLE\n"
            << "\tClustering from " << args["rmsdmat"] << "\n"
            << "\tCutoff : " << cp.cutoff << "\n"
            << "\tTotal number of structures : " << num << "\n";
    if (cp.free) fout << "\tFree clustering performed\n";
    else fout << "\tFirst cluster forced to contain structure " << cp.force_ref << "\n";
    fout << "\tTotal number of clusters found : " << cp.number_cluster
            << "\nEND\n";
    fout << "CLUSTERPARAMETERS\n"
            << "#  structs   maxstruct        skip      stride\n"
            << setw(10) << cp.num
            << setw(12) << cp.maxstruct
            << setw(12) << cp.skip
            << setw(12) << cp.stride << "\n"
            << "#       t0          dt\n"
            << setw(10) << cp.t0
            << setw(12) << cp.dt << "\n"
            << "#   cutoff        free   precision\n"
            << setw(10) << cp.cutoff
            << setw(12) << ((cp.free) ? 1 : 0)
            << setw(12) << rint(log(double(cp.precision)) / log(10.0)) << "\n"
            << "# clusters\n"
            << setw(10) << cp.number_cluster << "\n"
            << "END\n";

    fout << "CLUSTER\n";
    fout << "#    clu  center    time    size\n";
    for (size_t i = 0, ito = cluster.size(); i != ito; ++i) {
      fout << setw(8) << i
              << setw(8) << central_member[i];
      if (central_member[i] == 0)
        fout << setw(8) << "ref";
      else
        fout << setw(8) << cp.t0 + (central_member[i] - 1) * cp.dt;
      fout << setw(8) << cluster[i].size() << endl;
    }
    fout << "END\n";
    fout << "MEMBERS\n";
    for (size_t i = 0, ito = cluster.size(); i != ito; ++i) {
      fout << setw(6) << i;
      for (size_t j = 0, jto = cluster[i].size(); j != jto; ++j) {
        if (j % 10 == 0 && j != 0) fout << "\n      ";
        fout << " " << setw(7)  << cluster[i][j];
      }
      fout << "\n";
    }
    if ((cluster[cluster.size() - 1].size()) % 10 == 0) fout << "\n";
    fout << "END\n";
  }

  // standard output
  {
    cout << "# Clustering from " << args["rmsdmat"] << "\n"
            << "# Cutoff: " << cp.cutoff << "\n"
            << "# Total number of structures: " << num << "\n";
    if (cp.free) cout << "# Free clustering performed\n";
    else cout << "# First cluster forced to contain reference structure\n";
    cout << "# Total number of clusters found: " << cp.number_cluster << "\n";

    cout << "#\n";
    cout << "#    clu    size\n";
    if (cp.free) cout << setw(8) << "ref";
    else cout << setw(8) << 0;
    cout << setw(8) << cluster[0].size() << endl;

    for (size_t i = 1, ito = cluster.size(); i != ito; ++i) {
      cout << setw(8) << i
              << setw(8) << cluster[i].size() << endl;
    }
  }
  return 0;
}

template<typename the_type>
void read_matrix(string const filename, vector< vector < the_type > > &matrix,
        bool const human, cluster_parameter & cp, int precision) {
  
  if (human) {
    gio::Ginstream gin(filename);
    if (!gin.stream()) {
      throw gromos::Exception("cluster", "Error opening rmsdmat file\n");
    }

    string sdum;

    gin.getline(sdum);
    if (sdum != "RMSDMAT")
      throw gromos::Exception("cluster", "Error while reading rmsdmat file\n"
            "expected RMSDMAT block, got " + sdum);
    gin.getline(sdum);
    istringstream is(sdum);
    if (!(is >> cp.num >> cp.skip >> cp.stride))
      throw gromos::Exception("cluster", "Error while reading rmsdmat file\n"
            "could not read number of structures\n");
    if (cp.num < 0)
      throw gromos::Exception("cluster", "Error while reading rmsdmat file\n"
            "read negative number structures\n");
    if (cp.num > pow(2.0, 8.0 * sizeof (the_type))) {
      std::cerr << "number of structures: " << cp.num << "\n"
              << "maximum number possible: " << pow(2.0, 8.0 * sizeof (the_type))
              << std::endl;
      throw gromos::Exception("cluster", "GROMOS96 ERROR: number of "
              "structures is too large for data type. "
              "Change typedef in program and recompile");
    }

    if (cp.skip < 0)
      throw gromos::Exception("cluster", "Error while reading rmsdmat file\n"
            "read negative number of skipped structures\n");
    if (cp.stride < 0)
      throw gromos::Exception("cluster", "Error while reading rmsdmat file\n"
            "read negative number of structure stride\n");
    gin.getline(sdum);
    is.clear();
    is.str(sdum);
    if (!(is >> cp.precision))
      throw gromos::Exception("cluster", "Error while reading rmsdmat file\n"
            "could not read precision\n");
    if (cp.precision < 0)
      throw gromos::Exception("cluster", "Error while reading rmsdmat file\n"
            "read negative precision\n");
    if (cp.precision != precision)
      throw gromos::Exception("cluster", "Error while reading rmsdmat file\n"
            "Matrix has different precision as given in @precision!\n");

    double cuttest = cp.cutoff * cp.precision / 10;
    if (fabs(double(int(cuttest)) - cuttest) > cp.cutoff / 100.0)
      throw gromos::Exception("cluster", "A cutoff with this precision "
            "requires a higher precision in the rmsd "
            "matrix. \nYes that means that you have to "
            "redo your matrix $%^$#$@$%!!");

    if (cp.maxstruct < 0) cp.maxstruct = cp.num;
    else cp.maxstruct++;

    int icutoff = int(cp.cutoff * cp.precision);

    matrix.resize(cp.maxstruct);

    int ii, jj, rmsd;

    for (int i = 0; i < cp.num; ++i) {
      if (i < cp.maxstruct) {
        //cout << "diagonal " << i << endl;

        matrix[i].push_back(i);
      }

      for (int j = i + 1; j < cp.num; ++j) {
        gin.getline(sdum);
        is.clear();
        is.str(sdum);
        if (!(is >> ii >> jj >> rmsd))
          throw gromos::Exception("cluster", "Error while reading rmsdmat file\n"
                "could not read line:" + is.str());
        if (ii != i || jj != j) {
          ostringstream os;
          os << "Now we are really confused!\n"
                  << i << " = " << ii << " or " << j << " = " << jj << "\n";
          throw gromos::Exception("cluster", os.str());
        }

        if (ii < cp.maxstruct && jj < cp.maxstruct && rmsd < icutoff) {
          matrix[ii].push_back(jj);
          if (ii > 0)
            matrix[jj].push_back(ii);
        }
      }
    }
    if (!(gin.getline(sdum)))
      throw gromos::Exception("cluster", "Error while reading rmsdmat file\n"
            "where is the END?");
    if (sdum.substr(0, 3) != "END")
      throw gromos::Exception("cluster", "Error while reading rmsdmat file\n"
            "Trailing data on file. Expected 'END', got"
            + sdum);
  } else {
    igzstream fin(filename.c_str(), ios::in | ios::binary);
    if (!fin) {
      throw gromos::Exception("cluster", "Error opening rmsdmat file\n");
    }

    if (!fin.read((char*) &cp.num, sizeof (int)))
      throw gromos::Exception("cluster", "Error while reading rmsdmat file\n"
            "could not read number of structures\n");
    if (cp.num < 0)
      throw gromos::Exception("cluster", "Error while reading rmsdmat file\n"
            "read negative number of structures\n");
    if (cp.num > pow(2.0, 8.0 * sizeof (the_type))) {
      std::cerr << "number of structures: " << cp.num << "\n"
              << "maximum number possible: " << pow(2.0, 8.0 * sizeof (the_type))
              << std::endl;
      throw gromos::Exception("cluster", "GROMOS96 ERROR: number of "
              "structures is too large for data type. "
              "Change typedef in program and recompile");
    }
    if (!fin.read((char*) &cp.skip, sizeof (int)))
      throw gromos::Exception("cluster", "Error while reading rmsdmat file\n"
            "could not read skip\n");
    if (cp.skip < 0)
      throw gromos::Exception("cluster", "Error while reading rmsdmat file\n"
            "read number of skipped structures\n");
    if (!fin.read((char*) &cp.stride, sizeof (int)))
      throw gromos::Exception("cluster", "Error while reading rmsdmat file\n"
            "could not read stride\n");
    if (cp.stride < 0)
      throw gromos::Exception("cluster", "Error while reading rmsdmat file\n"
            "read negative number of stride structures\n");

    if (!fin.read((char*) &cp.precision, sizeof (int)))
      throw gromos::Exception("cluster", "Error while reading rmsdmat file\n"
            "could not read precision\n");
    if (cp.precision < 0)
      throw gromos::Exception("cluster", "Error while reading rmsdmat file\n"
            "read negative precision\n");
    if (cp.precision != precision)
      throw gromos::Exception("cluster", "Error while reading rmsdmat file\n"
            "Matrix has different precision as given in @precision!\n");

    double cuttest = cp.cutoff * cp.precision / 10;
    if (fabs(double(rint(cuttest)) - cuttest) > cp.cutoff / 100.0)
      throw gromos::Exception("cluster", "A cutoff with this precision "
            "requires a higher precision in the rmsd "
            "matrix. \nYes that means that you have to "
            "redo your matrix $%^$#$@$%!!");
    if (cp.maxstruct < 0) cp.maxstruct = cp.num;
    else cp.maxstruct++;

    unsigned int icutoff = unsigned(cp.cutoff * cp.precision);

    matrix.resize(cp.maxstruct);

    if (cp.precision < 1e5) {
      typedef unsigned short ushort;

      ushort rmsd;

      for (int i = 0; i < cp.num; ++i) {
        if (i < cp.maxstruct)
          matrix[i].push_back(i);

        for (int j = i + 1; j < cp.num; ++j) {
          if (!fin.read((char*) &rmsd, sizeof (ushort)))
            throw gromos::Exception("cluster", "Error while reading rmsdmat file\n"
                  "file corrupt");
          if (i < cp.maxstruct && j < cp.maxstruct && rmsd < short(icutoff)) {
            matrix[i].push_back(j);
            if (i > 0)
              matrix[j].push_back(i);
          }
        }
      }
    } else {
      unsigned rmsd;

      for (int i = 0; i < cp.num; ++i) {
        if (i < cp.maxstruct)
          matrix[i].push_back(i);

        for (int j = i + 1; j < cp.num; ++j) {
          if (!fin.read((char*) &rmsd, sizeof (unsigned)))
            throw gromos::Exception("cluster", "Error while reading rmsdmat file\n"
                  "file corrupt");
          if (i < cp.maxstruct && j < cp.maxstruct && rmsd < icutoff) {
            matrix[i].push_back(j);
            if (i > 0)
              matrix[j].push_back(i);
          }
        }
      }
    } // if precision
  } // if human
}
