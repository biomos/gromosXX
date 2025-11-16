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
 * @file dfgrid.cc
 * Calculate a distancefield grid for given coordinates and write a combined coordinate file.
 */

/**
 * @page programs Program Documentation
 *
 * @anchor dfgrid
 * @section dfgrid calculate distancefield grid and write to coordinate file
 * @author @ref mp
 * @date 11-11-2016
 *
 * Program dfgrid calculates a distancefield grid analogously to the way it is
 * done in md++, where distancefield (DF) distances can be used as restraints and
 * reaction coordinates in path pulling methods (de Ruiter \& Oostenbrink. 
 * J. Chem. Theor. Comput., 9:883 – 892, 2012). The DF method is based on the 
 * mapping of distances from a target point on a grid, where grid points that 
 * overlap with the protein are penalized. As a result the shortest path for a 
 * ligand to the target point will never go through the protein, making it less
 * likely to get stuck in a dead end.
 *
 * Program dfgrid calculates the grid for any given snapshot, facilitating the 
 * visualization and analysis. It writes out a coordinate file which contains 
 * both the input coordinates and the distancefield grid, in addition to the 
 * target (zero-distance) point which will often be a virtual atom. When choosing 
 * pdb as output format, the distances on the grid will be written to the b-factor 
 * column for easy visualization.
 *
 * The use of either \@stride or \@frames to reduce the analysis to a few 
 * snapshots is recommended.
 *
 * The following flags are defined in the same way as in the md++ parameter and
 * distance restraints files: \@gridspacing, \@proteinoffset, \@proteincutoff 
 *  and \@smooth , see Gromos Manual Vol. 2-9.12. \@proteinatoms has the same
 * function as the correponding md++ parameter, but here the atom selection is
 * specified in the form of an atomspecifier.
 *
 * \@max allows to specify a maximum distance, grid points to which higher 
 * distances are mapped will not be written to the coordinate file. \@protect
 * protects grid points within a certain radius from the target point from 
 * being flagged as protein.
 * 
 * With the \@distatoms flag you can specify (virtual) atoms for which the 
 * df distance will be printed in standard output and for which the shortest
 * df path will be added to the output coordinates for visual inspection. If 
 * one is only interested in the distances, \@nogrid will prevent writing of
 * the coordinate file.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0> 
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc </td><td>&lt;boundary type&gt; [<gather method&gt;] </td></tr>
 * <tr><td> \@atom</td><td>&lt;(virtual) atom specifier for the target (zero-distance) point&gt; </td></tr>
 * <tr><td> [\@distatoms</td><td>&lt;(virtual) atom specifier for atoms for which to output the df distance&gt;] </td></tr>
 * <tr><td> \@gridspacing</td><td>&lt;grid spacing&gt; </td></tr>
 * <tr><td> \@proteinoffset</td><td>&lt;penalty (nm) for being in the protein&gt; </td></tr>
 * <tr><td> \@proteincutoff</td><td>&lt;cutoff to determine gridpoints within the protein&gt; </td></tr>
 * <tr><td> \@proteinatoms</td><td>&lt;atomspecifier for everything considered as protein&gt; </td></tr>
 * <tr><td> [\@smooth</td><td>&lt;number of rounds to smoothen the forces at the edge of the protein&gt;] </td></tr>
 * <tr><td> [\@max</td><td>&lt;maximum distance: do not write out grid points with higher distances&gt; (default: 1)]</td></tr>
 * <tr><td> [\@protect</td><td>&lt;radius around the target atom that will not be flagged as protein&gt;] </td></tr>
 * <tr><td> [\@outformat</td><td>&lt;output coordinates format&gt;] </td></tr>
 * <tr><td> [\@notimeblock</td><td>&lt;do not write timestep block&gt;] </td></tr>
 * <tr><td> [\@nogrid</td><td>&lt;do not write grid coordinate file&gt;] </td></tr> 
 * <tr><td> [\@time </td><td>&lt;time and dt&gt;] </td></tr> 
 * <tr><td> [\@stride</td><td>&lt;write every nth frame&gt; (default: 1)]</td></tr> 
 * <tr><td> [\@frames</td><td>&lt;select frames to write out, starts at 0&gt; (default: 0)]</td></tr>
 * <tr><td> \@traj</td><td>&lt;input trajectory files&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  dfgrid
    @topo      ex.top
    @pbc       r cog
    @atom      va(-1,1:1,4,5)
    @gridspacing 0.2
    @proteinoffset 15
    @proteincutoff 0.2
    @proteinatoms 5311
    @max 16
    @smooth 1
    @protect 0.4
    [@notimeblock]
    [@time     0 2]
    [@stride  3]
    [@frames 10 20 40]
    @outformat pdb
    @traj      ex.tr
 @endverbatim
 *
 * <hr>
 */

#include <cassert>
#include <cstdlib>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <set>
#include <string>
#include <cmath>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/args/OutformatParser.h"
#include "../src/gio/InG96.h"
#include "../src/gio/OutG96.h"
#include "../src/bound/Boundary.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/Bond.h"
#include "../src/gcore/Solvent.h"
#include "../src/gcore/SolventTopology.h"
#include "../src/gcore/Box.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/utils/SimplePairlist.h"
#include "../src/gio/InTopology.h"
#include "../src/gmath/Vec.h"
#include "../src/utils/groTime.h"
#include "../src/gromos/Exception.h"

using namespace gcore;
using namespace gio;
using namespace args;
using namespace bound;
using namespace std;

int neighbour(int i, int j, std::vector<int> &ngrid){
  int n;
  switch(j){
    case 0:
      // backward in x
      if(i % ngrid[0] == 0) return i -1 + ngrid[0];
      else return i-1;
    case 1:
      // forward in x
      if((i+1) % ngrid[0] ==0) return i + 1 - ngrid[0];
      else return i+1;
    case 2:
      // backward in y
      n = i - (i % ngrid[0]);
      if(n % (ngrid[0]*ngrid[1]) == 0) return i - ngrid[0] + (ngrid[0]*ngrid[1]);
      else return i - ngrid[0];
    case 3:
      // forward in y
      n = i + ngrid[0] - (i% ngrid[0]);
      if(n % (ngrid[0]*ngrid[1]) == 0) return i + ngrid[0] - (ngrid[0]*ngrid[1]);
      else return i + ngrid[0];
    case 4:
      // backward in z
      n = i - (i % (ngrid[0]*ngrid[1]));
      if(n==0) return i - (ngrid[0]*ngrid[1]) + ngrid[0]*ngrid[1]*ngrid[2];
      else return i - ngrid[0] * ngrid[1];
    case 5:
      // forward in z
      n = i + ngrid[0]*ngrid[1] - (i % (ngrid[0]*ngrid[1]));
      if(n==(ngrid[0]*ngrid[1]*ngrid[2])) return i + ngrid[0]*ngrid[1] - (ngrid[0]*ngrid[1]*ngrid[2]);
      else return i + ngrid[0]*ngrid[1];
    default:
      return 0;
  }
  return 0;
}

// get coordinates from gridpoint index
gmath::Vec get_coords(int idx, std::vector<int> &ngrid, double gridspacing, double boxK,double boxL,double boxM) {
  int nz = int(double(idx) / double(ngrid[0] * ngrid[1]));
  int ny = int(double(idx % (ngrid[0] * ngrid[1]) ) / double(ngrid[0]));
  int nx = (idx % (ngrid[0] * ngrid[1])) % ngrid[0];
  gmath::Vec pos(nx*gridspacing - boxK/2, ny*gridspacing - boxL/2,nz*gridspacing - boxM/2);
  return pos;
}


int main(int argc, char **argv) {

  Argument_List knowns;
  knowns << "topo" << "pbc" << "traj" << "atom" << "proteincutoff" << "proteinoffset"
         << "proteinatoms" << "gridspacing" << "smooth" << "outformat" << "stride" 
         << "notimeblock" << "time" << "max" << "protect" << "frames" << "distatoms"
         << "nogrid";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo     <molecular topology file>\n";
  usage += "\t@pbc      <boundary type> [<gather method>]\n";
  usage += "\t@atom     <(virtual) atom specifier for the target (zero-distance) point>\n";
  usage += "\t[@distatoms     <(virtual) atom specifier for atoms for which to output the df distance>]\n";
  usage += "\t@gridspacing   <grid spacing>\n";
  usage += "\t@proteinoffset <penalty for being in the protein>\n";
  usage += "\t@proteincutoff <cutoff to determine gridpoints within the protein>\n";
  usage += "\t@proteinatoms  <atomspecifier for everything to be penalized as protein> \n";
  usage += "\t[@max     <maximum distance: do not write out grid points with higher distances> (default: 1)]\n";
  usage += "\t[@nogrid  <do not write out grid coordinate file>]\n";
  usage += "\t[@smooth  <number of rounds to smoothen the forces at the edge of the protein>]\n";
  usage += "\t[@protect <radius around the target atom that will not be flagged as protein>]\n";
  usage += "\t[@outformat   <output coordinates format>]\n";
  usage += "\t[@notimeblock <do not write timestep block>]\n";
  usage += "\t[@time    <time and dt>]\n"; 
  usage += "\t[@stride  <write every nth frame> (default: 1)]\n"; 
  usage += "\t[@frames  <select frames to write out, starts at 0> (default: 0)]\n"; 
  usage += "\t@traj     <input trajectory files>\n";


  try {
    Arguments args(argc, argv, knowns, usage);

    // read topology
    InTopology it(args["topo"]);
    System sys(it.system());
    int nummol = sys.numMolecules();

    
    // create Time object
    utils::Time time(args);

    // get @notimeblock argument
    bool notimeblock = false;
    if (args.count("notimeblock") >= 0)
      notimeblock = true;
      
    
    if (args.count("stride") >= 0 && args.count("frames") >= 0) {
        throw gromos::Exception("dfgrid",
              "choose either @stride or @frames");
      }
    
    vector<int> frames;
    if (args.count("frames") >0) {
            for (Arguments::const_iterator it = args.lower_bound("frames");
                it != args.upper_bound("frames"); ++it) {
          int f = atoi(it->second.c_str());
          frames.push_back(f);
        }
    }

    // parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(sys, args);
    //parse gather method
    Boundary::MemPtr gathmethod = args::GatherParser::parse(sys, sys, args);


    // read stride
    int Stride = args.getValue<int>("stride", false, 1);
    int smooth = args.getValue<int>("smooth", false, 1);
    double protect = args.getValue<double>("protect", false, -1.0);
    int max = args.getValue<int>("max", false, 1);
    double gridspacing = args.getValue<double>("gridspacing", true, 0.2); 
    double proteinoffset = args.getValue<double>("proteinoffset", true, 15);
    double proteincutoff = args.getValue<double>("proteincutoff", true, 0);

    bool writegrid = true;
    if (args.count("nogrid") >=0) writegrid = false;
    
    // read in the (virtual) target atom
    utils::AtomSpecifier atom(sys);
    {
      Arguments::const_iterator iter = args.lower_bound("atom"),
              to = args.upper_bound("atom");
      atom.addSpecifier(iter->second);
      iter++;
      if (iter != to) std::cerr << "Only one atom specifier needed, only first one read\n";
      if (atom.size() == 0) {
        throw gromos::Exception("dfgrid",
              "no atoms selecte");
      } else if (atom.size() > 1) {
        throw gromos::Exception("dfgrid",
              "more than 1 (virtual) atom in selection, need exactly one");
      }        
    }
    
    // read in atoms for which to output df distance
    utils::AtomSpecifier distatoms(sys);
    for (Arguments::const_iterator iter = args.lower_bound("distatoms");
            iter != args.upper_bound("distatoms"); ++iter) {
        distatoms.addSpecifier(iter->second);
    }


    // put atoms considered as protein and penalized into an atom specifier
    utils::AtomSpecifier proteinatoms(sys);
    for (Arguments::const_iterator iter = args.lower_bound("proteinatoms");
            iter != args.upper_bound("proteinatoms"); ++iter) {
        proteinatoms.addSpecifier(iter->second);
    }

    
    // define input coordinates
    InG96 ic;
    // read first frame to get the boxsize
    ic.open(args.lower_bound("traj")->second);
    ic >> sys;
    ic.close();

    gcore::Box box=sys.box();

    double boxK=box.K()[0], boxL=box.L()[1], boxM=box.M()[2];
    if (box.ntb() != gcore::Box::rectangular) {
        throw gromos::Exception("dfgrid",
              "only rectangular boxshape supported");
    }    
    
    
    std::vector<int> ngrid;
    std::vector<double> distance;
    
    // vector storing where the shortest path to current grid point came from
    std::vector<int> from_point;
    
    // initialize the grid
    // this will only work for rectangular boxes  
    ngrid.resize(3);  
    ngrid[0] = int(boxK / gridspacing)+1;
    ngrid[1] = int(boxL / gridspacing)+1;
    ngrid[2] = int(boxM / gridspacing)+1;

    const int ndim=ngrid[0] * 
       ngrid[1] *
       ngrid[2];
  
    // we initialize the distances to something large
    distance.resize(ndim,4*(boxK + boxL + boxM));
    from_point.resize(ndim);
    
    // parse outformat
    string ext;
    OutCoordinates *oc = OutformatParser::parse(args, ext);
    
    if (args.count("distatoms") >0)  {
      std::cout << "# ";
      for (unsigned int i=0; i<distatoms.size(); i++) {
        std::cout << distatoms.toString(i)<<" ";
      }
      std::cout << "# in protein" << std::endl;
    }

    int numFrames=0;
    int skipFrame = 0;
    for (Arguments::const_iterator iter = args.lower_bound("traj");
            iter != args.upper_bound("traj"); ++iter) {

      ic.open(iter->second);
      //ic.select("ALL");

      // loop over all frames
      while (!ic.eof()) {
        if (!notimeblock) {
          ic >> sys >> time;
        } else {
          ic >> sys;
        }
        bool writeframe=false;
        for (unsigned int j = 0; j < frames.size(); ++j) {
          if (frames[j] == numFrames) {
            writeframe=true;
            break;
          }
        }

        if ((! skipFrame && frames.size()<=0) || writeframe) {
          (*pbc.*gathmethod)();
          ofstream os;
          ostringstream pdbName;
          pdbName << "grid"<< setw(5)<< setfill('0') << numFrames << ext;
          string file=pdbName.str();
          if (writegrid) {
            os.open(file.c_str());
            oc->open(os);
            oc->select("SOLUTE");


            std::ostringstream title;
            title << iter->second << ", Frame"<< numFrames <<"\nAR: gridpoints, VA: center atom, DA: distatoms" << endl;
            oc->writeTitle(title.str());


            if (!notimeblock) {
              oc->writeTimestep(time.steps(), time.time());
            }
          }
          
          // UPDATE GRID    

          std::set<int> visited;
          std::set<int> todo;

          // find the gridpoint that is nearest to atom_i. This will be our seed
          box = sys.box();
          double boxK=box.K()[0], boxL=box.L()[1], boxM=box.M()[2];
          if (box.ntb() != gcore::Box::rectangular) {
          throw gromos::Exception("dfgrid",
              "only rectangular boxshape supported");
          }
          int nx, ny, nz;
          Vec origin(0.0,0.0,0.0);
          gmath::Vec ppos = pbc->nearestImage(origin,atom.pos(0), box);

          nx = int(-(-ppos[0] - boxK/2)/gridspacing+0.5);
          ny = int(-(-ppos[1] - boxL/2)/gridspacing+0.5);
          nz = int(-(-ppos[2] - boxM/2)/gridspacing+0.5);
          gmath::Vec startpos(nx*gridspacing - boxK/2,ny*gridspacing - boxL/2,nz*gridspacing - boxM/2);
          int start=int(nx + ngrid[0]*ny + ngrid[0]*ngrid[1]*nz);
          double dist_to_i=(startpos-ppos).abs();

          // reinitialize all points to a high value (4*(box+box+box))
          // and flag all gridpoints that are within the protein
          for(unsigned int j=0; j < distance.size(); j++)
              distance[j] = 4*(boxK+boxL+boxM);

          std::set<int> protein;
          for (unsigned int a = 0; a < proteinatoms.size(); ++a){
            gmath::Vec ppos = pbc->nearestImage(origin, proteinatoms.pos(a), box);
            int nx_min=int(-( proteincutoff - ppos[0] - boxK/2)/gridspacing);
            int nx_max=int(-(-proteincutoff - ppos[0] - boxK/2)/gridspacing)+1;
            int ny_min=int(-( proteincutoff - ppos[1] - boxL/2)/gridspacing);
            int ny_max=int(-(-proteincutoff - ppos[1] - boxL/2)/gridspacing)+1;
            int nz_min=int(-( proteincutoff - ppos[2] - boxM/2)/gridspacing);
            int nz_max=int(-(-proteincutoff - ppos[2] - boxM/2)/gridspacing)+1;
          
            for(int ix=nx_min; ix < nx_max; ix++){
              for(int iy=ny_min; iy < ny_max; iy++){
                for(int iz=nz_min; iz < nz_max; iz++){
                  gmath::Vec gpos(ix*gridspacing - boxK/2, iy*gridspacing - boxL/2, iz*gridspacing - boxM/2); 
                  gmath::Vec d = gpos - ppos;
                  gmath::Vec e = startpos - gpos;

                  if(d.abs() < proteincutoff && e.abs() > protect){
                    int j=int(ix + ngrid[0]*iy + ngrid[0]*ngrid[1]*iz);
                    protein.insert(j);                    
                  }
                }
              }
            }
          }

          int current=start;
          distance[current] = dist_to_i;
  
          if(protein.count(current)) 
            distance[current]+=proteinoffset;
  
          while(visited.size() != distance.size()) {
    
            visited.insert(current);
            todo.erase(current);
    
            double currdist = distance[current];
    
            // loop over the six neighbouring gridpoints    
            for(unsigned int i=0; i<6; i++){
              double newdist = currdist + gridspacing;

              int j=neighbour(current, i, ngrid);

              // we don't revisit grids
              if(!visited.count(j)){
	  
                // check if it is in the protein
                if(protein.count(j) && !protein.count(current)){ 
                  newdist += proteinoffset;
                }
	
                if(newdist < distance[j]){
                  distance[j]=newdist;
                  from_point[j] = current;
                }
                // and j is a candidate for the next round
                todo.insert(j);
	          }
            }
            // find the next gridpoint that has not been visited yet
            int next=current;
            // ugly initialization of the minimum, distance was initialized to 2*(box+box+box)
            double min=4*(boxK+boxL+boxM)+1;
            for(std::set<int>::const_iterator iter=todo.begin(), to=todo.end();
                     iter != to; iter++){
                if(distance[*iter] < min){
                  next=*iter;
                  min=distance[*iter];
                }
            }
            current=next;
          }

          // to avoid the artificial force due to the protein, do one round of smoothening
          // the first layer of protein grid-points gets a new distance, based on
          // the nearest solvent (you should not be getting deeper 'into' the protein anyway)
          // this is dangerous if your protein is only one grid-point deep
          std::vector<Vec> ptemp;
          for(int s = 0; s < smooth; s++){
            std::set<int> remove_from_protein;
            // loop over the protein gridpoints
            for(std::set<int>::const_iterator iter=protein.begin(), to = protein.end();
	                                  iter!=to; ++iter){
	          int j=*iter;
	          gmath::Vec gpos = get_coords(j, ngrid,gridspacing,boxK,boxL,boxM);
              // find its neighbours
              for(unsigned int i=0; i<6; i++){
                int j=neighbour(*iter, i, ngrid);
	
                // check that it is in the solvent
                if(!protein.count(j)){
                  if(distance[*iter] > (distance[j] + gridspacing)){
                    // OK, this can be shortened
                    distance[*iter]=distance[j] + gridspacing;
                    remove_from_protein.insert(*iter);
                  }
                }
              }
            }
    
            for(std::set<int>::const_iterator iter=remove_from_protein.begin(), 
                      to = remove_from_protein.end(); iter!=to; ++iter)
                protein.erase(*iter);
            }

            // collect gridpoints with distance < max
            int s_numpoints=0;
            std::vector<Vec> s_gridpoints;
            std::vector<double> s_distances;
            for(unsigned int j=0; j< distance.size(); j++) {
              if (distance[j] < max) {
	            gmath::Vec gpos = get_coords(j, ngrid,gridspacing,boxK,boxL,boxM);
                s_gridpoints.push_back(gpos);
                s_distances.push_back(distance[j]);
                s_numpoints++;
              }
            }
            
            
          
            // for each atom in distatoms find nearest grid point and distance to it
            vector<int> distatom_idx;
            std::vector<std::vector<int> > minpaths(distatoms.size());
            if (args.count("distatoms")>0) {
              std::cout << numFrames << " ";
              ostringstream ostr;
              for (unsigned int a=0; a < distatoms.size(); a++) {
                gmath::Vec pos = pbc->nearestImage(origin,distatoms.pos(a), box);
                int x = int(-(-pos[0] - boxK/2)/gridspacing+0.5);
                int y = int(-(-pos[1] - boxL/2)/gridspacing+0.5);
                int z = int(-(-pos[2] - boxM/2)/gridspacing+0.5);
                gmath::Vec gridpos(x*gridspacing - boxK/2,y*gridspacing - boxL/2,z*gridspacing - boxM/2);
                int idx=int(x + ngrid[0]*y + ngrid[0]*ngrid[1]*z);
                double dist = distance[idx]+(gridpos-pos).abs();
                std::cout << dist << " ";
                if (protein.count(idx)) ostr << distatoms.toString(a) << " ";
                
                // retrace minimal distance path
                minpaths[a].push_back(idx);
                int prev=from_point[idx];
                while (prev != start) {
                  minpaths[a].push_back(prev);
                  prev=from_point[prev];
                }
                minpaths[a].push_back(prev);
              }
              std::cout << "# " << ostr.str() << std::endl;
            }
            
 
            if (writegrid) {
            // create a system which includes the center atom, distatoms and the grid  
            System gridsys=sys;
            MoleculeTopology mt;

            // add atoms to the new molecule's topology
            // VA atom for center
            int resnum = 1;
            int atomnum = 0;
            AtomTopology at_va;
            at_va.setName("VA"); 
            mt.addAtom(at_va);
            mt.setResNum(atomnum,resnum);
            atomnum++;
            mt.setResName(resnum,"VA");
            resnum++;

            // DA atoms for distatoms, each as separate residue
            AtomTopology at_da;
            at_da.setName("DA");
            for (unsigned int i=0; i<distatoms.size();i++) {
              mt.addAtom(at_da);
              mt.setResNum(atomnum,resnum);
              mt.setResName(resnum,"DA");
              resnum++;
              atomnum++;
            }

            // PA: path atoms
            AtomTopology at_pa;
            at_pa.setName("PA");
            for (unsigned int i=0; i<distatoms.size();i++) {
              for (unsigned int j=0; j<minpaths[i].size();j++) {
                mt.addAtom(at_pa);
                if (j!=0) {
                  Bond bond(atomnum-1,atomnum);
                  mt.addBond(bond);
                }
                mt.setResNum(atomnum,resnum);
                atomnum++;
              }
              mt.setResName(resnum,"PA");
              resnum++;
            }

            // AR atoms for gridpoints
            AtomTopology at;
            at.setName("AR");
            for (int i=0; i<s_numpoints; i++) {
              mt.addAtom(at);
              mt.setResNum(atomnum,resnum);
              atomnum++;
            }
            mt.setResName(resnum,"AR");
            resnum++;

            // create the molecule
            Molecule molecule = Molecule(mt);
            molecule.initPos();
            molecule.initBfac();
            gridsys.addMolecule(molecule); 

            // add centeratom coordinates
            atomnum=0;
            gridsys.mol(nummol).pos(atomnum) = atom.pos(0);
            atomnum++;

            // add distatom coordinates        
            for (unsigned int i=0; i<distatoms.size();i++) {
              gridsys.mol(nummol).pos(atomnum) = distatoms.pos(i);
              atomnum++;
            }        

            // add path coordinates
            for (unsigned int i=0; i<distatoms.size();i++) {
              for (unsigned int j=0; j<minpaths[i].size();j++) {
                int idx = minpaths[i][j];
                gmath::Vec gpos = get_coords(idx, ngrid,gridspacing,boxK,boxL,boxM);
                gridsys.mol(nummol).pos(atomnum) = gpos;
                atomnum++;
              } 
            }       

            // add grid coordinates
            for (int i=0; i<s_numpoints; i++) {
              gridsys.mol(nummol).pos(atomnum) = s_gridpoints[i];
              gridsys.mol(nummol).setBfac(atomnum, s_distances[i]);
              atomnum++;
            }

            *oc << gridsys;
            oc->close();
            os.close();
            }
        }
        numFrames++;
        skipFrame++;
        skipFrame %= Stride;
      }
      ic.close();
    }
  }  catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}
