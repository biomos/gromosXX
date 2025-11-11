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
 * @file edyn.cc
 * Performs an essential dynamics analysis
 */

/**
 * @page programs Program Documentation
 *
 * @anchor edyn
 * @section edyn Performs an essential dynamics analysis
 * @author @ref mk
 * @date 23.8.06
 *
 * Program edyn performs an essential dynamics analysis over a (list of)
 * trajectory file(s). The covariance matrix is calculated for the specified
 * atoms and diagonalised. The eigenvalues and eigenvectors are written
 * to file, as well as information about selected eigenvalues.
 * 
 * The trajectory file(s) are subsequently analysed as projections along the
 * eigenvalues. For all of the selected eigenvalues, the atomic
 * components of the eigenvalues and the time series of the projection
 * along the eigenvalue are written to file. In addition, pdb files are
 * written with coordinates of the specified atoms at the extreme
 * values of the projection along the eigenvalue. With the \@skip flag, the   
 * usually time-consuming
 * projections can be skipped and, in this case, only the covariance matrix,
 * the eigenvalues and the eigenvectors will be printed to file.
 *
 * Most of the output of this program is written to a selection of files:
 * <table border=0 cellpadding=0>
 * <tr><td>AVE.pdb </td><td>contains the average position of the specified
 * atoms</td></tr>
 * <tr><td>EIVAL.out </td><td>contains the eigenvalues of the covariance
 * matrix</td></tr>
 * <tr><td>EIVEC.out </td><td>contains the eigenvectors of the covariance
 * matrix</td></tr>
 * <tr><td>EIFLUC.out </td><td>contains the fluctuation along the
 * eigenvectors</td></tr>
 * <tr><td>ESSDYN.out </td><td>contains the averages, fluctuations, minimum and
 * maximum values of the projections along the eigenvectors</td></tr>
 * </table>
 *
 * In addition, several files are written out for each selected eigenvalue, x:
 *
 * <table border=0 cellpadding=0>
 * <tr><td>EVCOMP_x.out </td><td>contains the atomic contributions to the
 * eigenvector</td></tr>
 * <tr><td>EVPRJ_x.out </td><td>contains the time series of the projection of
 * the trajectory along the eigenvector</td></tr>
 * <tr><td>PRJMAX_x.pdb </td><td>contains coordinates of the selected atoms,
 * displaced from the average positions along the eigenvector to the maximum
 * value of the observed projection</td></tr>
 * <tr><td>PRJMIN_x.pdb </td><td>contains coordinates of the selected atoms,
 * displaced from the average positions along the eigenvector to the minimum
 * value of the observed projection</td></tr>
 * </table>
 * <p>
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; [&lt;gather method&gt;] </td></tr>
 * <tr><td> \@atoms</td><td>&lt;@ref AtomSpecifier "atoms" to be considered&gt; </td></tr>
 * <tr><td> \@ref</td><td>&lt;reference coordinates&gt; </td></tr>
 * <tr><td> [\@eigenvalues</td><td>&lt;list of eigenvalues for which data is written&gt;] </td></tr>
 * <tr><td> \@traj</td><td>&lt;trajectory file(s)&gt; </td></tr>
 * <tr><td> [\@skip</td><td>&lt;skip the (time-consuming) projections&gt;] </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  edyn
    @topo         ex.top
    @pbc          r
    @atoms        1:CA
    @ref          exref.coo
    @eigenvalues  1 2 3 4 5 6 7 8 9 10 20 50 100
    @traj         ex.tr
 @endverbatim
 *
 * <hr>
 */

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <iostream>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/utils/Rmsd.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/fit/Reference.h"
#include "../src/fit/RotationalFit.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/gmath/Matrix.h"
#include "../src/gmath/Vec.h"
#include "../src/gromos/Exception.h"

using namespace std;
using namespace gmath;
using namespace gcore;
using namespace gio;
using namespace utils;
using namespace bound;
using namespace args;
using namespace fit;

void writePdb(const char *name, AtomSpecifier &atoms);

int main(int argc, char **argv){

  Argument_List knowns; 
  knowns << "topo" << "traj" << "atoms" << "pbc" << "ref" 
         << "eigenvalues" << "skip";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo         <molecular topology file>\n";
  usage += "\t@pbc          <boundary type>\n";
  usage += "\t@atoms        <atoms to be considered>\n";
  usage += "\t@ref          <reference coordinates>\n";
  usage += "\t[@eigenvalues <list of eigenvalues for which data is written>]\n";
  usage += "\t@traj         <trajectory files>\n";
  usage += "\t[@skip         (skip the (time-consuming) projections)]\n";

  try{
    Arguments args(argc, argv, knowns, usage);

    // read topology
    InTopology it(args["topo"]);
    System refSys(it.system());

    // System for calculation
    System sys(refSys);

    // read atoms to take into account
    AtomSpecifier atoms(sys);
    for(Arguments::const_iterator it=args.lower_bound("atoms"); 
	it!=args.upper_bound("atoms"); ++it){
      atoms.addSpecifier(it->second);
    }
    if(atoms.size()==0)
      throw gromos::Exception("edyn", "No atoms specified!");
    //Bruno: I added the option to skip the projections. Sometimes one is
    //only interested in the covariance or eigenvectors and the projections
    //are time consuming.
    // Check if skip
    bool skip = false;
    if(args.count("skip") >=0) skip = true;



    // sort atom numbers
    atoms.sort();

    int NDIM=atoms.size()*3;

    // read list of eigenvectors that we are interested in
    // the covariance matrix will have dimensions atoms.size()*3 to 
    // atoms.size()*3
    // So the maximum number of eigenvalues is atoms.size()*3
    // If the user wants higher values, warn him
    std::vector<int> sel;
    for(Arguments::const_iterator it=args.lower_bound("eigenvalues");
	it!=args.upper_bound("eigenvalues"); ++it){
      int t=atoi(it->second.c_str())-1;
      if(t<NDIM)
	sel.push_back(atoi(it->second.c_str())-1);
      else
	cout << "You specified " << atoms.size() << " atoms.\nThis leads "
	     << "to a covariance matrix with dimensions " << NDIM << "x"
	     << NDIM << "\nSo we will have at most " << NDIM
	     << " eigenvalues.\nIgnoring request for eigenvalue " << t+1 
	     << endl << endl;
    }
    cout << "Selected " << sel.size() << " eigenvalues" << endl;
    
    // read reference coordinates...
    InG96 ic;

    if(args.count("ref")==1)
      ic.open(args["ref"]);
    else
      ic.open(args["traj"]);
    
    ic >> refSys;
    ic.close();

    // Parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(refSys, args);
    //parse gather method
    Boundary::MemPtr gathmethod = args::GatherParser::parse(sys,refSys,args);
    // gather reference system
    (*pbc.*gathmethod)();
    delete pbc;
    
    Reference ref(&refSys);
    ref.addAtomSpecifier(atoms);

    // Parse boundary conditions for sys
    pbc = BoundaryParser::boundary(sys, args);

    RotationalFit rf(&ref);


    int numFrames = 0;
    int size = atoms.size();
    std::vector<Vec> pvec(size);
    std::vector<Vec> avpos(size, Vec(0.0,0.0,0.0));
    Matrix cov(size*3,size*3,0.0);
    
    cout << "reading trajectory..."<< endl;
    for(Arguments::const_iterator iter=args.lower_bound("traj");
	iter!=args.upper_bound("traj"); ++iter){
      ic.open(iter->second);
      // loop over all frames
      while(!ic.eof()){
	numFrames++;
	ic >> sys;	
	(*pbc.*gathmethod)();
	rf.fit(&sys);
	
	// calculate average positions, put coords into one array
	for (int i=0; i<size;++i) {
	  pvec[i]=atoms.pos(i);
	  avpos[i]+=atoms.pos(i);
	}  

	//build one triangle of the covariance matrix
	// Chris: This is only the <r*r> part
	for (int ii=0; ii<size;++ii){
	  for (int jj=0; jj<=ii;++jj){
	    for(int iii=0; iii<3; ++iii){
	      for(int jjj=0; jjj<3; ++jjj){
		cov(3*ii+iii,3*jj+jjj) += (pvec[ii][iii]*pvec[jj][jjj]);
	      }
	    }
	  }
	}
      }
      ic.close();
    }
    
    //average positions, write them out
    for (int i=0; i<size;++i) {
      avpos[i] /= numFrames;
      atoms.pos(i)=avpos[i];
    }
    
    char outFI[]="AVE";
    ostringstream o;
    o <<outFI<<".pdb";
    writePdb(o.str().c_str(),atoms);
    
    //check matrix dimension vs. total number of frames   
    cout << "building matrix..."<<endl;  
    if (numFrames < NDIM) {
      cout << "The number of dimensions (" << NDIM 
	   << ") is larger than the total number of frames (" 
	   << numFrames << ")!\n"
	   << "This may lead to poor results.\n" << endl;}
    
    // now finish the covariance matrix by subtracting the <r><r> part
    for (int ii=0;ii<size;++ii){
      for (int jj=0;jj<=ii;++jj){
	for(int iii=0; iii<3; ++iii){
	  for(int jjj=0; jjj<3; ++jjj){
	    cov(3*ii+iii,3*jj+jjj) = 
	      cov(3*ii+iii,3*jj+jjj)/numFrames -(avpos[ii][iii]*avpos[jj][jjj]);
	  }
	}
      }
    }  
    
    //calculate trace of the symmetric matrix
    double tcov=0;
    for (int z=0;z<NDIM;++z){
      tcov+=cov(z,z);
    }
    if (tcov < 0){ 
      throw gromos::Exception("edyn", 
	    " trace of the covariance matrix is negative. "
	    "Cannot go on. Something might be wrong with your trajectory.\n");
    }
    
    //build up the other triangle of the cov matrix
    for (int ii=0;ii<NDIM;++ii){
      for (int jj=0;jj <= ii;++jj){
	cov(jj,ii)=cov(ii,jj);
      }
    }

    // now print out the covariance matrix
    ofstream ocov("COVAR.out");
    ocov << "# Covariance matrix of " << NDIM << " x " << NDIM << endl;
    ofstream ocov2("COVATOM.out");
    ocov2 << "# Covariance matrix reduced to atomic correlations\n";
    
    for(int i=0; i < NDIM; ++i){
      for(int j=0; j< NDIM; ++j){
	ocov << setw(15) << cov(i,j) << endl;
      }
    }
    for(int i=0; i< NDIM; i+=3){
      for (int j=0; j<NDIM; j+=3){
	double fac=0.0;
	
	for(int k=0; k< 3; ++k){
	  fac += cov(i+k,j+k);
	}
	ocov2 << setw(15) << fac << endl;
      }
    }
    
    ocov.close();
    ocov2.close();
    

    cout << "diagonalizing matrix..." << endl;
    
    //diagonalize the matrix
    std::vector<double> eigen(NDIM, 0.0);
    cov.diagonaliseSymmetric(&eigen[0]);
    
    //calculate trace of the diagonalized matrix
    double  tdcov=0;
    for (int z=0;z<NDIM;++z){	 
      tdcov+=eigen[z];
    }
    if (tdcov < 0){
      throw gromos::Exception("edyn", 
            " trace of the diagonalized matrix "
	    "is negative. Cannot go on. Something might "
            "be wrong with your trajectory.\n");}
    
    
    //compare traces
    if (abs(tdcov-tcov) > (0.01*(tdcov+tcov))) {
      throw gromos::Exception("edyn", " trace of "
	   "the covariance and the diagonalized matrix "
           "deviate too much. Cannot go on. Something went "
           "wrong during the diagonalization. Check your "
           "trajectory.\n");}

    //spit out eigenvalues         
    ofstream oeig;
    oeig.open("EIVAL.out");
    oeig << "Eigenvalues" << endl;
    for (int i=0;i<NDIM;++i){         
      oeig <<  i+1 << " " << eigen[i] << endl;
    }
    oeig.close();
    
    //spit out relative fluctuations
    ofstream orel;
    orel.open("EIFLUC.out");       
    orel << "Relative Fluctuations of the Eigenvalues" << endl;
    double refl=0;
    for (int i=0;i<NDIM;++i){
      refl+=(eigen[i]/tdcov);
      orel << i+1 << " " << refl << endl;
    }
    orel.close();
    
    //spit out eigenvectors
    ofstream oeiv;
    oeiv.open("EIVEC.out");
    oeiv << "Eigenvectors" << endl;  
    for (int ii=0, x=0;ii<=NDIM;++ii){
      for (int jj=0;jj<NDIM;++jj){
	double eivec0 = cov(jj,ii);
	oeiv.setf(ios::right, ios::adjustfield);     
	oeiv << setw(13) << eivec0 << " ";
	x+=1;
	if (x == 6){oeiv << endl;x=0;}
      }
    }
    oeiv.close();

    if(skip == false) {
      //eigenvector components of the selected eigenvectors
      for(unsigned int i = 0; i < sel.size(); ++i) {
        ostringstream out;
        ofstream outf;
        out << "EVCOMP_" << sel[i] + 1 << ".out";
        outf.open(out.str().c_str());
        for(int j = 0; j < size; ++j) {
          double tmp = 0;
          for(int k = 0; k < 3; ++k) {
            tmp = tmp + cov((3 * j + k), sel[i]) * cov((3 * j + k), sel[i]);
          }
          tmp = sqrt(tmp);
          outf << atoms.gromosAtom(j) + 1 << "  " << tmp << endl;
        }
        outf.close();
      }

      //prepare for next loop
      std::vector<double> minprj(sel.size(), 100000.0),
              maxprj(sel.size(), -100000.0),
              avs(sel.size(), 0.0),
              avsq(sel.size(), 0.0),
              sig(sel.size(), 0.0);

      //put selected EV's in separate Matrix
      cout << "putting EV's in EIG\n";
      Matrix eig(NDIM, sel.size(), 0.0);
      for(unsigned int i = 0; i < sel.size(); ++i) {
        for(int j = 0; j < NDIM; ++j) {
          eig(j, i) = cov(j, sel[i]);
        }
      }


      //start loop
      numFrames = 0;
      for(Arguments::const_iterator iter = args.lower_bound("traj");
              iter != args.upper_bound("traj"); ++iter) {
        ic.open(iter->second);

        // loop over all frames
        while(!ic.eof()) {
          numFrames++;
          ic >> sys;
          (*pbc.*gathmethod)();
          rf.fit(&sys);

          //substract average from frame
          for(int i = 0; i < size; ++i) {
            pvec[i] = atoms.pos(i) - avpos[i];
          }

          //  Matrix mproj (numFrames,13,0.0);
          for(unsigned int i = 0; i < sel.size(); ++i) {
            double proj = 0.0;
            for(int j = 0; j < size; ++j) {
              for(int k = 0; k < 3; ++k) {
                proj += pvec[j][k] * eig(3 * j + k, i);
              }
            }
            proj = -proj;
            int f = sel[i];
            ofstream outfile;
            ostringstream ou;
            ou << "EVPRJ_" << f + 1 << ".out";
            outfile.open(ou.str().c_str(), ofstream::app);
            outfile << numFrames << ' ' << proj << endl;
            outfile.close();

            maxprj[i] = ((proj)>(maxprj[i]) ? (proj) : (maxprj[i]));
            minprj[i] = ((proj)<(minprj[i]) ? (proj) : (minprj[i]));

            //distribution properties
            avs[i] += proj;
            avsq[i] += proj*proj;
          }
        }
      }


      // Chris: what is this about? It seems to print a linear increase from
      //        minprj (- abit) to maxprj (+abit) in 60 steps?
      //        It seems to be a halfhearted attempt to write a distribution
      //        of the projection. Let's skip it.
      /*
      double dx=0.0;
      for (int i=0;i<sel.size();++i){
        ofstream disout;
        ostringstream di;
        char disOut[]="DXPRJ";
        di <<disOut <<"_"<<sel[i]+1;
        disout.open(di.str().c_str());
        gmath::Distribution dis(minprj[i], maxprj[i], 50);
      
        for (int z=0;z<60;z++){
          dx=(z-4)*((maxprj[i]-minprj[i])/50.0)+minprj[i];
          disout << dx << endl;
        }
        disout.close();
      }
       */

      //determine averages and write out the extreme structures
      ofstream output;
      output.open("ESSDYN.out");
      output << "# Projection of trajectory along the eigenvectors\n";
      output << "#\n";
      output << "#  EV"
              << setw(15) << "average"
              << setw(15) << "fluctuation"
              << setw(15) << "min. proj."
              << setw(15) << "max. proj."
              << endl;

      for(unsigned int i = 0; i < sel.size(); ++i) {
        avs[i] = avs[i] / numFrames;
        avsq[i] = avsq[i] / numFrames;
        // Chris: this does not make sense, it should probably be
        //sig[i]=sqrt((avsq[i]-avs[i]));//*(avsq[i]-avs[i]));
        sig[i] = sqrt(avsq[i] - avs[i] * avs[i]);

        output << setw(4) << sel[i] + 1
                << setw(15) << avs[i]
                << setw(15) << sig[i]
                << setw(15) << minprj[i]
                << setw(15) << maxprj[i]
                << endl;

        // obtain max and min coordinates
        for(int j = 0; j < size; ++j) {
          for(int k = 0; k < 3; ++k) {
            atoms.pos(j)[k] = avpos[j][k] + eig(3 * j + k, i) * maxprj[i];
          }
        }
        char outFile[] = "PRJMAX";
        ostringstream out;
        out << outFile << "_" << sel[i] + 1 << ".pdb";
        writePdb(out.str().c_str(), atoms);

        for(int j = 0; j < size; ++j) {
          for(int k = 0; k < 3; ++k) {
            atoms.pos(j)[k] = avpos[j][k] + eig(3 * j + k, i) * minprj[i];
          }
        }
        char outF[] = "PRJMIN";
        ostringstream ou;
        ou << outF << "_" << sel[i] + 1 << ".pdb";
        writePdb(ou.str().c_str(), atoms);

      }
      output.close();
    }
            
  }
  
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

void writePdb(const char *name, AtomSpecifier &atoms){ 
  ofstream out;
  out.open(name);
  out.setf(ios::fixed, ios::floatfield);
  out.setf(ios::unitbuf);
  out.precision(3);
  for (unsigned int h=0; h<atoms.size(); ++h){
    int res=atoms.resnum(h);
    out << "ATOM";
    out.setf(ios::right, ios::adjustfield);
    out << setw(7) << atoms.gromosAtom(h)+1;
    out.setf(ios::left, ios::adjustfield);
    out << "  " <<setw(4) << atoms.name(h);
    out << setw(4) << atoms.resname(h);
    out.setf(ios::right, ios::adjustfield);
    out << setw(5) << res +1 << "    "
	<< setw(8) << atoms.pos(h)[0]*10
	<< setw(8) << atoms.pos(h)[1]*10
	<< setw(8) << atoms.pos(h)[2]*10
	<< "  1.00  0.00" << endl;
  }
  out.close();
}
