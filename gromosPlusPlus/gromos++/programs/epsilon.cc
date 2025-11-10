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
 * @file epsilon.cc
 * Calculate the dielectric properties of the system
 */

/**
 * @page programs Program Documentation
 *
 * @anchor epsilon
 * @section epsilon Calculate the dielectric properties of the system
 * @author @ref dp
 * @date 16-04-14
 *
 * Program epsilon calculates the dielectric properties of the system. For systems containing only neutral molecules, it estimates the relative dielectric permittivity, 
 * @f$\epsilon(0)@f$, of simulation box from a Kirkwood - Fr&ouml;hlich type of
 * equation, as derived by Neumann [Mol. Phys. 50, 841 (1983)]
 *
 * @f[ (\epsilon(0) - 1) \frac{2\epsilon_{RF}+1}{2\epsilon_{RF}+\epsilon(0)} = \frac{<\vec{M}^2> - <\vec{M}>^2}{3\epsilon_0 Vk_BT} @f]
 *
 * where @f$\vec{M}@f$ is the total dipole moment of the system, 
 * @f$\epsilon_0@f$ is the dielectric permittivity of vacuum,
 * @f$\epsilon_{RF}@f$ is a reaction-field epsilon value, @f$V@f$ is the volume
 * and @f$k_BT@f$ is the absolute temperature multiplied by the Boltzmann 
 * constant. 
 *
 * For systems containing ionic species, the total dipole moment is split into a rotational
 * @f$\vec{M_d}@f$ and a translational part @f$\vec{M_j}@f$:
 * @f[\vec{M_d} = \sum_{m=1}^{N_m}\sum_{a=1}^{N_{m,a}}q_{m,a}(\vec{r}_{m,a}-\vec{r}_{cm,m})@f]
 * @f[\vec{M_j} = \sum_{m=1}^{N_m}q_{m}\vec{r}_{cm,m}@f]
 * where  @f$N_m@f$ is the total number of molecules, @f$N_{m,a}@f$ is the number
 * of atoms of the molecule m and cm stands for the center of mass.
 *
 * The generalized frequency-dependent dielectric constant can be decomposed
 * into the following contributions:
 *
 * @f[<\vec{M_d}^2>@f]
 * autocorrelation functions @f[<\vec{M_d}(\tau)\vec{M_d}(\tau+t)>@f]@f[<\vec{J}(\tau)\vec{J}(\tau+t)>@f]
 * as well as the cross term @f[<\vec{M_d}(\tau)\vec{J}(\tau+t)>@f]
 * where @f$\vec{J} = \frac{d\vec{M_j}}{dt} = \sum_{m=1}^{N_m}q_{m}\vec{v}_{cm,m}@f$.
 *
 * Using fit functions one can calculate both the
 * frequency-dependent dielectric response and the
 * static dielectric constant of the system.
 * For more details see [Schr&ouml;der & Steinhauser, J. Chem. Phys. 132, 244109 (2010); doi. 10.1063/1.3432620]
 *
 * Practically, to calculate the static dielectric constant
 * the contribution of the cross term can be neglected, while
 * the contribution from JJ autocorrelation can be calculated
 * from @f$<\Delta\vec{M_j}^2>@f$ (the mean square displacement
 * of @f$\vec{M_j}@f$) using the Einstein relation.
 * For more details see [Schr&ouml;der, J. Chem. Phys. 135, 024502 (2011); doi 10.1063/1.3601750]
 *
 * Note that @f$<\Delta\vec{M_j}^2>@f$ has to be calculated from an unfolded trajectory.
 * For that reason, the first frame of the trajectory is gathered using
 * the gbond (to avoid broken molecules) and the rest with the gtime method.

 * The program outputs relative permittivity (epsilon) calculated exclusively
 * from the @f$<\vec{M_d}^2>@f$ contribution. For systems containing ionic species,
 * @f$<\Delta\vec{M_j}^2>@f$ is calculated and written in Mj2.out,
 * from which the translational contribution to the relative permittivity can be calculated.
 * See fit_Mj2.py script in gromos++/examples/ directory.
 *
 * Optionally, the program can calculate @f$<\vec{M_d}(\tau)\vec{M_d}(\tau+t)>@f$ 
 * and @f$<\vec{J}(\tau)\vec{J}(\tau+t)>@f$ autocorrelation functions as well as
 * the @f$<\vec{M_d}(\tau)\vec{J}(\tau+t)>@f$ cross term (MdMd.out, JJ.out and MdJ.out).
 * Note that velocity trajectory files have to be provided to calculate the JJ and MdJ contributions.
 * The translational contribution to the relative permittivity can be calculated from
 * JJ.out. See fit_JJ.py script in gromos++/examples/ directory.
 *
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; [&lt;gather method&gt;] </td></tr>
 * <tr><td> [\@time</td><td>&lt;@ref utils::Time "time and dt"&gt;] </td></tr>
 * <tr><td> [\@e_rf</td><td>&lt;reaction field epsilon&gt;] </td></tr>
 * <tr><td> \@temp</td><td>&lt;temperature&gt; </td></tr>
 * <tr><td> \@traj</td><td>&lt;trajectory files&gt; </td></tr>
 * <tr><td> [\@traj_vel</td><td>&lt;velocity trajectory files&gt;] </td></tr>
 * <tr><td> [\@omega</td><td>&lt;enable omega-dependent calculation&gt;] </td></tr>
 * <tr><td> [\@MdJ</td><td>&lt;enable cross-term MdJ calculation (usually neglected)&gt;] </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  epsilon
    @topo  ex.top
    @pbc   r
    [@time  0 0.2]
    @e_rf  61
    @temp  300
    @traj  ex.tr
 @endverbatim
 *
 * <hr>
 */

#include <cstdlib>
#include <sstream>
#include <string>
#include <vector>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <cassert>
#include <algorithm>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/fit/Reference.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/Box.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/fit/PositionUtils.h"
#include "../src/gmath/Vec.h"
#include "../src/gmath/Physics.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/utils/groTime.h"
#include "../src/gromos/Exception.h"

using namespace std;
using namespace fit;
using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;

// declaring some variables
double sum_vol=0, vcorr=1.0;
vector<Vec> vec_Md;
vector<Vec> vec_Mj;
vector<Vec> vec_Mj_p;
Vec sum_Md(0.0,0.0,0.0);
vector<double> vec_time,vec_time_v;
double sum_Md2 = 0.0;
int numFrames=0,numFrames_v=0;
vector<double> mol_net_ch;
vector<int> mol_ions;
vector<double> sol_net_ch;
vector<int> sol_ions;
bool ions_insys = false;

// we also need these factors    
double e_rf, f;

//greatest common divisor
double GCD(double a, double b){
	while (a>0.00001 and b>0.00001){
		if (a>b)
			a=a-b;
		else
			b=b-a;
	}
	return a+b;
}
// calculate & update
void calc_update(System sys,utils::Time time,utils::AtomSpecifier neut_mol){
//calculate
	sys.box().update_triclinic();
        double vol=vcorr*sys.box().K_L_M();
	sum_vol+=vol;
	
	// calculate Md, Mj and Mj_p
	Vec Md(0.0,0.0,0.0);
	Vec Mj(0.0,0.0,0.0);

	//loop over neutral molecules
	for (unsigned int a=0;a<neut_mol.size();a++)
		Md+=(neut_mol.pos(a))*neut_mol.charge(a);

	// loop over ions
	utils::AtomSpecifier temp_atoms(sys);
	Vec temp_com;
	// solute
	int m;
        for(unsigned int ch_m=0; ch_m<mol_ions.size(); ++ch_m) {
		m=mol_ions[ch_m];
		temp_atoms.clear();
		temp_atoms.addMolecule(m);
		temp_com=PositionUtils::com(sys, temp_atoms);
		for(unsigned int a=0; a < temp_atoms.size(); a++) {
			Md+=(temp_atoms.pos(a)-temp_com)*temp_atoms.charge(a);
		}
		Mj+=temp_com*mol_net_ch[ch_m];
	}

	// solvent
	int temp_sol_atom=0;
	string temp_spec;
	ostringstream ostrs;
	unsigned int ch_s;
        for(int s=0; s<sys.numSolvents(); ++s) {
		ch_s=find(sol_ions.begin(), sol_ions.end(), s)-sol_ions.begin();
		if (ch_s != sol_ions.size()){
			for (int m=0; m<(sys.sol(s).numAtoms()/sys.sol(s).topology().numAtoms()); ++m) {
				temp_atoms.clear();
				ostrs.str("");
				ostrs << temp_sol_atom+1;
				temp_spec="s:"+ostrs.str()+"-";
				temp_sol_atom+=sys.sol(s).topology().numAtoms();
				ostrs.str("");
				ostrs << temp_sol_atom;
				temp_spec+=ostrs.str();
				temp_atoms.addSpecifier(temp_spec);
				temp_com=PositionUtils::com(sys, temp_atoms);
				for(unsigned int a=0; a < temp_atoms.size(); a++) {
					Md+=(temp_atoms.pos(a)-temp_com)*temp_atoms.charge(a);
				}
				Mj+=temp_com*sol_net_ch[ch_s];
			}
		}
	}
//update
        vec_Md.push_back(Md);
        vec_Mj.push_back(Mj);
	vec_time.push_back(time.time());

	numFrames++;
        sum_Md2+= Md.abs2();

	// calculate the average sum_Md2
	//double Md2_avg=sum_Md2 / numFrames;

	// calculate the current estimate of eps
	sum_Md+=Md;
	double fluc=sum_Md2 / numFrames - (sum_Md / numFrames).abs2();
	double fac, a, b;
	fac = f * sum_vol/numFrames;
        a = 2 * e_rf * fluc + fac;
	b= -fluc + fac;
	double eps = a/b;

//output
	cout << time
	     << setw(15) << setprecision(8) << Md.abs()
	     << setw(15) << setprecision(8) << sum_Md2 / numFrames
	     << setw(15) << setprecision(8) << (sum_Md / numFrames).abs2()
	     << setw(15) << setprecision(8) << eps
	     << endl;
}

void calc_update_Mj_p(System sys,utils::Time time,utils::AtomSpecifier neut_mol){
	Vec Mj_p(0.0,0.0,0.0);
	// loop over ions
	utils::AtomSpecifier temp_atoms(sys);
	Vec temp_com_v;
	// solute
	int m;
        for(unsigned int ch_m=0; ch_m<mol_ions.size(); ++ch_m) {
		m=mol_ions[ch_m];
		temp_atoms.clear();
		temp_atoms.addMolecule(m);
		temp_com_v=PositionUtils::com_v(sys, temp_atoms);
		Mj_p+=temp_com_v*mol_net_ch[ch_m];
	}

	// solvent
	int temp_sol_atom=0;
	string temp_spec;
	ostringstream ostrs;
	int ch_s;
        for(int s=0; s<sys.numSolvents(); ++s) {
		ch_s=find(sol_ions.begin(), sol_ions.end(), s)-sol_ions.begin();
		if (ch_s != (int)sol_ions.size()){
			for (int m=0; m<(sys.sol(s).numAtoms()/sys.sol(s).topology().numAtoms()); ++m) {
				temp_atoms.clear();
				ostrs.str("");
				ostrs << temp_sol_atom+1;
				temp_spec="s:"+ostrs.str()+"-";
				temp_sol_atom+=sys.sol(s).topology().numAtoms();
				ostrs.str("");
				ostrs << temp_sol_atom;
				temp_spec+=ostrs.str();
				temp_atoms.addSpecifier(temp_spec);
				temp_com_v=PositionUtils::com_v(sys, temp_atoms);
				Mj_p+=temp_com_v*sol_net_ch[ch_s];
			}
		}
	}
//update
        vec_Mj_p.push_back(Mj_p);
	numFrames_v++;
	vec_time_v.push_back(time.time());
}


int main(int argc, char **argv){
  Argument_List knowns;
  knowns << "topo" << "pbc" << "temp" << "e_rf" << "traj" << "time" << "traj_vel" << "omega" << "MdJ";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo <molecular topology file>\n";
  usage += "\t@pbc    <boundary type> [<gather method>]\n";
  usage += "\t@temp   <temperature>\n";
  usage += "\t[@e_rf   <reaction field epsilon>]\n";
  usage += "\t@traj   <trajectory files>\n";
  usage += "\t[@time   <time and dt> - applies to both traj and traj_vel]\n";
  usage += "\t[@traj_vel   <velocity trajectory files>]\n";
  usage += "\t[@omega  <enable omega-dependent calculation>]\n";
  usage += "\t[@MdJ    <enable cross-term MdJ calculation (usually neglected)>]\n";
  
 
  try{
// reading arguments
    Arguments args(argc, argv, knowns, usage);

    // read the temperature
    double temp = args.getValue<double>("temp");

    //  read topology
    InTopology it(args["topo"]);
    System sys(it.system());
    System refSys(it.system());

    // read e_rf
    e_rf= args.getValue<double>("e_rf", false, 1.0);
    	 
    // parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(sys, args);
    Boundary::MemPtr gathmethod = args::GatherParser::parse(sys,refSys,args);
    // parse gather method

    Arguments::const_iterator iter=args.lower_bound("pbc");
    if(iter!=args.upper_bound("pbc"))
      if(iter->second[0]=='t') vcorr=0.5;

    // parse omega and MdJ
    bool omega_flag = false;
    bool MdJ_flag = false;
    if (args.count("omega") >= 0)
      omega_flag = true;
    if (args.count("MdJ") >= 0)
      MdJ_flag = true;

    // get the @time argument
    utils::Time time(args);
  
    f=3.0 * gmath::physConst.get_eps0() * gmath::physConst.get_boltzmann() * temp*(2*e_rf + 1.0);

    // declaration of some temp variables

    // define input coordinate and reading in the first frame
    InG96 ic;
    iter=args.lower_bound("traj");
    ic.open((iter->second).c_str());
    ic.select("ALL");
    ic >> sys >> time;
    ic.close();

    // create atomspecifier for neutral molecules
    utils::AtomSpecifier neut_mol(sys);

    // find ions and their net charges - solute molecules
     for(int m=0; m<sys.numMolecules(); ++m) {
	double tc=0.0;
	for(int a=0; a< sys.mol(m).numAtoms(); ++a)
	  tc+=sys.mol(m).topology().atom(a).charge();
	if(tc!=0.0){
          ions_insys = true;
          mol_ions.push_back(m);
	  mol_net_ch.push_back(tc);
        }
	else {neut_mol.addMolecule(m);}
      }
    // find ions and their net charges - solvent molecules
     for(int m=0; m<sys.numSolvents(); ++m) {
	double tc=0.0;
	for(int a=0; a< sys.sol(m).topology().numAtoms(); ++a)
	  tc+=sys.sol(m).topology().atom(a).charge();
	if(tc!=0.0){
          sol_ions.push_back(m);
	  sol_net_ch.push_back(tc);
        }
      }

    if (sol_ions.empty())
      neut_mol.addSpecifier("s:a");
    else{
      int temp_sol_atom=0;
      string temp_spec;
      ostringstream ostrs;
      for(int s=0; s<sys.numSolvents(); ++s) {
        if (!(find(sol_ions.begin(), sol_ions.end(), s) != sol_ions.end())){
          ostrs.str("");
          ostrs << temp_sol_atom+1;
          temp_spec="s:"+ostrs.str()+"-";
          temp_sol_atom+=sys.sol(s).numAtoms();
          ostrs.str("");
          ostrs << temp_sol_atom;
          temp_spec+=ostrs.str();
          neut_mol.addSpecifier(temp_spec);
        }
	else{
          temp_sol_atom+=sys.sol(s).numAtoms();
	  if (sys.sol(s).numAtoms()!=0)
	          ions_insys = true;
	}
      }
    }

    if (ions_insys){
// there are ions in the system
	// taking care of the gathering method - we have to gather the first frame with gbond and the rest with gtime
        Arguments::const_iterator temp_args_iter=args.lower_bound("pbc");
	if(temp_args_iter!=args.upper_bound("pbc")){
		temp_args_iter++;
		if (string(temp_args_iter->second)!="2" && temp_args_iter->second!="gtime")
			throw gromos::Exception("epsilon","gtime (2) gathering method has to be used when the system contains ions");
	}
	else{
		throw gromos::Exception("epsilon","gtime (2) gathering method has to be used when the system contains ions");
	}
	// gather with gbond
	gathmethod = &Boundary::gatherbond;
	(*pbc.*gathmethod)();
	// now reset refsys
	pbc->setReferenceSystem(sys);
	gathmethod = &Boundary::gathertime;
	cout << "# exceptionally, the first frame gathered with gbond and the rest with gtime\n";
    }

// calculations
    // write a title

    cout << "#\n#"
	     << setw(14) << "time"
	     << setw(15) << "dipole"
	     << setw(15) << "<dipole^2>"
	     << setw(15) << "<dipole>^2"
	     << setw(15) << "epsilon" << endl;

    // loop over all trajectories
    iter=args.lower_bound("traj");
    Arguments::const_iterator to=args.upper_bound("traj");
    for(; iter!=to; ++iter){
      
      // open file
      ic.open((iter->second).c_str());
      ic.select("ALL");
      
      // loop over single trajectory
      while(!ic.eof()){
	ic >> sys >> time;
	(*pbc.*gathmethod)();
	calc_update(sys,time,neut_mol);
      }
      ic.close();
    }

// final data
    cout << "#\n#"
	     << setw(14) << "2e_rf"
	     << setw(15) << "<dipole^2>"
	     << setw(15) << "<dipole>^2"
	     << setw(20) << "3e_0<V>kT(2e_rf+1)"
	     << endl;
    cout << "#"
	     << setw(14) << setprecision(8) << 2*e_rf
	     << setw(15) << setprecision(8) << sum_Md2 / numFrames
	     << setw(15) << setprecision(8) << (sum_Md / numFrames).abs2()
	     << setw(20) << setprecision(8) << f * sum_vol/numFrames
	     << endl;

// velocity analysis
    if (args.count("traj_vel") >= 0){
	utils::Time time(args);
	System sys(it.system());
    // loop over all trajectories
	InG96 ic_v;
	Arguments::const_iterator iter=args.lower_bound("traj_vel");
	Arguments::const_iterator to=args.upper_bound("traj_vel");
	for(; iter!=to; ++iter){
	// open file
		ic_v.open((iter->second).c_str());
		ic_v.select("ALL");
		// loop over single trajectory
		while(!ic_v.eof()){
			ic_v >> sys >> time;
			calc_update_Mj_p(sys,time,neut_mol);
		}
		ic_v.close();
	}
    }

// post analysis
    if (ions_insys or omega_flag or MdJ_flag){
	//check if dt equal throughout
	double dt=vec_time[1]-vec_time[0];
	for (int t_step=2;t_step<numFrames; ++t_step) {
		if ((vec_time[t_step]-vec_time[t_step-1]-dt)>0.00001 or (vec_time[t_step]-vec_time[t_step-1]-dt)<-0.00001)
			throw gromos::Exception("epsilon","dt not equal throughout the trajectory");
	}
    }

    if (ions_insys){
	ofstream Mj2_out;
	Mj2_out.open("Mj2.out");
	Mj2_out << "# Time series of the mean square displacement of Mj" << endl;
	for (int d_t=0;d_t<numFrames; ++d_t) {
		double temp_Mj2=0.0;
		double dt=vec_time[d_t]-vec_time[0];
		int st_fr=0;
		for (;st_fr<(numFrames-d_t); ++st_fr){
			temp_Mj2+=(vec_Mj[st_fr+d_t]-vec_Mj[st_fr]).abs2();
		}
		Mj2_out << setw(15) << setprecision(8) << dt
			<< setw(15) << setprecision(8) << temp_Mj2/st_fr
			<< endl;
	    }
	Mj2_out.close();
    }

// if omega dependent analysis is enabled
    if (omega_flag){
	ofstream MdMd_out;
	MdMd_out.open("MdMd.out");
	MdMd_out << "#"
	     << setw(14) << "time"
	     << setw(20) << "<Md(tau)Md(tau+t)>" << endl;
	for (int d_t=0;d_t<numFrames; ++d_t) {
		double temp_MdMd=0.0;
		double dt=vec_time[d_t]-vec_time[0];
		for (int st_fr=0;st_fr<(numFrames-d_t); ++st_fr){
			temp_MdMd+=vec_Md[st_fr].dot(vec_Md[st_fr+d_t]);
		}
		MdMd_out << setw(15) << setprecision(8) << dt
			<< setw(20) << setprecision(8) << temp_MdMd/(numFrames-d_t)
			<< endl;
	    }
	MdMd_out.close();
    }

    if (args.count("traj_vel") >= 0){
	//check if dt equal throughout
	double dt=vec_time_v[1]-vec_time_v[0];
	for (int t_step=2;t_step<numFrames_v; ++t_step) {
		if ((vec_time_v[t_step]-vec_time_v[t_step-1]-dt)>0.00001 or (vec_time_v[t_step]-vec_time_v[t_step-1]-dt)<-0.00001)
			throw gromos::Exception("epsilon","dt not equal throughout the velocity trajectory");
	}
	ofstream JJ_out;
	JJ_out.open("JJ.out");
	JJ_out << "#"
	     << setw(14) << "time"
	     << setw(20) << "<J(tau)J(tau+t)>" << endl;
	for (int d_t=0;(vec_time_v[d_t]-vec_time_v[0])<30; ++d_t) {
		double temp_JJ=0.0;
		dt=vec_time_v[d_t]-vec_time_v[0];
		for (int st_fr=0;st_fr<(numFrames_v-d_t); ++st_fr){
			temp_JJ+=vec_Mj_p[st_fr].dot(vec_Mj_p[st_fr+d_t]);
		}
		JJ_out << setw(15) << setprecision(8) << dt
			<< setw(20) << setprecision(8) << temp_JJ/(numFrames_v-d_t)
			<< endl;
	}
	if (MdJ_flag){
		if (vec_time[0]!=vec_time_v[0])
			throw gromos::Exception("epsilon","the initial time (the first frame) of the trc and trv files not the same");
		dt=vec_time[1]-vec_time[0];
		double dt_v=vec_time_v[1]-vec_time_v[0];
		double gcd=GCD(dt,dt_v);
		int step = dt_v/gcd;
		int step_v = dt/gcd;
		ofstream MdJ_out;
		MdJ_out.open("MdJ.out");
		MdJ_out << "#"
		     << setw(14) << "time"
		     << setw(20) << "<Md(tau)J(tau+t)>" << endl;
		for (int d_t=0;d_t<numFrames_v; ++d_t) {
			double temp_MdJ=0.0;
			dt=vec_time_v[d_t]-vec_time[0];
			int st_fr=0;
			for (;st_fr<numFrames and (st_fr*step_v/step)<(numFrames_v-d_t);st_fr=st_fr+step){
				int st_fr_v=st_fr*step_v/step;
				temp_MdJ+=vec_Md[st_fr].dot(vec_Mj_p[st_fr_v+d_t]);
			}
			MdJ_out << setw(15) << setprecision(8) << dt
				<< setw(20) << setprecision(8) << temp_MdJ/(st_fr/step)
				<< endl;
		    }
		MdJ_out.close();
	}
    }
    else
	if (MdJ_flag)
		throw gromos::Exception("epsilon","MdJ not calculated - only in combination with @traj_vel");
  }
  
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

