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
 * @file ext_ti_ana.cc
 * predict free energy derivatives at non-simulated lambda values
 * to specific properties
 */

/**
 * @page programs Program Documentation
 *
 * @anchor ext_ti_ana
 * @section ext_ti_ana predict free energy derivatives at non-simulated lambda values
 * @author @ref adr @ref mp
 * @date 16. 3. 2016
 *
 * 
 * From a simulation at a single coupling parameter @f$\lambda@f$ program ext_ti_ana 
 * predicts free energy derivatives @f$\frac{\partial H}{\partial \lambda} @f$ over a 
 * range of @f$\lambda@f$ values. To do this it requires that the terms described in 
 * de Ruiter & Oostenbrink [JCTC 2016, 12, 4476-4486] have been precalculated during 
 * the simulation and written to the energy trajectories (gromos md++ block 
 * PRECALCLAM).
 *
 * In GROMOS the coupling parameter of different interaction types x, @f$\mu_x@f$, 
 * can be set individually from a fourth order polynomial of the global coupling 
 * parameter @f$\lambda@f$ (gromos md++ block LAMBDAS): 
 * @f[\mu_x = a_x \lambda^4 + b_x \lambda^3 + c_x \lambda^2 + d_x \lambda + e_x. @f]
 *
 * ext_ti_ana can make predictions for other combinations of the coefficients  
 * (a,b,c,d,e) than the ones used in the simulation. The interaction properties (x) 
 * that can be given individual @f$\lambda@f$ dependencies are:
 * - slj: lennard jones softness
 * - scrf: coulomb reaction-field softness 
 * - lj: lennard jones
 * - crf: coulomb reaction-field
 * - bond: bond
 * - ang: angle
 * - impr: improper dihedral angle
 * - dih: dihedral angle
 * - kin: kinetic (this one is not implemented in the LAMBDAS block of Gromos md++)
 * - disres: distance restraint // Betty
 * - dihres: dihedral angle restraint // Betty
 * - disfld: distance field interaction // Betty
 * To facilitate the evaluation of many sets of coefficients for slj and scrf, these 
 * can be given in an input file of the following format:
 *
 * @verbatim
TITLE
..
END
SLJ
#   uniquelabel  a        b        c        d        e
    1           -0.4      0.4     -0.3      1.3      0.0
    2           -0.4      0.4     -0.2      1.2      0.0
    3           -0.4      0.4     -0.1      1.1      0.0

END
SCRF
#   uniquelabel  a        b        c        d        e
    1            0.0      0.0      0.7      0.3      0.0
    2            0.0      0.0     -0.1      0.3      0.0
END
 @endverbatim

 * Predictions will be made for any combination of every given slj with every given 
 * scrf.
 *
 * The coefficients used in the simulation are specified by the \@lamX_sim flags or 
 * read from the LAMBDAS block in a gromos input parameter file specified by \@imd. 
 * The coefficients we want to predict for are specified by the \@lamX flags. Also the 
 * temperature of the simulation, the simulated @f$\lambda@f$ (\@slam) and its exponent
 * (\@NLAMs) and the parameters of the PRECALCLAM block (\@nrlambdas, \@minlam, 
 * \@maxlam) can be read from this imd file. If parameters are specified explicitly as 
 * input flags they will always overwrite the corresponding values read from the imd
 * file.
 *
 * If the flag \@countframes is set, the number of contributing frames for each 
 * predicted @f$\lambda@f$ is appended as an additional column to the output. This 
 * number is evaluated as the number of snapshots for which the difference in the 
 * predicted energy and the simulated energy is less than the free energy difference 
 * between the two states, as calculated using the perturbation formula.
 *
 * Error estimates can be calculated using bootstrapping. A random set of data points 
 * of the size of the original set will be chosen and the predictions made. This is 
 * repeated for as many bootstrap replicates as requested. The standard deviation over
 * the bootstrap replicates is reported as a bootstrap error.
 * 
 * The predicted TI curves from several simulations at different @f$\lambda@f$ values 
 * can be combined using @ref ext_ti_merge .
 * 
 * Finally, program ext_ti_ana can write out time series of energies at alternative 
 * value of @f$\lambda@f$ to be used for free-energy estimates with Bennett's 
 * acceptance ratio (BAR). If option @bar_data is used without further input 
 * parameters, this data is written out for all predicted @f$\lambda@f$ values, or if 
 * further input parameters are given it is only written for the selected values of 
 * @f$\lambda@f$. Program @ref bar can be used to estimate free-energy differences 
 * from these files
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@en_files   </td><td>&lt;energy files&gt; </td></tr>
 * <tr><td> \@fr_files   </td><td>&lt;free energy files&gt; </td></tr>
 * <tr><td> \@library    </td><td>&lt;library for block information&gt; </td></tr>
 * <tr><td> \@temp       </td><td>&lt;simulation temperature&gt; </td></tr>
 * <tr><td> \@nrlambdas  </td><td>&lt;no. of precalculated lambdas&gt; </td></tr>
 * <tr><td> \@minlam     </td><td>&lt;minimum precalculated lambda&gt; </td></tr>
 * <tr><td> \@maxlam     </td><td>&lt;maximum precalculated lambda&gt; </td></tr>
 * <tr><td> \@slam       </td><td>&lt;lambda value of simulation&gt; </td></tr>
 * <tr><td> \@NLAMs      </td><td>&lt;lambda exponent of simulation&gt; </td></tr>
 * <tr><td> [\@lam&lt;X&gt;_sim  </td><td>&lt;a b c d e; coefficients for individual lambda dependence, default: 0 0 0 1 0&gt;] </td></tr>
 * <tr><td></td><td> (where &lt;X&gt; is one of: slj, scrf, lj, crf, bond, ang, impr, dih, kin, disres, dihres, disfld)</td></tr> // Betty
 * <tr><td> [\@imd       </td><td>&lt;gromos input parameter file>] </td></tr>
 * <tr><td> [\@NLAMp     </td><td>&lt;lambda exponent value to predict for, default: \@NLAMs&gt;] </td></tr>
 * <tr><td> [\@lam&lt;X&gt;      </td><td>&lt;a b c d e, coefficients for individual lambda dependence, default: \@lam&lt;X&gt;_sim&gt;] </td></tr>
 * <tr><td></td><td> (where &lt;X&gt; is one of: slj, scrf, lj, crf, bond, ang, impr, dih, kin, disres, dihres, disfld)</td></tr> // Betty
 * <tr><td> [\@slj_scrf_file</td><td> &lt;file with sets of slj and scrf lambda coefficients&gt;] </td></tr>
 * <tr><td> [\@no_&lt;X&gt;    </td><td>&lt;exclude &lt;X&gt; free energy derivative contribution&gt;] </td></tr>
 * <tr><td></td><td> (where &lt;X&gt; is one of: lj, crf, bond, ang, impr, dih, kin, disres, dihres, disfld)</td></tr> // Betty
 * <tr><td> [\@countframes</td><td>&lt;count nr of contributing frames for each plam&gt;] </td></tr>
 * <tr><td> [\@pmin      </td><td>&lt;min index of prediction&gt;] </td></tr>
 * <tr><td> [\@pmax      </td><td>&lt;max index of prediction&gt;] </td></tr>
 * <tr><td> [\@bootstrap </td><td>&lt;no. of bootstraps&gt;] </td></tr>
 * <tr><td> [\@outdir    </td><td>&lt;directory to write output to&gt;] </td></tr>
 * <tr><td> [\@lam_precision</td><td> &lt;lambda value precision in outfiles, default: 2&gt;] </td></tr>
 * <tr><td> [\@bar_data   </td><td>&lt;print energies to be used for BAR (not reweighted)&gt;] </td></tr>
 * <tr><td> [\@dHdl_data  </td><td>&lt;print dHdl values to be used for MBAR (not reweighted)&gt;] </td></tr>
 * <tr><td> [\@verbose   </td><td>&lt;print used parameters to file header&gt;] </td></tr>
 * <tr><td> [\@cpus   </td><td>&lt;number of omp threads, default: 1&gt;] </td></tr>
 * </table>
 *
 *@bar_data      <print energies to be used for BAR (not reweighted)>
 * Example:
 * @verbatim
  ext_ti_ana
    @en_files       ex.tre
    @fr_files       ex.trg
    @library        ene_ana.lib
    @temp           298.15
    @nrlambdas      101
    @minlam         0
    @maxlam         1
    @slam           0.3
    @NLAMs          1
    @lambond_sim    0 1 -3 3 0    
    @lamcrf_sim     0 1 0 0 0
    @lambond        -0.4 0.4 -0.3 1.3 0.0
    @slj_scrf_file  file.txt
    @bootstrap      100
    @outdir         datafolder
    @bar_data       0 0.5 1.0
    @verbose
   @endverbatim
 *
 * <hr>
 */
#include <cassert>
#include <cctype>
#include <ctime>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <cstdlib>
#include <stdlib.h>
#include <string>
#include <vector>
#include <map>
#include <cmath>

#include <sys/stat.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#ifdef OMP
#include <omp.h>
#endif

#include "mk_script.h"
#include "../src/args/Arguments.h"
#include "../src/gio/Ginstream.h"
#include "../src/gmath/Stat.h"
#include "../src/gmath/Physics.h"
#include "../src/utils/EnergyTraj.h"
#include "../src/gromos/Exception.h"

using namespace std;
using namespace args;
using namespace gio;
using namespace gcore;

struct LambdaStruct{
  int idx;
  double lam;
  double dlam;
  int idx_low;
  int idx_high;
  double wlow;
  double whigh;
};

template < typename T >
inline std::ostream& operator << (std::ostream& os, const std::vector<T>& v) 
{
    os << "[";
    for (typename std::vector<T>::const_iterator ii = v.begin(); ii != v.end(); ++ii)
    {
        os << " " << *ii;
    }
    os << " ]";
    return os;
} 


void set_standards(utils::EnergyTraj &e);
void read_library(string name, utils::EnergyTraj& e);
void read_slj_scrf_file(string name, vector<vector<LambdaStruct> > &SLJ, 
                         vector<vector<LambdaStruct> > &SCRF, 
                         vector<vector<int> > &slabels, 
                         double minlam, double lam_step, int nrlambdas, 
                         vector<double> &p_LAM, 
                         vector<vector<double> > &SLJ_abcde, 
                         vector<vector<double> > &SCRF_abcde);
bool check_etraj_version(utils::EnergyTraj &etrj, Ginstream &gin_en, Ginstream &gin_fr, bool do_free_energy_files, bool do_energy_files, std::string library);

LambdaStruct getSingleLambda(vector<double> abcde, double minlam, double lam_step,
                      int nrlambdas, double slam, vector<double> &p_LAM, std::string lambda_type="");

vector<LambdaStruct> getLambdas(vector<double> abcde, double minlam, double lam_step,
                   int nrlambdas, vector<double> &p_LAM,std::string lambda_type="");
                   
unsigned int factorial(unsigned int n) 
{
    if (n == 0)
       return 1;
    return n * factorial(n - 1);
}

bool file_exists(string file);
int dir_exists(const char *path);

std::map<std::string, vector<double> > read_lambdas(std::vector<ilambdas::lambint> &lambints, std::vector<std::string>  &lambdas_types);

int main(int argc, char **argv){
  
  Argument_List knowns; 
  knowns << "en_files" << "fr_files" << "library" << "nrlambdas" << "minlam" << "maxlam" 
         << "slam" << "NLAMs" << "NLAMp" << "temp" 
         << "slj_scrf_file" << "outdir" << "imd" 
         << "lamslj" << "lamlj" << "lamscrf" << "lamcrf" << "lamkin" << "lambond" << "lamang" << "lamimpr" << "lamdih" << "lamdisres" << "lamdihres" << "lamdisfld" // Betty
         << "lamslj_sim" << "lamlj_sim" << "lamscrf_sim" << "lamcrf_sim" << "lamkin_sim" 
         << "lambond_sim" << "lamang_sim" << "lamimpr_sim" << "lamdih_sim" << "lamdisres_sim" << "lamdihres_sim" << "lamdisfld_sim" // Betty
         << "no_lj" << "no_crf" << "no_kin" << "no_bond" << "no_ang" << "no_dih" << "no_disres" << "no_dihres" << "no_disfld" // Betty
         << "no_impr" << "pmin" << "pmax" << "bootstrap" << "countframes" << "lam_precision" << "verbose" 
         << "cpus" << "bar_data" << "dHdl_data";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@en_files       <energy files>\n";
  usage += "\t@fr_files       <free energy files>\n";
  usage += "\t@library        <library for block information> \n";
  usage += "\t@temp           <simulation temperature>\n";
  usage += "\t@nrlambdas      <no. of precalculated lambdas>\n";
  usage += "\t@minlam         <minimum precalculated lambda>\n";
  usage += "\t@maxlam         <maximum precalculated lambda>\n";
  usage += "\t@slam           <lambda value of simulation>\n";
  usage += "\t@NLAMs          <lambda exponent of simulation>\n";
  usage += "\t[@lam<X>_sim    <a b c d e; coefficients for individual lambda dependence, default: 0 0 0 1 0>]\n";
  usage += "\t                (where <X> is one of: slj, scrf, lj, crf, bond, ang, impr, dih, kin, disres, dihres, disfld)\n"; // Betty
  usage += "\t[@imd           <gromos input parameter file>] \n";
  usage += "\t                (flags temp, nrlambdas, minlam, maxlam, slam, NLAMs and lamX_sim can be\n";
  usage += "\t                skipped if they can be read from the provided imd file, from PRECALCLAM, \n";
  usage += "\t                PERTURBATION and MULTIBATH blocks.)\n";
  usage += "\t[@NLAMp         <lambda exponent value to predict for, default: @NLAMs >]\n";
  usage += "\t[@lam<X>          <a b c d e, coefficients for individual lambda dependence, default: @lam<X>_sim>]\n";
  usage += "\t                (where <X> is one of: slj, scrf, lj, crf, bond, ang, impr, dih, kin, disres, dihres, disfld)\n"; // Betty
  usage += "\t[@slj_scrf_file < file with sets of slj and scrf lambda coefficients>]\n";
  usage += "\t[@no_<X>        <exclude <X> free energy derivative contribution>]\n";
  usage += "\t                (where <X> is one of: lj, crf, bond, ang, impr, dih, kin, disres, dihres, disfld)\n"; // Betty
  usage += "\t[@countframes   <count nr of contributing frames for each plam>]\n";
  usage += "\t[@pmin          <min index of prediction>]\n";
  usage += "\t[@pmax          <max index of prediction>]\n";
  usage += "\t[@bootstrap     <no. of bootstraps>]\n";
  usage += "\t[@outdir        <directory to write output to>]\n";
  usage += "\t[@cpus           <number of threads> Default: 1]\n";
  usage += "\t[@lam_precision <lambda value precision in outfiles, default: 2>]\n";
  usage += "\t[@verbose       <print used parameters to file header>]\n";
  usage += "\t[@bar_data      <print energies to be used for BAR (not reweighted)>]\n";
  usage += "\t[@dHdl_data     <print dHdl values to be used for MBAR (not reweighted)>]\n";

  try{
    #ifdef OMP
    double start_total = omp_get_wtime();
    #endif
    Arguments args(argc, argv, knowns, usage);
    // check whether we have both energy and free energy files
    if(args.count("en_files")<=0 || args.count("fr_files")<=0)
      throw gromos::Exception("ext_ti_ana", "en_files and fr_files need to be specified:\n"+usage);

    // require a library 
    if(args.count("library") <=0)
      throw gromos::Exception("ext_ti_ana", "no library file specified:\n"+usage);
    
    // read a library file?
    string library="";
    {
      Arguments::const_iterator iter=args.lower_bound("library"), 
	to=args.upper_bound("library");
      if(iter!=to){
	library=iter->second;
	++iter;
      }
    }

    bool slj_scrf_file = false;
    if(args.count("slj_scrf_file")>0) slj_scrf_file = true;

    bool verbose = false;
    if(args.count("verbose")>=0) verbose = true;

    int bootstrap = args.getValue<int>("bootstrap",false, 0);
    // do we calculate error estimates?
    bool no_error=true;
    if(bootstrap) no_error = false;
    else no_error = true;

    int lam_precision = args.getValue<int>("lam_precision",false, 2);

    //we will do some numerical comparisons between lambda values that may lead to problems
    double lam_epsilon = pow(10, -(lam_precision+2));
    
    unsigned int num_cpus = args.getValue<unsigned int>("cpus",false, 1);

    // which contributions do we include?
    bool no_lj = true;
    if(args.count("no_lj")<0) no_lj = false;
    bool no_crf = true;
    if(args.count("no_crf")<0) no_crf = false;
    bool no_kin = true;
    if(args.count("no_kin")<0) no_kin = false;
    bool no_bond = true;
    if(args.count("no_bond")<0) no_bond = false;
    bool no_ang = true;
    if(args.count("no_ang")<0) no_ang = false;
    bool no_impr = true;
    if(args.count("no_impr")<0) no_impr = false;
    bool no_dih = true;
    if(args.count("no_dih")<0) no_dih = false;
    /** Betty
     * perturbed special ineractions
     */
    bool no_disres = true;
    if(args.count("no_disres")<0) no_disres = false;
    bool no_dihres = true;
    if(args.count("no_dihres")<0) no_dihres = false;
    bool no_disfld = true;
    if(args.count("no_disfld")<0) no_disfld = false;
    // Betty

    bool countframes;
    if(args.count("countframes")==0) countframes = true;
    else countframes = false;

    string outdir = args.getValue<string>("outdir",false, ".");
    // create directory if it doesn't exist
    int de = dir_exists(outdir.c_str());
    if (de == 0) { // outdir doesn't exist
      const int dir_err = mkdir(outdir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
      if (-1 == dir_err)
      {
        throw gromos::Exception("ext_ti_ana", "Error creating directory!\n"); 
      }
    } else if (de == 2) { // there is already a file of this name which is not a directory
      throw gromos::Exception("ext_ti_ana", "Error creating directory, file exists!\n"); 
    }
    
    // read sim. lambdas coefficients from the flags; a default 0 0 0 1 0 is set 
    // if the flag is not given, but will be re-set later if corresponding lines
    // occur in the imd file
    map< std::string, vector<double> > abcde, abcde_sim;
    // mapping the LAMBDAS block ntli numbers to our flags;
    // bond(1), angle(2), dihedral(3), improper(4), vdw(5), vdw_soft(6),
    //      crf(7), crf_soft(8), distanceres(9), distancefield(10), t
    //      dihedralres(11), mass(12) 
    // kin is not there, I put it in position 0
    std::string tmp[] = {"kin", "bond", "ang", "dih", "impr", "lj", "slj", "crf", "scrf", "disres", "dihres", "disfld"}; // Betty
    std::vector<std::string> lambdas_types (tmp, tmp+(sizeof(tmp)/sizeof(std::string)));
    
    for (unsigned int i=0; i<lambdas_types.size(); i++) {
       // ref: by default 0 0 0 1 0
       vector<double> coeffs_sim = args.getValues<double>("lam"+lambdas_types[i]+"_sim",5,false, Arguments::Default<double>() << 0 << 0 << 0 << 1 << 0);
       abcde_sim[lambdas_types[i]] = coeffs_sim;
    }

    // the following values can be given as flags or read from an imd file
    // the parameters from the imd file are overwritten if the corresponding 
    // flag is also given
    double temp, minlam, maxlam, slam, NLAMs;
    int nrlambdas;
    string s_input;
    if (args.count("imd")>=0) {
      s_input = args.getValue<string>("imd", true, ""); // I have to give a default here, otherwise I get: std::logic_error
                                                        // what():  basic_string::_S_construct null not valid
      if (file_exists(s_input)) {
        Ginstream imd(s_input);
        input gin;
        imd >> gin;
        if (gin.perturbation.found) {
          slam=gin.perturbation.rlam;
          NLAMs=gin.perturbation.nlam;
        }
        if (gin.multibath.found) {
          temp=gin.multibath.temp0[0];
        }
        if (gin.precalclam.found) {
          minlam=gin.precalclam.minlam;
          maxlam=gin.precalclam.maxlam;
          nrlambdas=gin.precalclam.nrlam;
        }
        if (gin.lambdas.found) {    
          // read coefficients from imd LAMBDAS block    
          std::map<std::string, vector<double> > lambdas_map = read_lambdas(gin.lambdas.lambints, lambdas_types);

          // use those coefficients only if the corresponding argument flag does not exist
          for (std::map<std::string, vector<double> >::const_iterator it=lambdas_map.begin(), to=lambdas_map.end();
                it != to; it++) {
            if (args.count("lam"+it->first+"_sim") < 0)
              abcde_sim[it->first] = args.getValues<double>("lam"+it->first+"_sim",5,false, it->second);
          }
        }
        
        cerr << "# Read parameters from: " << s_input << endl;
        
        // if the following flags are given, the corresponding values read from 
        //the imd will be overwritten
        string tmp[] = {"slam", "NLAMs", "temp", "nrlambdas", "minlam", "maxlam"};
        std::vector<string> imdflags (tmp, tmp+sizeof(tmp)/sizeof(string));
        ostringstream oss;
        bool warn=false;
        oss << "# WARNING: ";
        for (std::vector<string>::const_iterator it=imdflags.begin(), to=imdflags.end(); it != to; it++) {
          if (args.count(*it)>=0) {
            oss << " @" << *it; 
            warn=true;
          }
        }
        oss << " flag(s) overwrite values from the imd file\n";
        if (warn) cerr << oss.str();
        
        // overwrite if parameters were given as input flag
        temp = args.getValue<double>("temp", false, temp);
        nrlambdas = args.getValue<int>("nrlambdas", false, nrlambdas);
        minlam = args.getValue<double>("minlam", false, minlam);
        maxlam = args.getValue<double>("maxlam", false, maxlam);
        slam = args.getValue<double>("slam", false, slam);
        NLAMs = args.getValue<double>("NLAMs", false, NLAMs);   
      }
      else throw gromos::Exception("ext_ti_ana", "@imd: file does not exist: "+s_input+"\n");
    } else {
      temp = args.getValue<double>("temp", true);
      nrlambdas = args.getValue<int>("nrlambdas", true);
      minlam = args.getValue<double>("minlam", true);
      maxlam = args.getValue<double>("maxlam", true);
      slam = args.getValue<double>("slam", true);
      NLAMs = args.getValue<double>("NLAMs", true);  
    }
   
    double  NLAMp = args.getValue<double>("NLAMp", false, NLAMs); 
    double lam_step = (maxlam - minlam)/(nrlambdas-1);

    // fill vector of lambdas for which we have pre-calculated terms
    vector<double> p_LAM(1); 
    p_LAM.resize(nrlambdas);
    for(int p=0; p<nrlambdas; p++){
         p_LAM[p] = p *lam_step + minlam;
    }    

    // now that we know what p_LAM values we have, we can check for any bar_data and dHdl_data
    bool bar_data = false;
    bool dHdl_data = false;
    std::vector<double> bar_lam;
    std::map<int, int> bar_lam_ind;
    std::vector<double> dHdl_lam;
    std::map<int, int> dHdl_lam_ind;

    if(args.count("bar_data")>=0){
      bar_data = true;
      if(args.count("bar_data")>0){
	bar_lam = args.getValues<double>("bar_data", args.count("bar_data"),false);
	std::sort(bar_lam.begin(), bar_lam.end());
	
	for(unsigned int i=0; i<bar_lam.size(); i++){
	  for(unsigned int p = 0; p < p_LAM.size(); p++){
	    if(abs(bar_lam[i] - p_LAM[p]) < lam_epsilon){
	      bar_lam_ind[p] = i;
	    }
	  }
	}
	if(bar_lam.size() != bar_lam_ind.size())
	  throw gromos::Exception("ext_ti_ana", "not all bar_data values are available - check precision?");
      }
    }

    if(args.count("dHdl_data")>=0){
      dHdl_data = true;
      if(args.count("dHdl_data")>0){
	dHdl_lam = args.getValues<double>("dHdl_data", args.count("dHdl_data"),false);
	std::sort(dHdl_lam.begin(), dHdl_lam.end());

	for(unsigned int i=0; i<dHdl_lam.size(); i++){
	  for(unsigned int p = 0; p < p_LAM.size(); p++){
	    if(abs(dHdl_lam[i] - p_LAM[p]) < lam_epsilon){
	      dHdl_lam_ind[p] = i;
	    }
	  }
	}
	if(dHdl_lam.size() != dHdl_lam_ind.size())
	  throw gromos::Exception("ext_ti_ana", "not all dHdl_data values are available - check precision?");
      }
    }

    // prepare Lambdas
    std::map<std::string, LambdaStruct> reflambdas;
    std::map<std::string, std::vector<LambdaStruct> > lambdas;
    for (unsigned int i=0; i<lambdas_types.size(); i++) {
       std::string lambda_type = lambdas_types[i];
       reflambdas[lambda_type] = getSingleLambda(abcde_sim[lambda_type], minlam, lam_step, nrlambdas, slam, p_LAM, lambda_type);
       // if no flag given, by default the lambda coefficients for 
       // prediction are the same as the ref ones
       abcde[lambda_type] = args.getValues<double>("lam"+lambda_type,5,false, abcde_sim[lambda_type]);
       lambdas[lambda_type] = getLambdas(abcde[lambda_type], minlam, lam_step, nrlambdas, p_LAM, lambda_type);
    } 

    // for slj and scrf we can read multiple sets from a file: 
    vector<vector<LambdaStruct> > SLJ, SCRF;
    vector<vector<double> > SLJ_abcde, SCRF_abcde;
    vector<vector<int> > slabels;
    if(slj_scrf_file){
      cerr << "# slj_scrf_file specified" << endl;
      if(args.count("lamslj")>=0 || args.count("lamscrf")>=0){
        throw gromos::Exception("ext_ti_ana", "specify slj_scrf_file OR lamslj/lamscrf"); 
      }
      string slj_scrf_filename = args["slj_scrf_file"];
      read_slj_scrf_file(slj_scrf_filename, SLJ, SCRF, slabels, minlam, lam_step, nrlambdas, p_LAM, SLJ_abcde, SCRF_abcde) ;

    } else {
      cerr << "# single slj and scrf value" << endl;
      //individual lambdas: slj
      SLJ.push_back(getLambdas(abcde["slj"], minlam, lam_step, nrlambdas, p_LAM, "slj") );
      SLJ_abcde.push_back(abcde["slj"]);
      //individual lambdas: scrf
      SCRF.push_back(getLambdas(abcde["scrf"], minlam, lam_step, nrlambdas, p_LAM, "scrf") );
      SCRF_abcde.push_back(abcde["scrf"]);
    }
    
    //cerr << "finished reading lambdas" << endl;

    LambdaStruct Rslj = reflambdas["slj"];
    LambdaStruct Rscrf = reflambdas["scrf"];
    LambdaStruct Rlj = reflambdas["lj"];
    LambdaStruct Rcrf = reflambdas["crf"];
    LambdaStruct Rkin = reflambdas["kin"];
    LambdaStruct Rbon = reflambdas["bond"];
    LambdaStruct Rang = reflambdas["ang"];
    LambdaStruct Rdih = reflambdas["dih"];
    LambdaStruct Rimp = reflambdas["impr"];
    /** Betty
     * perturbed special ineractions
     */
    LambdaStruct Rdisres = reflambdas["disres"];
    LambdaStruct Rdihres = reflambdas["dihres"];
    LambdaStruct Rdisfld = reflambdas["disfld"];
    // Betty

    //Lambdas slj = lambdas["slj"];
    //Lambdas scrf = lambdas["scrf"];
    vector<LambdaStruct> lj = lambdas["lj"];
    vector<LambdaStruct> crf = lambdas["crf"];
    vector<LambdaStruct> kin = lambdas["kin"];
    vector<LambdaStruct> bon = lambdas["bond"];
    vector<LambdaStruct> ang = lambdas["ang"];
    vector<LambdaStruct> dih = lambdas["dih"];
    vector<LambdaStruct> imp = lambdas["impr"];
    /** Betty
     * perturbed special ineractions
     */
    vector<LambdaStruct> disres = lambdas["disres"];
    vector<LambdaStruct> dihres = lambdas["dihres"];
    vector<LambdaStruct> disfld = lambdas["disfld"];
    // Betty

    // determine the plam indeces to be used
    int pmin = args.getValue<int>("pmin", false, 0);
    int pmax = args.getValue<int>("pmax", false, p_LAM.size());
    int nr_plam = pmax-pmin; 
  
    int num_scrf=SCRF.size();
    int num_slj=SLJ.size();
    
    // declare some vectors which we will need later
    vector<vector<vector<gmath::Stat<double> > > > x( num_slj, vector<vector<gmath::Stat<double> > > (num_scrf, vector<gmath::Stat<double> > (nr_plam)));
    vector<vector<vector<gmath::Stat<double> > > > vyvr( num_slj, vector<vector<gmath::Stat<double> > > (num_scrf, vector<gmath::Stat<double> > (nr_plam)));
    vector<vector<vector<gmath::Stat<double> > > > expvyvr( num_slj, vector<vector<gmath::Stat<double> > > (num_scrf, vector<gmath::Stat<double> > (nr_plam)));
    vector<vector<vector<gmath::Stat<double> > > > Xexpvyvr( num_slj, vector<vector<gmath::Stat<double> > > (num_scrf, vector<gmath::Stat<double> > (nr_plam)));

    // for BAR we need the reference energy and the energies at the other lambda values
    int nr_bar_lam=0;
    if(bar_data){
      if(bar_lam.size()==0) nr_bar_lam = nr_plam;
      else nr_bar_lam = bar_lam.size();
    }
    gmath::Stat<double> E_sim;
    vector<vector<vector<gmath::Stat<double> > > > E_bar( num_slj, vector<vector<gmath::Stat<double> > > (num_scrf, vector<gmath::Stat<double> > (nr_bar_lam)));
    
    // DHDL
    int nr_dHdl_lam=0;
    if(dHdl_data){
      if(dHdl_lam.size()==0) nr_dHdl_lam = nr_plam;
      else nr_dHdl_lam = dHdl_lam.size();
    }
    vector<vector<vector<gmath::Stat<double> > > > E_dHdl( num_slj, vector<vector<gmath::Stat<double> > > (num_scrf, vector<gmath::Stat<double> > (nr_dHdl_lam)));    

    double kBT = gmath::physConst.get_boltzmann() * temp;
    
    // define an energy trajectory
    utils::EnergyTraj etrj;

    // read topology for the mass and the number of molecules
    double mass=0;
    etrj.addConstant("MASS",mass);
    
    // learn about the variable names how they map to the elements
    read_library(library, etrj);
    
    // define two input streams
    Ginstream gin_en;
    Ginstream gin_fr;
    bool do_energy_files     =(args.count("en_files")>0);
    bool do_free_energy_files=(args.count("fr_files")>0);
    
    Arguments::const_iterator it_en=args.lower_bound("en_files"),
      to_en=args.upper_bound("en_files"),
      it_fr=args.lower_bound("fr_files"),
      to_fr=args.upper_bound("fr_files");


    //int cont=0, en_cont=0, fr_cont=0;
    if(do_energy_files) {
      gin_en.open(it_en->second.c_str()); 
      ++it_en; 
      //en_cont=1;
    }
    
    if(do_free_energy_files) {
      gin_fr.open(it_fr->second.c_str());
      ++it_fr;
      //fr_cont=1;
    }

    //cont=en_cont+fr_cont;
    bool version_checked = false;
    
    // starting with the calculations    
    while (true) {
      // version number
      if (!version_checked) {
        version_checked = check_etraj_version(etrj, gin_en, gin_fr, do_energy_files, do_free_energy_files, library);
      }
    
      // read the next frame of the energy trajectory
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
      // read the next frame of the free energy trajectory
      if(do_free_energy_files){
	int end_of_file=etrj.read_frame(gin_fr, "FRENERTRJ");
	if(end_of_file){
	  if(it_fr!=to_fr){
	    gin_fr.close();
	    gin_fr.open(it_fr->second.c_str());
	    ++it_fr;
	    //try again...
	    etrj.read_frame(gin_fr, "FRENERTRJ");
	  }
	  else
	    break;
	}
      }

      int eblockindex = etrj.return_blockindex("PRECALCLAM");
      int feblockindex = etrj.return_blockindex("FREEPRECALCLAM");
      vector<vector<double> > se;
      vector<vector<double> > sfe;
      int eDihblockindex = etrj.return_blockindex("ABDIH");
      int feDihblockindex = etrj.return_blockindex("FREEABDIH");
      vector<double> seDih;
      vector<double> sfeDih;

      // read data
      se = etrj.return_block(eblockindex);
      sfe = etrj.return_block(feblockindex);
      seDih = etrj.return_line(eDihblockindex,0);
      sfeDih = etrj.return_line(feDihblockindex,0);

      // se[l][0] : A_e_lj
      // se[l][1] : B_e_lj
      // se[l][2] : A_e_crf
      // se[l][3] : B_e_crf
      // se[l][4] : AB_e_bond
      // se[l][5] : AB_e_angle
      // se[l][6] : AB_e_improper
      // sfe[l][0] : A_de_lj
      // sfe[l][1] : B_de_lj
      // sfe[l][2] : A_de_crf
      // sfe[l][3] : B_de_crf
      // sfe[l][4] : AB_de_bond
      // sfe[l][5] : AB_de_angle
      // sfe[l][6] : AB_de_improper

      // etrj has a different number of precalculated values than
      // we requested
      if (se.size() != p_LAM.size()) {
          ostringstream os;
          os << "Expected "<< p_LAM.size() << " precalculated values but found "<< se.size() << " in the energy trajectory\n"; 
          throw gromos::Exception("ext_ti_ana", os.str());
      }
      if (sfe.size() != p_LAM.size())  {
          ostringstream os;
          os << "Expected "<< p_LAM.size() << " precalculated values but found "<< se.size() << " in the free energy trajectory\n"; 
          throw gromos::Exception("ext_ti_ana", os.str());
      }
      
      // calculate the reference energies
      double pow1_Rlj = pow((1-Rlj.lam),NLAMs);
      double pow_Rlj = pow(Rlj.lam,NLAMs);
      double pow1_Rcrf = pow((1-Rcrf.lam),NLAMs);
      double pow_Rcrf = pow(Rcrf.lam,NLAMs);
      double E_s = 0.0;

      E_s += pow1_Rlj * (Rslj.wlow*se[Rslj.idx_low][0] + Rslj.whigh*se[Rslj.idx_high][0]) 
           + pow_Rlj  * (Rslj.wlow*se[Rslj.idx_low][1] + Rslj.whigh*se[Rslj.idx_high][1]); //lj
      E_s += pow1_Rcrf * (Rscrf.wlow*se[Rscrf.idx_low][2] + Rscrf.whigh*se[Rscrf.idx_high][2])
           + pow_Rcrf  * (Rscrf.wlow*se[Rscrf.idx_low][3] + Rscrf.whigh*se[Rscrf.idx_high][3]); //crf
      E_s += (Rkin.wlow*se[Rkin.idx_low][4] + Rkin.whigh*se[Rkin.idx_high][4]); //kin
      E_s += (Rbon.wlow*se[Rbon.idx_low][5] + Rbon.whigh*se[Rbon.idx_high][5]); //bonds
      E_s += (Rang.wlow*se[Rang.idx_low][6] + Rang.whigh*se[Rang.idx_high][6]); //angle
      E_s += (Rimp.wlow*se[Rimp.idx_low][7] + Rimp.whigh*se[Rimp.idx_high][7]); //impdih
      /** Betty
       * perturbed special interactions
       */
      E_s += (Rdisres.wlow*se[Rdisres.idx_low][8] + Rdisres.whigh*se[Rdisres.idx_high][8]);
      E_s += (Rdihres.wlow*se[Rdihres.idx_low][9] + Rdihres.whigh*se[Rdihres.idx_high][9]);
      E_s += (Rdisfld.wlow*se[Rdisfld.idx_low][10] + Rdisfld.whigh*se[Rdisfld.idx_high][10]);
      //Betty
      E_s += (1-Rdih.lam) * seDih[0] + Rdih.lam * seDih[1]; //dih
      
      if(bar_data) E_sim.addval(E_s);

      // for all plam, calculate the energies
      for(int p=pmin; p<pmax; p++){

        // precalculate some constants
        double pow1_lj = pow((1-lj[p].lam),NLAMp);
        double pow_lj = pow(lj[p].lam,NLAMp);
        double pow1_crf = pow((1-crf[p].lam),NLAMp);
        double pow_crf = pow(crf[p].lam,NLAMp);

        // calculate energies (with interpolation between lambda_i)
        double E = 0.0;
        E += (kin[p].wlow*se[kin[p].idx_low][4] + kin[p].whigh*se[kin[p].idx_high][4]); //kin
        E += (bon[p].wlow*se[bon[p].idx_low][5] + bon[p].whigh*se[bon[p].idx_high][5]); //bonds
        E += (ang[p].wlow*se[ang[p].idx_low][6] + ang[p].whigh*se[ang[p].idx_high][6]); //angle
        E += (imp[p].wlow* se[imp[p].idx_low][7] + imp[p].whigh*se[imp[p].idx_high][7]); //impdih
        /** Betty
	 * perturbed special interactions
	 */
	E += (disres[p].wlow* se[disres[p].idx_low][8] + disres[p].whigh*se[disres[p].idx_high][8]);
	E += (dihres[p].wlow* se[dihres[p].idx_low][9] + dihres[p].whigh*se[dihres[p].idx_high][9]);
	E += (disfld[p].wlow* se[disfld[p].idx_low][10] + disfld[p].whigh*se[disfld[p].idx_high][10]);
	// Betty
        E += (1-dih[p].lam) * seDih[0] + dih[p].lam * seDih[1]; //dih

        // and their derivatives
        double dE = 0.0;
        double dpow1_lj = NLAMp * pow((1-lj[p].lam),(NLAMp-1));
        double dpow_lj = NLAMp * pow(lj[p].lam,(NLAMp-1));
        double dpow1_crf = NLAMp * pow((1-crf[p].lam),(NLAMp-1));
        double dpow_crf = NLAMp * pow(crf[p].lam,(NLAMp-1));

        if(!no_kin)     dE += (kin[p].wlow*sfe[kin[p].idx_low][4] + kin[p].whigh*sfe[kin[p].idx_high][4]) *kin[p].dlam;
        if(!no_bond)   dE += (bon[p].wlow*sfe[bon[p].idx_low][5] + bon[p].whigh*sfe[bon[p].idx_high][5]) *bon[p].dlam;
        if(!no_ang)   dE += (ang[p].wlow*sfe[ang[p].idx_low][6] + ang[p].whigh*sfe[ang[p].idx_high][6]) *ang[p].dlam;
        if(!no_impr)  dE += (imp[p].wlow*sfe[imp[p].idx_low][7] + imp[p].whigh*sfe[imp[p].idx_high][7]) *imp[p].dlam;
        /** Betty
	 * perturbed special interactions
	 */
	if(!no_disres)  dE += (disres[p].wlow*sfe[disres[p].idx_low][8] + disres[p].whigh*sfe[disres[p].idx_high][8]) *disres[p].dlam;
	if(!no_dihres)  dE += (dihres[p].wlow*sfe[dihres[p].idx_low][9] + dihres[p].whigh*sfe[dihres[p].idx_high][9]) *dihres[p].dlam;
	if(!no_disfld)  dE += (disfld[p].wlow*sfe[disfld[p].idx_low][10] + disfld[p].whigh*sfe[disfld[p].idx_high][10]) *disfld[p].dlam;
	// Betty
        if(!no_dih) dE += sfeDih[0]; 

        //calculate lj and crf contributions for all slj and scrf parameterizations
        double Emult;
        double dEmult;
        for(unsigned int i=0; i<SLJ.size(); i++){
          for(unsigned int j=0; j<SCRF.size(); j++){
          // calculate energies (with interpolation between lambda_i)
          Emult = E;
          Emult += pow1_lj * (SLJ[i][p].wlow*se[SLJ[i][p].idx_low][0] + SLJ[i][p].whigh*se[SLJ[i][p].idx_high][0])
             + pow_lj  * (SLJ[i][p].wlow*se[SLJ[i][p].idx_low][1] + SLJ[i][p].whigh*se[SLJ[i][p].idx_high][1]); //lj
          Emult += pow1_crf * (SCRF[j][p].wlow*se[SCRF[j][p].idx_low][2] + SCRF[j][p].whigh*se[SCRF[j][p].idx_high][2])
             + pow_crf  * (SCRF[j][p].wlow*se[SCRF[j][p].idx_low][3] + SCRF[j][p].whigh*se[SCRF[j][p].idx_high][3]); //crf

          dEmult = dE;
          if(!no_lj) dEmult += (-dpow1_lj * (SLJ[i][p].wlow*se[SLJ[i][p].idx_low][0] + SLJ[i][p].whigh*se[SLJ[i][p].idx_high][0]) 
                            + dpow_lj * (SLJ[i][p].wlow*se[SLJ[i][p].idx_low][1] + SLJ[i][p].whigh*se[SLJ[i][p].idx_high][1]) ) * lj[p].dlam 
                        + (pow1_lj * (SLJ[i][p].wlow*sfe[SLJ[i][p].idx_low][0] + SLJ[i][p].whigh*sfe[SLJ[i][p].idx_high][0]) 
                          - pow_lj * (SLJ[i][p].wlow*sfe[SLJ[i][p].idx_low][1] + SLJ[i][p].whigh*sfe[SLJ[i][p].idx_high][1]) ) * SLJ[i][p].dlam;
          if(!no_crf) dEmult += (-dpow1_crf * (SCRF[j][p].wlow*se[SCRF[j][p].idx_low][2] + SCRF[j][p].whigh*se[SCRF[j][p].idx_high][2]) 
                            + dpow_crf * (SCRF[j][p].wlow*se[SCRF[j][p].idx_low][3] + SCRF[j][p].whigh*se[SCRF[j][p].idx_high][3]) ) * crf[p].dlam 
                        + (pow1_crf * (SCRF[j][p].wlow*sfe[SCRF[j][p].idx_low][2] + SCRF[j][p].whigh*sfe[SCRF[j][p].idx_high][2]) 
                          - pow_crf * (SCRF[j][p].wlow*sfe[SCRF[j][p].idx_low][3] + SCRF[j][p].whigh*sfe[SCRF[j][p].idx_high][3]) ) * SCRF[j][p].dlam;

          // prepare for reweighting
          double diff = -(Emult-E_s)/(gmath::physConst.get_boltzmann() * temp);
          x[i][j][p-pmin].addval(dEmult);
          double expdiff = exp(diff);
          vyvr[i][j][p-pmin].addval(diff);
          expvyvr[i][j][p-pmin].addval(expdiff);
          Xexpvyvr[i][j][p-pmin].addval(dEmult * expdiff);

	  if(bar_data) {
	    if(bar_lam.size()==0) E_bar[i][j][p-pmin].addval(Emult);
	    else{
	      if(bar_lam_ind.count(p-pmin))
		E_bar[i][j][bar_lam_ind[p-pmin]].addval(Emult);
	    }
	    
	  }
	  if(dHdl_data) {
	    if(dHdl_lam.size()==0) E_dHdl[i][j][p-pmin].addval(dEmult);
	    else{
	      if(dHdl_lam_ind.count(p-pmin))
		E_dHdl[i][j][dHdl_lam_ind[p-pmin]].addval(dEmult);
	    }

	  }
	  
          }
        }//i
      }//p
    }//frames    
        
    #ifdef OMP
        if (num_cpus > SLJ.size()){
            if(int(SLJ.size()) > omp_get_max_threads())
                num_cpus = omp_get_max_threads();
            else
                num_cpus = SLJ.size();
            cerr << "# Number of threads > number of trajectory files: not feasible. Corrected to " 
                 << num_cpus << " threads." << endl;
        }

        if ((int)num_cpus > omp_get_max_threads()){
            cerr << "# You specified " << num_cpus << " number of threads. There are only " 
                 << omp_get_max_threads() << " threads available." << endl;
            num_cpus = omp_get_max_threads();
        } 
        omp_set_num_threads(num_cpus); //set the number of cpus for the parallel section
        if (num_cpus > 1)  cerr << "# Number of threads: " << num_cpus << endl;
    #else
        if(num_cpus != 1)
            throw gromos::Exception("hbond","@cpus: Your compilation does not support multiple threads. Use --enable-openmp for compilation.\n\n" + usage);
    #endif

    //initialize random number generator for bootstraps
    int ti = time(NULL);
    #ifdef OMP
    double totaltime=0;
    double bootstrap_time=0;
    #pragma omp parallel for reduction(+:totaltime,bootstrap_time) //firstprivate(sys, time) 
    #endif
    // loop over all trajectories
    for(unsigned int i=0; i<SLJ.size(); i++){ 
        //#ifdef OMP
        //#pragma omp critical
        //#endif
        //cerr << "# SLJ: " << i << endl;
        #ifdef OMP
        double start_tot = omp_get_wtime();
        #endif

        unsigned int seed=ti+i;
      
        abcde["slj"] = SLJ_abcde[i];
        for(unsigned int j=0; j<SCRF.size(); j++){
            abcde["scrf"] = SCRF_abcde[j];

            ostringstream oss;
            oss << "#SLAM " << slam << endl;  
            if (verbose) {
                oss << "# NLAMs=" << NLAMs
                    << ", NLAMp=" << NLAMp
                    << ", T=" << temp << endl;
            
                oss << "#----- lambda coefficients: ------------\n";
                oss << "#        coeffs_sim    coeffs_prediction" << endl;
                for (unsigned int ii=0; ii<lambdas_types.size(); ii++) {
                    oss << "# " << setw(4) <<lambdas_types[ii]<< ": "
                        << abcde_sim[lambdas_types[ii]] << " | " 
                        << abcde[lambdas_types[ii]] << std::endl;
                }
                oss << "#---------------------------------------\n";
            } else {
                // if we read multiple slj-scrf combinations from a file 
                // it will be good to always have the used combination in
                // the header of the output
                oss << "# " << setw(4) << "slj" << ": "<< abcde_sim["slj"] 
                    << " | " << abcde["slj"] << std::endl;
                oss << "# " << setw(4) << "scrf" << ": "<< abcde_sim["scrf"] 
                    << " | " << abcde["scrf"] << std::endl;
            }
            oss << "#plam "<< setw(15)<<"pred_dHdl";
            if (!no_error) oss << "" << setw(15) << "err";
            if (countframes) oss << " " << setw(15) << "contrib.frames";
            oss << "\n";
      
            // for the file with the bar data
	    ofstream oss_bar;
	    if(bar_data){
		
	      ostringstream os_bar;
	      if (slabels.size() == 0) {
		// we only have a single slj/scrf combination
		os_bar << outdir << "/bar_data_l" << setprecision(lam_precision) << fixed << slam << ".dat";
	      } else {
		os_bar << outdir << "/bar_data_l" << setprecision(lam_precision) << fixed << slam << "_slj" << slabels[0][i] << "_scrf" << slabels[1][j] << ".dat";
	      }
	      oss_bar.open(os_bar.str().c_str());

	      oss_bar << "#      nr lam_p" << endl;
	      oss_bar << setw(15) << nr_bar_lam << endl;
	      oss_bar << "#         lam_s" << setw(19) << "lam_p ..." << endl;
	      oss_bar << setw(15) << setprecision(lam_precision) << fixed << slam;

	      for(int p=0; p<nr_bar_lam; p++){
		if(bar_lam.size() ==0)
		  oss_bar << setw(15) << setprecision(lam_precision) << fixed << p_LAM[p+pmin];
		else
		  oss_bar << setw(15) << setprecision(lam_precision) << fixed << bar_lam[p];
	      }
	      oss_bar << endl;
	      oss_bar << "#           E_s" << setw(19) << "E_p(lam_p) ..." << endl;	      
	      for(int n = 0; n < E_sim.n(); n++){
		
		oss_bar << setw(15) << setprecision(8) << E_sim.val(n);
		oss_bar << ' ';

		for(int p=0; p<nr_bar_lam; p++){
		  oss_bar << setw(15) << setprecision(8) << E_bar[i][j][p].val(n);
		  oss_bar << ' ';
		}
		oss_bar << endl;
	      }
	    }

            // for the file with the dHdl data
	    ofstream oss_dHdl;
	    if(dHdl_data){

	      ostringstream os_dHdl;
	      if (slabels.size() == 0) {
		// we only have a single slj/scrf combination
		os_dHdl << outdir << "/dHdl_data_l" << setprecision(lam_precision) << fixed << slam << ".dat";
	      } else {
		os_dHdl << outdir << "/dHdl_data_l" << setprecision(lam_precision) << fixed << slam << "_slj" << slabels[0][i] << "_scrf" << slabels[1][j] << ".dat";
	      }
	      oss_dHdl.open(os_dHdl.str().c_str());

	      oss_dHdl << "# lam_s " << setw(15) << setprecision(lam_precision) << fixed << slam << endl;
	      oss_dHdl << "#      nr lam_p" << endl;
	      oss_dHdl << setw(15) << nr_dHdl_lam << endl;
	      oss_dHdl << "#         lam_p" << setw(19) << endl;

	      for(int p=0; p<nr_dHdl_lam; p++){
		if(dHdl_lam.size() ==0)
		  oss_dHdl << setw(15) << setprecision(lam_precision) << fixed << p_LAM[p+pmin];
		else
		  oss_dHdl << setw(15) << setprecision(lam_precision) << fixed << dHdl_lam[p];
	      }

	      oss_dHdl << endl;
	      oss_dHdl << "#          dE_p ..." << setw(19) << endl;

	      for(int n = 0; n < E_sim.n(); n++){
		for(int p=0; p<nr_dHdl_lam; p++){
		  oss_dHdl << setw(15) << setprecision(8) << E_dHdl[i][j][p].val(n);
		  oss_dHdl << ' ';
		}
		oss_dHdl << endl;
	      }
	    }
	    
            // for each plam
            for(int p=0; p<nr_plam; p++){
      
                oss << setprecision(lam_precision) << fixed << p_LAM[p+pmin] ;
                oss.precision(6);
        
                // calculate the final value
                int sign = 0;
                double lnXexpave = gmath::Stat<double>::lnXexpave(x[i][j][p], vyvr[i][j][p], sign);
                double lnexpave = vyvr[i][j][p].lnexpave();
                double final = exp(lnXexpave - lnexpave) * sign;
                oss << " " << setw(15) << final;

                // calculate statistical uncertainty if required
                if(bootstrap) {
                    #ifdef OMP
                    double start_int=omp_get_wtime();
                    #endif
                    int length = x[i][j][p].n();
                    gmath::Stat<double> final_boot;
                    for(int ii=0;ii<bootstrap;ii++){
                        gmath::Stat<double> x_new;
                        gmath::Stat<double> vyvr_new;

                        for(int jj=0;jj<length;jj++){
                            double r;
                            #ifdef OMP
                            #pragma omp critical
                            #endif
                            r=rand_r(&seed);
                            int index = int((r*length)/RAND_MAX);
                            x_new.addval(x[i][j][p].val(index));
                            vyvr_new.addval(vyvr[i][j][p].val(index));
                        }
                        double lnXexpave = gmath::Stat<double>::lnXexpave(x_new, vyvr_new, sign);
                        double lnexpave = vyvr_new.lnexpave();
                        final_boot.addval(exp(lnXexpave - lnexpave) * sign);
                    }
                    double boot_std = sqrt(gmath::Stat<double>::covariance(final_boot,final_boot));
                    oss << " " << setw(15) << boot_std;
                    
                    #ifdef OMP
                    bootstrap_time += omp_get_wtime()-start_int;
                    #endif
                }
            
                // count nr of contributing frames if required
                if(countframes){
                    int count = 0;
                    double deltag = -kBT * lnexpave;
                    for(int ii=0; ii<vyvr[i][j][p].n(); ii++){
                        if(-kBT*vyvr[i][j][p].val(ii) <= (deltag+0.001)){
                            count = count + 1;
                        }
                    }
                    oss << " "<< setw(15) << count;
                }
            oss  << endl;
            }//p
            ostringstream os;
            if (slabels.size() == 0) {
                // we only have a single slj/scrf combination
                os << outdir << "/predicted_l" << setprecision(lam_precision) << fixed << slam << ".dat";
            } else {
                os << outdir << "/predicted_l" << setprecision(lam_precision) << fixed << slam << "_slj" << slabels[0][i] << "_scrf" << slabels[1][j] << ".dat";
            }
            ofstream fout(os.str().c_str());
            fout << oss.str();
            fout.close();
        }//j(SCRF)
        #ifdef OMP
        totaltime += omp_get_wtime()-start_tot;
        #endif
    }//i(SLJ)
    #ifdef OMP
    cout.precision(2);
    //cout << "# Preparation time: " << prep_time << " s" << endl;
    //cout << "# Total calc time: " << bootstrap_time << " s" << endl;
    cout << "# Total CPU time: \t" << totaltime << " s" << endl;
    cout << "### Total real time: \t"<< omp_get_wtime()-start_total << " s" << endl;
    #endif
  }
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

LambdaStruct getSingleLambda(vector<double> abcde, double minlam, double lam_step,
                    int nrlambdas, double slam, vector<double> &p_LAM, std::string lambda_type)
{
  LambdaStruct newlambdas;
  double lam = abcde[0]*pow(slam,4) + abcde[1]*pow(slam,3) + abcde[2]*pow(slam,2) + abcde[3]*slam + abcde[4];
  double dlam = 4*abcde[0]*pow(slam,3) + 3*abcde[1]*pow(slam,2) + 2*abcde[2]*slam + abcde[3];

  double maxlam = minlam + (nrlambdas-1)*lam_step;
  if (lambda_type != "lj" && lambda_type != "crf" && lambda_type != "dih") {
      if (lam > maxlam || lam < minlam) {
             ostringstream os;
             os << "at lambda=" << slam << ": with lambda coefficients " << abcde 
                << "\n            we get final lambda values outside the precalculated"
                << "\n            range [" << minlam << "-"<< maxlam
                <<"], interpolation not possible.";
             throw gromos::Exception("ext_ti_ana", os.str());
      }
  }

  newlambdas.idx = int(round((lam-minlam)/lam_step));

  int idx_low = floor((lam-minlam)/lam_step);
  int idx_high = ceil((lam-minlam)/lam_step);

  if( idx_low == idx_high){
    // lam is in our precalclam list
    newlambdas.wlow = 0;
    newlambdas.whigh= 1;
  }
  else{
    // not in the list, determine weights for interpolation
    newlambdas.wlow = (lam-p_LAM[idx_high])/(p_LAM[idx_low]-p_LAM[idx_high]);
    newlambdas.whigh = (lam-p_LAM[idx_low])/(p_LAM[idx_high]-p_LAM[idx_low]);
  }

//  cerr << idx_low << " " << idx_high << " " << newlambdas.wlow << " " << newlambdas.whigh << endl;

  newlambdas.idx_low = idx_low;
  newlambdas.idx_high = idx_high;
  newlambdas.lam = lam;
  newlambdas.dlam = dlam;

  return newlambdas;
}

vector <LambdaStruct> getLambdas(vector<double> abcde, double minlam, double lam_step,
                   int nrlambdas, vector<double> &p_LAM, std::string lambda_type)
{
  vector<LambdaStruct> newlambdas(p_LAM.size());  

  for(unsigned int p=0; p<p_LAM.size(); p++){
    newlambdas[p]=getSingleLambda(abcde, minlam, lam_step, nrlambdas, p_LAM[p], p_LAM, lambda_type);
  }

  return newlambdas;
}

void set_standards(utils::EnergyTraj &e)
{  
  e.addConstant("BOLTZ", gmath::physConst.get_boltzmann());
}

void read_slj_scrf_file(string name, vector<vector<LambdaStruct> >  &SLJ, 
                        vector<vector<LambdaStruct> >  &SCRF,  vector<vector<int> > &slabels, 
                        double minlam, double lam_step, int nrlambdas, vector<double> &p_LAM, 
                        vector<vector<double> > &SLJ_abcde, vector<vector<double> > &SCRF_abcde)
{
  Ginstream gin;
  vector<string> buffer;
  
  string slj_scrf_file_format;
  std::vector<int> slj_labels, scrf_labels;
  
  try{
    gin.open(name);
  }
  catch (gromos::Exception ex){
      cout << ex.what() << endl;
      throw gromos::Exception("read_slj_scrf_file", "failed to open slj_scrf file " +name);
  }

  while(!gin.stream().eof()){
    gin.getblock(buffer);
    if(gin.stream().eof()) break;
    if(buffer[buffer.size()-1].find("END")!=0)
      throw gromos::Exception("read_slj_scrf_file", "file " + gin.name() +
                              " is corrupted. No END in "+buffer[0]+
                              " block. Got\n" + buffer[buffer.size()-1]);

    if(buffer[0] == "SLJ"){
      //cerr << "# reading SLJ block" << endl;
      for(unsigned int i=1; i< buffer.size()-1; i++){
        istringstream linestream(buffer[i]);
        int nr;
        vector<double> abcde(5);
        linestream >> nr;
        for(int iii = 0; iii < 5; ++iii){ linestream >> abcde[iii]; }
          SLJ.push_back(getLambdas(abcde, minlam,lam_step, nrlambdas, p_LAM, "slj"));
          SLJ_abcde.push_back(abcde);
          slj_labels.push_back(nr);
      }
    } else if (buffer[0] == "SCRF"){
      for(unsigned int i=1; i< buffer.size()-1; i++){
        istringstream linestream(buffer[i]);
        int nr;
        vector<double> abcde(5);
        linestream >> nr;
        for(int iii = 0; iii < 5; ++iii){ linestream >> abcde[iii]; }
          SCRF.push_back(getLambdas(abcde, minlam,lam_step, nrlambdas, p_LAM, "scrf"));
          SCRF_abcde.push_back(abcde);
          scrf_labels.push_back(nr);
      }
    } else {
      throw gromos::Exception("read_slj_scrf_file", "read unknown block:"
                              +buffer[0]+
                              "\n known are: SLJ, SCRF\n");
        
    }
  }
  slabels.push_back(slj_labels);
  slabels.push_back(scrf_labels);
   
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
      throw gromos::Exception("ene_ana", "Library file " + gin.name() +
			      " is corrupted. No END in "+buffer[0]+
			      " block. Got\n"
			      + buffer[buffer.size()-1]);
    string sdum;
    
    if(buffer[0]=="ENERTRJ" || buffer[0]=="FRENERTRJ"){
      for(unsigned int i=1; i< buffer.size()-1; i++){
	e.addBlock(buffer[i], buffer[0]);
	
      }
    }
    
    if (buffer[0] == "ENEVERSION") {
      // This is a ENEVERSION block that was not put directly after the
      // TITLE block. This is no big deal for the library, it would
      // however lead to an error in a trajectory.
      // Shall we disallow this?

      // Just to make sure there are not two blocks
      if (gin.has_version() || e.has_version()) {
        throw gromos::Exception("ene_ana", "Library file " + gin.name() +
            " is corrupted. Two ENEVERSION blocks found.");
      }
      string ver;
      gio::concatenate(buffer.begin() + 1, buffer.end()-1, ver);
      ver.erase( std::remove_if( ver.begin(), ver.end(), ::isspace ), ver.end() );
      e.set_version(ver);
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
  if (gin.has_version()) {
    e.set_version(gin.version());
  }
}

bool file_exists(string file){
    ifstream infile(file.c_str()); //destructor closes the file
    return infile.is_open();
}

int dir_exists(const char *path)
{
    struct stat info;

    if(stat( path, &info ) != 0)
        // path does not exists yet
        return 0;
    else if(info.st_mode & S_IFDIR)
        // is a directory
        return 1;
    else
        // exists but not a directory
        return 2;
}

std::map<std::string, vector<double> > read_lambdas(std::vector<ilambdas::lambint> &lambints, std::vector<std::string>  &lambdas_types) {      
  std::cerr <<  "# WARNING: LAMBDAS block: if there is more than one set of coefficients\n"
            <<  "#      for the same property (different lambda dependence for different\n"
            <<  "#      energy groups), ext_ti_ana will use only the first set!\n";

  std::map<std::string, vector<double> > lambdas_map;
  bool warn=false;
  int counts[lambdas_types.size()];
  for (unsigned int t = 0; t < lambdas_types.size(); t++) {
    counts[t]=0;
  }
  
  for (unsigned int i=0; i < lambints.size(); i++) {
    double tmp[] = {lambints[i].ali, lambints[i].bli, lambints[i].cli, lambints[i].dli, lambints[i].eli};
    std::vector<double> coeffs (tmp, tmp + sizeof(tmp)/sizeof(double));

    // t=0 is the kinetic part, which is not defined in gromosXX LAMBDAS block
    for (unsigned int t = 1; t < lambdas_types.size(); t++) {
      if (lambints[i].ntli == int(t)) {
        counts[t]++;

        // check if the coefficients are the same in all occurrences of this lambdas type
        if (lambdas_map.find(lambdas_types[t]) != lambdas_map.end()) {
        
          if (coeffs != lambdas_map[lambdas_types[t]]) { 
              warn=true; 
              //throw gromos::Exception("ext_ti_ana", "LAMBDAS block: lam_"+lambdas_types[t] + " has different coefficients for different energy groups");
          }
          
        } else {
          lambdas_map[lambdas_types[t]] = coeffs;
        }
      }
    }
  }    
  if (warn)  std::cerr <<  "# WARNING: LAMBDAS block: different coefficients for different energy groups\n";
  
  return lambdas_map;
}

bool check_etraj_version(utils::EnergyTraj &etrj, Ginstream &gin_en, Ginstream &gin_fr, 
                          bool do_energy_files, bool do_free_energy_files, 
                          std::string library) {
        std::string lib_version="NONE", en_version="NONE", fr_version="NONE";
        if (etrj.has_version()) {
          lib_version=etrj.get_version();
          cerr << "         Library " << library << " version: "
               << lib_version << std::endl;
                   
          if (do_energy_files) {
            if (gin_en.has_version()) {
              en_version=gin_en.version();
              if (!etrj.version_match(gin_en.version())) {
                cerr << "WARNING: Version number mismatch!\n";
              } else {
                cerr << "MESSAGE: Version number check successful!\n";
              }
            } else {
              cerr << "WARNING: Version number missing!\n";
            }
              cerr << "         Energy Trajectory " << gin_en.name() << " version: "
                   << en_version << std::endl;
          }
          if (do_free_energy_files) {
            if (gin_fr.has_version()) {
              fr_version=gin_fr.version();
              if (!etrj.version_match(gin_fr.version())) {
                cerr << "WARNING: Version number mismatch!\n";
              } else {
                cerr << "MESSAGE: Version number check successful!\n";
              }
            } else {
              cerr << "WARNING: Version number missing!\n";
            }
              cerr << "         Energy Trajectory " << gin_fr.name() << " version: "
                   << fr_version << std::endl;
          }
        } else {
          cerr << "WARNING: Version number missing!\n"
               << "         Library " << library << " version: "
               << lib_version << std::endl;
          if (do_energy_files) {
            cerr << "         Energy Trajectory " << gin_en.name() << " version: ";
            if (gin_en.has_version()) {
              cerr << gin_en.version() << std::endl;
            } else {
              cerr << "NONE" << std::endl;
            }
          }
          if (do_free_energy_files) {
            cerr << "         Energy Trajectory " << gin_fr.name() << " version: ";
            if (gin_fr.has_version()) {
              cerr << gin_fr.version() << std::endl;
            } else {
              cerr << "NONE" << std::endl;
            }
          }
        }
        
        return true;
}
