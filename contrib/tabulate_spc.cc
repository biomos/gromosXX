/**
 * @file tabulate_spc.cc
 * creates a header file containing potential energy and force tables for SPC
 */

/**
 * @page programs Program Documentation
 *
 * @anchor tabulate_spc
 * @section tabulate_spc tabulate forces and energies for SPC
 * @author @ref ns
 * @date 21.10.2008
 *
 * Program tabulate_spc is used to create a table containing forces and energies
 * for SPC. The program calculates values of squared oxygen-oxygen distances in 
 * SPC water on a grid. For every grid point three tables (O-O, O-H, H-H), containing
 * the energies, forces and their linear interpolation gradient, are calculated
 * and written to the output.
 *
 * The number of grid points for short- and longrange interaction can be controlled
 * using the \@table argument. 5000 points for shortrange and 2000 points for 
 * longrange have been tested and give good accuracy. Especially the shortrange 
 * table size should be reduced to increase performance. This will slighlty 
 * increase the error for energies and forces.
 *
 * The tables are calculated on a range of @f$r^2@f$ that is defined by the solvent
 * diameter and the cutoffs. The latter can be modified using the \@cutoff argument.
 *
 * Reaction field data (permittivity and kappa) can be controlled using the 
 * \@reaction_field argument.
 *
 * The program creates a header file. You should use this header file to replace
 * spc_table.h to apply your settings.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@table</td><td>&lt;table size: shortrange, longrange&gt; </td></tr>
 * <tr><td> \@cutoff</td><td>&lt;cutoffs: shortrange, longrange&gt; </td></tr>
 * <tr><td> \@reaction_field</td><td>&lt;reaction field parameters: epsilon, kappa&gt; </td></tr>
 * </table>
 *
 */

#include <map>
#include <sstream>
#include <stdheader.h>

#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <simulation/simulation.h>
#include <configuration/configuration.h>

#include <interaction/interaction_types.h>
#include <math/periodicity.h>
#include <math/volume.h>

#include <interaction/nonbonded/pairlist/pairlist.h>
#include <interaction/nonbonded/interaction/storage.h>
#include <interaction/nonbonded/interaction/nonbonded_parameter.h>

#include <interaction/nonbonded/interaction/nonbonded_term.h>
#include <interaction/nonbonded/interaction/latticesum.h>
#include <configuration/mesh.h>
#include <interaction/nonbonded/interaction/nonbonded_innerloop.h>
#include <locale>

#include <io/argument.h>
#include <util/usage.h>
#include <map>

using namespace std;

int main(int argc, char** argv) {
  util::Known knowns;
  knowns << "table" << "cutoff" << "reaction_field";
  string usage = string("#") + argv[0] + "\n";
  usage += "\t@table          <shortrange_size longrange_size>\n";
  usage += "\t@cutoff         <shortrange longrange>\n";
  usage += "\t@reaction_field <epsilon kappa>\n";

  io::Argument args;
  if (args.parse(argc, argv, knowns)){
    std::cerr << usage << std::endl;
    return 1;
  }
  
  /////////////////////////////////////////////////
  // general parameters for reaction field and SPC
  /////////////////////////////////////////////////
  double shortrange_cutoff = 0.8;
  double longrange_cutoff = 1.4;
  {
    io::Argument::const_iterator it = args.lower_bound("cutoff"),
        to = args.upper_bound("cutoff");
    if (it != to) {
      if (args.count("cutoff") != 2) {
        cerr << "cutoff needs two numbers." << endl;
        return 1;
      }
      
      istringstream arg1(it->second);
      if (!(arg1 >> shortrange_cutoff)) {
        cerr << "Cutoff must be numeric.\n";
        return 1;
      }
      ++it;
      istringstream arg2(it->second);
      if (!(arg2 >> longrange_cutoff)) {
        cerr << "Cutoff must be numeric.\n";
        return 1;
      }
    }
  }
  
  const double solvent_diameter = 0.25; // distance of H-H + move withing 5 steps
  const double c12 = 2.634129E-6;
  const double c6 = 2.617346E-3;
  const double qOO = -0.82 * -0.82;
  const double qOH = -0.82 * 0.41;
  const double qHH = 0.41 * 0.41;
  double solvent_permittivity = 61.0;
  double solvent_kappa = 0.0;
  {
    io::Argument::const_iterator it = args.lower_bound("reaction_field"),
        to = args.upper_bound("reaction_field");
    if (it != to) {
      if (args.count("reaction_field") != 2) {
        cerr << "reaction_field needs two numbers." << endl;
        return 1;
      }
      
      istringstream arg1(it->second);
      if (!(arg1 >> solvent_permittivity)) {
        cerr << "epsilon must be numeric.\n";
        return 1;
      }
      ++it;
      istringstream arg2(it->second);
      if (!(arg2 >> solvent_kappa)) {
        cerr << "kappa must be numeric.\n";
        return 1;
      }
    }
  }
  math::four_pi_eps_i = 138.9354;
  
  unsigned int shortrange_table_size = 5000;
  unsigned int longrange_table_size = 2000;
  {
    io::Argument::const_iterator it = args.lower_bound("table"),
        to = args.upper_bound("table");
    if (it != to) {
      if (args.count("table") != 2) {
        cerr << "table needs two numbers." << endl;
        return 1;
      }
      
      istringstream arg1(it->second);
      if (!(arg1 >> shortrange_table_size)) {
        cerr << "table size must be numeric.\n";
        return 1;
      }
      ++it;
      istringstream arg2(it->second);
      if (!(arg2 >> longrange_table_size)) {
        cerr << "table size must be numeric.\n";
        return 1;
      }
    }
  }
  // write the parameters to the headerfile for checking
  cout.precision(13);
  cout.flags(ios::scientific);
  cout << "#ifdef CHECK_PARAM" << endl
       << "const unsigned int shortrange_table = " << shortrange_table_size << ";" << endl
       << "const unsigned int longrange_table = " << longrange_table_size << ";" << endl
       << "const double shortrange_cutoff = " << setw(20) << shortrange_cutoff << ";" << endl
       << "const double longrange_cutoff = " << setw(20) << longrange_cutoff << ";" << endl
       << "const double solvent_diameter = " << setw(20) << solvent_diameter << ";" << endl
       << "const double c12 = " << setw(20) << c12 << ";" << endl
       << "const double c6 = " << setw(20) << c6 << ";" << endl
       << "const double qOO = " << setw(20) << qOO << ";" << endl
       << "const double qOH = " << setw(20) << qOH << ";" << endl
       << "const double qHH = " << setw(20) << qHH << ";" << endl
       << "const double solvent_permittivity = " << setw(20) << solvent_permittivity << ";" << endl
       << "const double solvent_kappa = " << setw(20) << solvent_kappa << ";" << endl
       << "const double four_pi_eps_i = " << setw(20) << math::four_pi_eps_i << ";" << endl
       << "#endif" << endl; 
  
  /////////////////////////////////////////////////
  simulation::Simulation sim;
  sim.param().nonbonded.rf_epsilon = solvent_permittivity;
  sim.param().nonbonded.rf_cutoff = longrange_cutoff;
  sim.param().nonbonded.rf_kappa = solvent_kappa;
  sim.param().force.interaction_function = simulation::lj_crf_func;
  
  interaction::Nonbonded_Term term;
  term.init(sim);  
  /////////////////////////////////////////////////
  // let's tabulate the short range first
  /////////////////////////////////////////////////
  
  cout << "#ifdef SHORTRANGE" << endl;
  
  double data_start = 0.0, data_end = shortrange_cutoff + solvent_diameter;
  
  double data_range = (data_end - data_start) * (data_end - data_start);
  // take square because r**2 is tabulated
  double step = data_range / shortrange_table_size; 
  
  math::Vec r(0.0), r_next(0.0, 0.0, 0.0);
  double x = 0.0;
  
  cout << "static const unsigned int table_size = " << shortrange_table_size << ";" << endl
       << "static const double data_range = " << data_range << ";" << endl
       << endl;
  cout << "static const lj_crf_values OO_table[table_size] = {" << endl;
  for(unsigned int i = 0; i < shortrange_table_size; ++i, x += step) {
    r(0) = sqrt(x);
    r_next(0) = sqrt(x + step);
    double e_lj = 0.0, e_crf = 0.0, f = 0.0;
    double e_lj_next = 0.0, e_crf_next = 0.0, f_next = 0.0;
    if (x) // avoid x == 0.0
      term.lj_crf_interaction(r, c6, c12, qOO, f, e_lj, e_crf);
    
    term.lj_crf_interaction(r_next, c6, c12, qOO, f_next, e_lj_next, e_crf_next);
    
    cout << "   { " << setw(20) << e_lj << "," << setw(20) << e_lj_next - e_lj << ","
         << setw(20) << e_crf << "," << setw(20) << e_crf_next - e_crf << ","
         << setw(20) << f << "," << setw(20) << f_next - f << "}";
    if (i != shortrange_table_size - 1) {
      cout << ", ";
    }
    cout << endl;
  }
  cout << "};" << endl;
  cout << "static const crf_values OH_table[table_size] = {" << endl;
  x = 0.0;
  for(unsigned int i = 0; i < shortrange_table_size; ++i, x += step) {
    r(0) = sqrt(x);
    r_next(0) = sqrt(x + step);
    double e_lj = 0.0, e_crf = 0.0, f = 0.0;
    double e_lj_next = 0.0, e_crf_next = 0.0, f_next = 0.0;
    if (x)
      term.lj_crf_interaction(r, 0.0, 0.0, qOH, f, e_lj, e_crf);
    term.lj_crf_interaction(r_next, 0.0, 0.0, qOH, f_next, e_lj_next, e_crf_next);
    cout << "   { " << setw(20) << e_crf << "," << setw(20) << e_crf_next - e_crf << 
            "," << setw(20) << f << "," << setw(20) << f_next - f << "}";
    if (i != shortrange_table_size - 1) {
      cout << ", ";
    }
    cout << endl;
  }
  cout << "};" << endl;
  cout << "static const crf_values HH_table[table_size] = {" << endl;
  x = 0.0;
  for(unsigned int i = 0; i < shortrange_table_size; ++i, x += step) {
    r(0) = sqrt(x);
    r_next(0) = sqrt(x + step);
    double e_lj = 0.0, e_crf = 0.0, f = 0.0;
    double e_lj_next = 0.0, e_crf_next = 0.0, f_next = 0.0;
    if (x)
      term.lj_crf_interaction(r, 0.0, 0.0, qHH, f, e_lj, e_crf);
    term.lj_crf_interaction(r_next, 0.0, 0.0, qHH, f_next, e_lj_next, e_crf_next);
    cout << "   { " << setw(20) << e_crf << "," << setw(20) << e_crf_next - e_crf << 
            "," << setw(20) << f << "," << setw(20) << f_next - f << "}";
    if (i != shortrange_table_size - 1) {
      cout << ", ";
    }
    cout << endl;
  }
  cout << "};" << endl;
  cout << "#endif" << endl;
  cout << "#ifdef LONGRANGE" << endl;
  /////////////////////////////////////////////////
  // let's tabulate the longrange now
  /////////////////////////////////////////////////
  data_start = shortrange_cutoff - solvent_diameter;
  data_end = longrange_cutoff + solvent_diameter;
  
  data_range = data_end * data_end - data_start * data_start;
  // take square because r**2 is tabulated
  step = data_range / longrange_table_size; 
  
  double r2_start = data_start * data_start;
  
  x = r2_start;
  
  cout << "static const unsigned int table_size = " << longrange_table_size << ";" << endl
       << "static const double data_start = " << r2_start << ";" << endl
       << "static const double data_range = " << data_range << ";" << endl
       << endl;
  cout << "static const lj_crf_values OO_table[table_size] = {" << endl;
  for(unsigned int i = 0; i < longrange_table_size; ++i, x += step) {
    r(0) = sqrt(x);
    r_next(0) = sqrt(x + step);
    double e_lj = 0.0, e_crf = 0.0, f = 0.0;
    double e_lj_next = 0.0, e_crf_next = 0.0, f_next = 0.0;
    
    term.lj_crf_interaction(r, c6, c12, qOO, f, e_lj, e_crf);
    term.lj_crf_interaction(r_next, c6, c12, qOO, f_next, e_lj_next, e_crf_next);
    
    cout << "   { " << setw(20) << e_lj << "," << setw(20) << e_lj_next - e_lj << ","
         << setw(20) << e_crf << "," << setw(20) << e_crf_next - e_crf << ","
         << setw(20) << f << "," << setw(20) << f_next - f << "}";
    if (i != longrange_table_size - 1) {
      cout << ", ";
    }
    cout << endl;
  }
  cout << "};" << endl;
  cout << "static const crf_values OH_table[table_size] = {" << endl;
  x = r2_start;
  for(unsigned int i = 0; i < longrange_table_size; ++i, x += step) {
    r(0) = sqrt(x);
    r_next(0) = sqrt(x + step);
    double e_lj = 0.0, e_crf = 0.0, f = 0.0;
    double e_lj_next = 0.0, e_crf_next = 0.0, f_next = 0.0;

    term.lj_crf_interaction(r, 0.0, 0.0, qOH, f, e_lj, e_crf);
    term.lj_crf_interaction(r_next, 0.0, 0.0, qOH, f_next, e_lj_next, e_crf_next);
    
    cout << "   { " << setw(20) << e_crf << "," << setw(20) << e_crf_next - e_crf << 
            "," << setw(20) << f << "," << setw(20) << f_next - f << "}";
    if (i != longrange_table_size - 1) {
      cout << ", ";
    }
    cout << endl;
  }
  cout << "};" << endl;
  cout << "static const crf_values HH_table[table_size] = {" << endl;
  x = r2_start;
  for(unsigned int i = 0; i < longrange_table_size; ++i, x += step) {
    r(0) = sqrt(x);
    r_next(0) = sqrt(x + step);
    double e_lj = 0.0, e_crf = 0.0, f = 0.0;
    double e_lj_next = 0.0, e_crf_next = 0.0, f_next = 0.0;

    term.lj_crf_interaction(r, 0.0, 0.0, qHH, f, e_lj, e_crf);
    term.lj_crf_interaction(r_next, 0.0, 0.0, qHH, f_next, e_lj_next, e_crf_next);
    
    cout << "   { " << setw(20) << e_crf << "," << setw(20) << e_crf_next - e_crf << 
            "," << setw(20) << f << "," << setw(20) << f_next - f << "}";
    if (i != longrange_table_size - 1) {
      cout << ", ";
    }
    cout << endl;
  }
  cout << "};" << endl;
  cout << "#endif" << endl;
}

