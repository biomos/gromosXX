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
 * @file fit_rdc.cc
 * fit (an) alignment tensor(s) to structures and experimental RDCs, and back calculate RDCs
 */

/**
 * @page contrib Contrib Program Documentation
 *
 * @anchor fit_rdc
 * @section fit_rdc fit a(n ensemble of) structure(s) to RDC(s)
 * @author jra, lnw
 * @date 30.08.2013, improved mid-2015
 *
 * PROGRAM DESCRIPTION
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td> &lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; [&lt;gathermethod&gt;] (don't use glist, gltime or gref as there is no list); </td></tr>
 * <tr><td> \@fitspec</td><td>&lt;file containing data to fit to (the rdc data from the RDCRESSPEC block)&gt; </td></tr>
 * <tr><td> [\@bcspec</td><td>&lt;file containing data to backcalculate (if different to \@fitspec)&gt]; </td></tr>
 * <tr><td> [\@framespec</td><td>&lt;frames to consider for the fit: ALL (default), EVERY or SPEC&gt]; </td></tr>
 * <tr><td> [\@frames</td><td>&lt;frames to consider for the fit&gt; (if framespec EVERY or SPEC)] </td></tr>
 * <tr><td> \@traj</td><td>&lt;trajectory files&gt; </td></tr>
 * <tr><td> [\@fit</td><td>&lt;fitting method, either LLS (default) or SVD&gt;] </td></tr>
 * <tr><td> [\@verbose</td><td>&lt;print more (OFF(default)/NO/ON/YES)&gt;] </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
 * fit_rdc
 *    @topo       ex.top
 *    @pbc        r
 *    @fitspec    ex.fit
 *    @traj       ex.trj
 @endverbatim
 *
 * <hr>
 */

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/bound/Boundary.h"
#include "../src/gcore/System.h"
#include "../src/gio/Ginstream.h"
#include "../src/gio/InG96.h"
#include "../src/gio/InTopology.h"
#include "../src/utils/Neighbours.h"
#include "../src/utils/PropertyContainer.h"
#include "../src/utils/RdcFuncs.h"
#include "../src/utils/debug.h"
#include "../src/utils/groTime.h"
#include "../src/gromos/Exception.h"

#if (__cplusplus > 199711L) // we have c++11 or newer
#include <chrono>
#include <ctime>
#endif

using namespace std;

using namespace args;
using namespace utils;


template <typename T, typename U>
ostream &operator<<(ostream &s, const pair<T, U> &v) {
  s << "<" << v.first << "," << v.second << ">";
  return s;
}
#define container_output(container)                                              \
  template <typename T> ostream &operator<<(ostream &s, const container<T> &v) { \
    s << "{";                                                                    \
    for (typename container<T>::const_iterator x(v.begin()); x != v.end();) {    \
      s << *x;                                                                   \
      if (++x != v.end())                                                        \
        s << ",";                                                                \
    }                                                                            \
    s << "}";                                                                    \
    return s;                                                                    \
  }
container_output(vector);

int main(int argc, char **argv) {
  // print the time and command for future reference
#if (__cplusplus > 199711L) // we have c++11 or newer
  std::chrono::time_point<std::chrono::system_clock> now = std::chrono::system_clock::now();
  std::time_t now_time = std::chrono::system_clock::to_time_t(now);
  std::cout << "time: " << std::ctime(&now_time) ;
#endif
  cout <<  "command: ";
  for (int i=0; i<argc; i++) cout << argv[i] << " ";
  cout << endl;

  Argument_List knowns;
  knowns << "topo"
         << "pbc"
         << "fitspec"
         << "bcspec"
         << "framespec"
         << "frames"
         << "traj"
         << "fit"
         << "verbose";

  string usage = "# " + string(argv[0]) + "\n";
  usage += "\t@topo        <molecular topology file>\n";
  usage += "\t@pbc         <boundary type> [<gathermethod>] (don't use glist, gltime or gref as there is no list)\n";
  usage += "\t@fitspec     <file containing data to fit to (the rdc data from the RDCRESSPEC block)>\n";
  usage += "\t[@bcspec     <file containing data to backcalculate> (if different to @fitspec)]\n";
  usage += "\t[@framespec  <frames to consider for the fit: ALL (default), EVERY or SPEC>]\n";
  usage += "\t[@frames     <frames to consider for the fit> (if framespec EVERY or SPEC)]\n";
  usage += "\t@traj        <trajectory files>\n";
  usage += "\t[@fit        <fitting method, either LLS (default) or SVD>]\n";
  usage += "\t[@verbose    <print more (OFF(default)/NO/ON/YES)>]\n";

  const int n_ah = 5;
  const double ps_inv2s_inv = 1e12;

  try {
    Arguments args(argc, argv, knowns, usage);

    // -- read 'topo', the topology
    gio::InTopology in_topo(args["topo"]);
    gcore::System sys(in_topo.system());

    // -- parse 'framespec',  and 'frames' ie, which frames of the trajectory do we want?
    string framespec = "ALL";
    vector<int> frames;
    if (args.count("framespec") > 0) {
      framespec = args["framespec"];
      if (framespec != "ALL" && framespec != "EVERY" && framespec != "SPEC") {
        throw gromos::Exception("fit_rdc", "framespec format " + framespec + " unknown.");
      }
      if (framespec == "EVERY" || framespec == "SPEC") {
        for (Arguments::const_iterator it = args.lower_bound("frames"), to = args.upper_bound("frames"); it != to; ++it) {
          const int foo = atoi(it->second.c_str());
          frames.push_back(foo);
        }
        if (frames.size() == 0) {
          throw gromos::Exception("fit_rdc", "if you give EVERY or SPEC you have to use @frames as well");
        }
        if (frames.size() != 1 && framespec == "EVERY") {
          throw gromos::Exception("fit_rdc", "if you give EVERY you have to give exactly one number with @frames");
        }
      }
      if (framespec == "SPEC") {
        std::sort(frames.begin(), frames.end());
      } // sort frames, so processing the last entry means we're done
    }

    // only used in case of 'ALL' or 'EVERY'
    // store calculated and exp D values, such that averaged Q can be calculated
    // <frame<#fit-groups<#bc-groups<# rdcs in bc group<D,Dexp>>>>
    vector<vector<vector<vector<pair<double, double> > > > > D_calc_exp;

    DEBUG(5, "timepoints: " << frames)

    // -- parse 'pbc' that is the boxtype/ the gathering method
    bound::Boundary::MemPtr gathmethod = args::GatherParser::parse(sys, sys, args); // the second 'sys' is supposed to be refsys, but there is no refsys
    bound::Boundary *boundary = BoundaryParser::boundary(sys, args);

    // -- parse 'fitspec', the RDCRESSPEC block
    gio::Ginstream fitspec_stream(args["fitspec"]);
    vector<string> buffer, buffer_fitrdc, buffer_groups; // will be reused for bc
    // first search for RDCRESSPEC block
    bool found_rdc = false, found_groups = false;
    while (!fitspec_stream.stream().eof() && (!found_rdc || !found_groups)) {
      fitspec_stream.getblock(buffer);
      if (buffer[0] == "RDCRESSPEC") {
        found_rdc = true;
        buffer_fitrdc = buffer;
      }
      if (buffer[0] == "RDCGROUPS") {
        found_groups = true;
        buffer_groups = buffer;
      }
    }
    fitspec_stream.close();
    if (!found_rdc) {
      throw gromos::Exception("fit_rdc", "RDC fit file does not contain an RDCRESSPEC block!");
    }
    if (!found_groups) {
      throw gromos::Exception("fit_rdc", "RDC fit file does not contain an RDCGROUPS block!");
    }
    if (buffer_fitrdc[buffer_fitrdc.size() - 1].find("END") != 0) {
      throw gromos::Exception("fit_rdc", "RDC fit file is corrupted. No END in " + buffer_fitrdc[0] + " block. Got\n" + buffer_fitrdc[buffer_fitrdc.size() - 1]);
    }
    if (buffer_groups[buffer_groups.size() - 1].find("END") != 0) {
      throw gromos::Exception("fit_rdc", "RDC fit file is corrupted. No END in " + buffer_groups[0] + " block. Got\n" + buffer_groups[buffer_groups.size() - 1]);
    }

    // read in the RDC data for fitting
    DEBUG(7, "rdc buffer: " << buffer_fitrdc)
    rdcdata_t fit_dat_flat = read_rdc(buffer_fitrdc, sys, /* bool fit=*/true);
    DEBUG(5, "fit_dat_flat: " << fit_dat_flat)

    // read alignment blocks
    DEBUG(7, "groups_buffer: " << buffer_groups)
    vector<vector<unsigned int> > fit_groups = read_groups(buffer_groups, fit_dat_flat.size());
    DEBUG(5, "fit_groups: " << fit_groups)

    // split rdcs according to groups
    vector<rdcdata_t> fit_dat;
    for (unsigned int i = 0; i < fit_groups.size(); i++) {
      fit_dat.push_back(rdcdata_t());
      for (unsigned int j = 0; j < fit_groups[i].size(); j++) {
        fit_dat[i].push_back(fit_dat_flat[fit_groups[i][j] - 1]);
      }
    }
    DEBUG(5, "grouped rdcs: " << fit_dat)

    vector<unsigned int> n_rdc_fit(fit_groups.size(), 0);
    for (unsigned int i = 0; i < fit_groups.size(); i++) {
      n_rdc_fit[i] = fit_groups[i].size();
    }
    DEBUG(5, "n rdc fit: " << n_rdc_fit)

    // -- parse 'bcspec', check if we have different data to back-calculate to
    bool backcalc = false;
    vector<unsigned int> n_rdc_bc;
    rdcdata_t bc_dat_flat;
    vector<rdcdata_t> bc_dat;
    if (args.count("bcspec") > 0) { // we back-calc to a different data set
      backcalc = true;

      // read in the data to back-calculate from file
      gio::Ginstream bcspec_stream(args["bcspec"]);
      vector<string> buffer_bcrdc;
      // first search for RDCRESSPEC block
      bool found_rdc = false, found_groups = false;
      while (!bcspec_stream.stream().eof() && (!found_rdc || !found_groups)) {
        bcspec_stream.getblock(buffer);
        if (buffer[0] == "RDCRESSPEC") {
          found_rdc = true;
          buffer_bcrdc = buffer;
        }
        if (buffer[0] == "RDCGROUPS") {
          found_groups = true;
          buffer_groups = buffer;
        }
      }
      bcspec_stream.close();
      if (!found_rdc) {
        throw gromos::Exception("fit_rdc", "RDC BC file does not contain an RDCRESSPEC block!");
      }
      if (!found_groups) {
        throw gromos::Exception("fit_rdc", "RDC BC file does not contain an RDCGROUPS block!");
      }
      if (buffer_bcrdc[buffer_bcrdc.size() - 1].find("END") != 0) {
        throw gromos::Exception("fit_rdc", "RDC BC file is corrupted. No END in " + buffer_bcrdc[0] + " block. Got\n" + buffer_bcrdc[buffer_bcrdc.size() - 1]);
      }
      if (buffer_groups[buffer_groups.size() - 1].find("END") != 0) {
        throw gromos::Exception("fit_rdc", "RDC BC file is corrupted. No END in " + buffer_groups[0] + " block. Got\n" + buffer_groups[buffer_groups.size() - 1]);
      }

      DEBUG(7, "rdc buffer: " << buffer_bcrdc)
      bc_dat_flat = read_rdc(buffer_bcrdc, sys, /*bool fit=*/false);
      DEBUG(5, "bc_dat_flat: " << bc_dat_flat)

      // read alignment blocks in bc
      DEBUG(7, "group buffer: " << buffer_groups)
      vector<vector<unsigned int> > bc_groups = read_groups(buffer_groups, bc_dat_flat.size());
      DEBUG(5, "groups: " << bc_groups)

      // split rdcs according to groups
      for (unsigned int i = 0; i < bc_groups.size(); i++) {
        bc_dat.push_back(rdcdata_t());
        for (unsigned int j = 0; j < bc_groups[i].size(); j++) {
          bc_dat[i].push_back(bc_dat_flat[bc_groups[i][j] - 1]);
        }
      }
      DEBUG(5, "grouped rdcs: " << fit_dat)

      n_rdc_bc.resize(bc_groups.size(), 0);
      for (unsigned int i = 0; i < bc_groups.size(); i++) {
        n_rdc_bc[i] = bc_groups[i].size();
      }
      DEBUG(5, "n rdc bc: " << n_rdc_bc)
    } else { // we back-calc to the same data set we fitted the tensor to previously
      bc_dat = fit_dat;
      n_rdc_bc = n_rdc_fit;
    }

    // -- parse 'fit'
    bool do_lls = true, do_svd = false;
    if (args.count("fit") == 1) {
      if (args["fit"] == "LLS" || args["fit"] == "lls") {
      } else if (args["fit"] == "SVD" || args["fit"] == "svd") {
        do_lls = false;
        do_svd = true;
      } else {
        throw gromos::Exception("fit_rdc", "Invalid fitting method (valid values are 'LLS' and 'SVD')");
      }
    }

    // -- parse 'verbose'
    bool verbose = false;
    if (args.count("verbose") == 1) {
      if (args["verbose"] == "NO" || args["verbose"] == "no" || args["verbose"] == "OFF" || args["verbose"] == "off") {
      } else if (args["verbose"] == "YES" || args["verbose"] == "yes" || args["verbose"] == "ON" || args["verbose"] == "on") {
        verbose = true;
      } else {
        throw gromos::Exception("fit_rdc", "Invalid verbosity choice (valid are OFF/NO/ON/YES)')");
      }
    }

    // ---------------------------------------------------------------------------------------

    if(args.count("traj") <= 0) 
      throw gromos::Exception("fit_rdc", "a trajectory is required");

    // start at -1 to get times right. this is for the framespec.
    int current_frame = -1;
    int current_frame_of_traj = -1;
    int frames_done = 0; // required only for 'spec'
    int to_skip = -1;    // unused skips/stride steps of the last trajectory

    // loop over all trajectories
    for (Arguments::const_iterator it = args.lower_bound("traj"), to = args.upper_bound("traj"); it != to; ++it) {

      // open file
      gio::InG96 cartesian_traj((it->second).c_str());
      cartesian_traj.select("SOLUTE");

      current_frame_of_traj = -1;

      // loop over single trajectory
      while (!cartesian_traj.eof()) {
        if (framespec == "ALL") {
          current_frame++, current_frame_of_traj++;
        }
        if (framespec == "EVERY") {
          if (to_skip != -1) { // this is not the first traj and there is half a stride left from the last one
            cartesian_traj.skip(to_skip);
            current_frame += to_skip + 1;
            current_frame_of_traj += to_skip + 1;
            to_skip = -1;
          } else {
            if (current_frame != -1 && cartesian_traj.stride() == 1) { // set stride if not set yet, unless this is the first frame
              cartesian_traj.stride(frames[0]);
            }
            const int increment = current_frame_of_traj == -1 ? 1 : frames[0];
            current_frame += increment, current_frame_of_traj += increment;
          }
        }
        if (framespec == "SPEC") {
          if (to_skip != -1) {
            cartesian_traj.skip(to_skip);
            current_frame += to_skip + 1;
            current_frame_of_traj += to_skip + 1;
            to_skip = -1;
            frames_done++;
          } else { // skip to the next chosen frame (either from the last frame or the first frame)
            cartesian_traj.skip(frames[frames_done] - (current_frame == -1 ? -1 : (frames[frames_done - 1])) - 1);
            current_frame += frames[frames_done] - (current_frame == -1 ? -1 : frames[frames_done - 1]);
            frames_done++;
          }
        }

        Time time;
        cartesian_traj >> sys >> time;
        if (!sys.hasPos) {
          throw gromos::Exception("fit_rdc", "Unable to read POSITION(RED) block from trajectory file.");
        }
        if (cartesian_traj.stride_eof()) { // when skipping or striding beyond the end
          to_skip = cartesian_traj.skip(); // to do in the next traj
          current_frame -= to_skip + 1;    // and we counted to far previously
          frames_done--;                   // not done yet
          break;
        }

        // gather
        (boundary->*gathmethod)();

        if (framespec == "ALL" || framespec == "EVERY") {
          D_calc_exp.push_back(vector<vector<vector<pair<double, double> > > >()); // start next frame
        }

        for (unsigned int fit_group = 0; fit_group < fit_dat.size(); fit_group++) {
          if (framespec == "ALL" || framespec == "EVERY") {
            D_calc_exp.back().push_back(vector<vector<pair<double, double> > >()); // start next fit-group
          }

          for (unsigned int bc_group = 0; bc_group < bc_dat.size(); bc_group++) {
            if (framespec == "ALL" || framespec == "EVERY") {
              D_calc_exp.back().back().push_back(vector<pair<double, double> >()); // start next bc-group
            }

            cout << endl
                 << "# frame " << current_frame << " at" << time << "ps, block " << fit_group << endl;

            // but if !bc we only want fit_group=bc_group
            if (!backcalc)
              bc_group = fit_group;

            double coef_mat[n_rdc_fit[fit_group] * n_ah];
            for (unsigned int i = 0; i < n_rdc_fit[fit_group] * n_ah; i++) {
              coef_mat[i] = 0;
            }
            double coef_mat_bc_j[n_rdc_bc[bc_group] * n_ah];
            for (unsigned int i = 0; i < n_rdc_bc[bc_group] * n_ah; i++) {
              coef_mat_bc_j[i] = 0;
            }
            double coef_mat_bc_k[n_rdc_bc[bc_group] * n_ah];
            for (unsigned int i = 0; i < n_rdc_bc[bc_group] * n_ah; i++) {
              coef_mat_bc_k[i] = 0;
            }

            vector<double> exp_fit_norm(n_rdc_fit[fit_group]); // of the data set we fit to: exp/dmax
            vector<double> exp_bc_norm(n_rdc_bc[bc_group]);    // of the data set we back calc to: exp/dmax
            for (unsigned int i = 0; i < n_rdc_fit[fit_group]; i++) {
              exp_fit_norm[i] = fit_dat[fit_group][i].exp / fit_dat[fit_group][i].dmax;
            }
            for (unsigned int i = 0; i < n_rdc_bc[bc_group]; i++) {
              exp_bc_norm[i] = bc_dat[bc_group][i].exp / bc_dat[bc_group][i].dmax; // yes, I know, but it's cheap and vectorised
            }

            // calculate the coefficients for the RDCs to fit to and the
            // inter-nuclear distances if fit_rij == ALL and find HA coords if
            // there are CA-HA RDCs
            DEBUG(10, "fit_group: " << fit_group)
            DEBUG(10, "fit_dat[fit_group].size(): " << fit_dat[fit_group].size())
            DEBUG(10, "fit_dat[fit_group]: " << fit_dat[fit_group])
            calc_coef_fit(sys, fit_dat[fit_group], coef_mat);

            // and for the RDCs to back-calculate (if different)
            if (backcalc) {
              calc_coef_bc(sys, bc_dat[bc_group], coef_mat_bc_j, coef_mat_bc_k);
            } else {
              std::copy(coef_mat, coef_mat + n_rdc_fit[fit_group] * n_ah, coef_mat_bc_j);
              // and coef_mat_bc_k may stay all zeros
            } // if not, the coefficients were already calculated

            // -- alignment tensor --
            vector<double> tensor5;
            if (do_lls) {
              tensor5 = lls_fit(sys, fit_dat[fit_group], coef_mat);
            } else if (do_svd) {
              // coef_mat gets modified but is not used again
              tensor5 = svd_fit(sys, fit_dat[fit_group], coef_mat);
            } else {
              assert(false);
            }
            DEBUG(7, "tensor: " << tensor5)

            vector<double> eigenvalues;
            vector<vector<double> > eigenvectors;
            diagonalize_tensor(tensor5, eigenvalues, eigenvectors);
            DEBUG(7, "eigenvalues: " << eigenvalues)
            DEBUG(7, "eigenvectors: " << eigenvectors)

            cout.precision(10);
            cout.setf(ios::right, ios::adjustfield);
            if (verbose) {
              cout << "# Alignment tensor" << endl;
              cout << "# Sxx = " << tensor5[0] << endl;
              cout << "# Syy = " << tensor5[1] << endl;
              cout << "# Szz = " << -tensor5[0] - tensor5[1] << endl;
              cout << "# Sxy = " << tensor5[2] << endl;
              cout << "# Sxz = " << tensor5[3] << endl;
              cout << "# Syz = " << tensor5[4] << endl;
            }

            cout << "# Eigenvalues" << endl;
            // absolute values of EV in ascending order
            cout << "# Sxx_d = " << eigenvalues[0] << endl;
            cout << "# Syy_d = " << eigenvalues[1] << endl;
            cout << "# Szz_d = " << eigenvalues[2] << endl;

            if (verbose) {
              cout << "# Eigenvectors   x_coor   y_coor   z_coor" << endl;
              cout << "# x-axis " << eigenvectors[0][0] << " " << eigenvectors[0][1] << " " << eigenvectors[0][2] << endl;
              cout << "# y-axis " << eigenvectors[1][0] << " " << eigenvectors[1][1] << " " << eigenvectors[1][2] << endl;
              cout << "# z-axis " << eigenvectors[2][0] << " " << eigenvectors[2][1] << " " << eigenvectors[2][2] << endl;
            }

            // -- alignment parameters --
            double Da = 0.5 * eigenvalues[2];
            double Dr = (eigenvalues[1] - eigenvalues[0]) / 3.0;
            double Aa = 2.0 * Da; // axial component
            double Ar = 2.0 * Dr; // rhombic component
            double R = Ar / Aa;   // rhombicity
            cout << "# Alignment parameters" << endl;
            cout << "# Da = " << Da << endl;
            cout << "# Dr = " << Dr << endl;
            cout << "# Aa (axial component): " << Aa << endl;
            cout << "# Ar (rhombic component): " << Ar << endl;
            cout << "# R (rhombicity): " << R << endl;

            // -- back-calculate the bc-RDCs using S
            vector<double> calc_bc_norm_j(n_rdc_bc[bc_group], 0);
            vector<double> calc_bc_norm_k(n_rdc_bc[bc_group], 0);
            for (unsigned int k = 0; k < n_rdc_bc[bc_group]; k++) {
              for (int h = 0; h < n_ah; h++) {
                calc_bc_norm_j[k] += coef_mat_bc_j[k * n_ah + h] * tensor5[h];
                calc_bc_norm_k[k] += coef_mat_bc_k[k * n_ah + h] * tensor5[h];
              }
            }

            // -- sum the back-calculated RDCs (this only matters for side-chain RDCs; for the rest, nothing happens)
            vector<double> calc_bc_norm(n_rdc_bc[bc_group]);
            for (unsigned int k = 0; k < n_rdc_bc[bc_group]; k++) {
              calc_bc_norm[k] = calc_bc_norm_j[k] + calc_bc_norm_k[k];
            }

            // -- un-normalise the RDCs
            vector<double> exp_bc(n_rdc_bc[bc_group]), calc_bc(n_rdc_bc[bc_group]); // experimental and back calculated rdc, both from the bc data set
            for (unsigned int k = 0; k < n_rdc_bc[bc_group]; k++) {
              exp_bc[k] = bc_dat[bc_group][k].exp;
              calc_bc[k] = calc_bc_norm[k] * bc_dat[bc_group][k].dmax;
              DEBUG(7, "calc_bc_norm: " << scientific << calc_bc_norm[k])
              DEBUG(7, "dmax: " << scientific << bc_dat[bc_group][k].dmax)
              DEBUG(7, "exp_bc: " << scientific << exp_bc[k])
              DEBUG(7, "calc_bc: " << scientific << calc_bc[k])
            }

            // -- calculate and print Q-values, R-values and rmsd
            // for normalised RDCs (according to Bax-2001)
            cout << "# Q-values for back-calculated data" << endl;
            const double Q_norm = calc_Q(calc_bc_norm, exp_bc_norm);
            cout << "# Q (Normalized) (PALES) = " << Q_norm << endl;
            // for un-normalised RDCs
            const double Q_raw = calc_Q(calc_bc, exp_bc);
            cout << "# Q (Raw) (PALES)  = " << Q_raw << endl;

            // R-values
            cout << "# R-values for back-calculated data, crystallographic" << endl;
            const double R_norm = calc_R(calc_bc_norm, exp_bc_norm);
            cout << "# R (Normalized) = " << R_norm << endl;
            const double R_raw = calc_R(calc_bc, exp_bc);
            cout << "# R (Raw) = " << R_raw << endl;

            // rmsd
            cout << "# RMSD for back-calculated data" << endl;
            const double RMSD_raw = calc_RMSD(calc_bc, exp_bc) * ps_inv2s_inv;
            cout << "# RMSD (Raw) = " << RMSD_raw << " Hz" << endl;

            // -- print the experimental and back-calculated RDCs --
            cout.precision(5);
            // title
            cout << setw(18) << "#   residue/atom 1"
                 << setw(18) << "    residue/atom 2"
                 << setw(18) << "    residue/atom 3"
                 << setw(12) << "D_exp"
                 << setw(12) << "D_calc"
                 << endl;

            // now print actual (back-calculated) RDCs
            for (unsigned int i = 0; i < n_rdc_bc[bc_group]; i++) {
              double calc_rdc = calc_bc[i] * ps_inv2s_inv;
              // k may be zero if not a side-chain NH RDC
              int katom, kres;
              string kname;
              if (bc_dat[bc_group][i].type == 5 || bc_dat[bc_group][i].type == 6) {
                kres = sys.mol(bc_dat[bc_group][i].mol).topology().resNum(bc_dat[bc_group][i].k) + 1;
                katom = bc_dat[bc_group][i].k + 1;
                kname = sys.mol(bc_dat[bc_group][i].mol).topology().atom(bc_dat[bc_group][i].k).name();
              } else {
                katom = 0;
                kname = "--";
                kres = 0;
                // abs value for HH RDCs
                if (bc_dat[bc_group][i].type == 7 || bc_dat[bc_group][i].type == 8) {
                  calc_rdc = std::abs(calc_rdc);
                }
              }

              cout << setw(6) << sys.mol(bc_dat[bc_group][i].mol).topology().resNum(bc_dat[bc_group][i].i) + 1
                   << setw(6) << sys.mol(bc_dat[bc_group][i].mol).topology().atom(bc_dat[bc_group][i].i).name()
                   << setw(6) << bc_dat[bc_group][i].i + 1
                   << setw(6) << sys.mol(bc_dat[bc_group][i].mol).topology().resNum(bc_dat[bc_group][i].j) + 1
                   << setw(6) << sys.mol(bc_dat[bc_group][i].mol).topology().atom(bc_dat[bc_group][i].j).name()
                   << setw(6) << bc_dat[bc_group][i].j + 1
                   << setw(6) << kres
                   << setw(6) << kname
                   << setw(6) << katom
                   << setw(12) << setprecision(4) << exp_bc[i] * ps_inv2s_inv
                   << setw(12) << setprecision(4) << calc_rdc
                   << endl;

              if (framespec == "ALL" || framespec == "EVERY") {
                D_calc_exp.back().back().back().push_back(make_pair(exp_bc[i] * ps_inv2s_inv, calc_rdc));
              }
            } // end loop over RDCs to print

            // but if !bc we only want fit_group=bc_group
            if (!backcalc)
              break;
          } // bc groups
        }   // fit groups

        if (framespec == "SPEC" && current_frame == frames.back())
          break;
      } // ends loop over this trj

      cartesian_traj.close();
    } // end loop over all trj

    // calculate and print averaged Q values for every fit group / bc group combination
    if (framespec == "ALL" || framespec == "EVERY") {
      DEBUG(5, "stored D-calc / D-exp: " << D_calc_exp)
      DEBUG(5, "dimensions of D_calc_exp: " << D_calc_exp.size() << ", " << D_calc_exp.front().size() << ", " << D_calc_exp.front().front().size())
      for (unsigned i = 0; i < D_calc_exp.front().front().size(); i++) {
        DEBUG(5, "dimensions of D_calc_exp (last level): " << D_calc_exp.front().front()[i].size())
      }
      const unsigned n_frames = D_calc_exp.size();
      vector<vector<vector<double> > > D_average;                // fit<bc<D>>
      for (unsigned i = 0; i < D_calc_exp.front().size(); i++) { // all fit groups in the first frame
        D_average.push_back(vector<vector<double> >());
        for (unsigned j = 0; j < D_calc_exp.front()[i].size(); j++) {
          D_average[i].push_back(vector<double>(D_calc_exp.front()[i][j].size(), 0));
        }
      }
      for (unsigned h = 0; h < n_frames; h++) {
        for (unsigned i = 0; i < D_average.size(); i++) {
          for (unsigned j = 0; j < D_average[i].size(); j++) {
            for (unsigned k = 0; k < D_average[i][j].size(); k++) {
              DEBUG(8, "ijk: " << i << ", " << j << ", " << k << " Dav,increment: " << D_average[i][j][k] << ", " << D_calc_exp[h][i][j][k].second / n_frames)
              D_average[i][j][k] += D_calc_exp[h][i][j][k].second / n_frames;
            }
          }
        }
      }
      DEBUG(5, "D_average: " << D_average)

      cout << endl
           << "Averaged Q values (all frames with equal weight), according to all combinations of fit groups and bc groups" << endl;
      vector<vector<double> > Q_average(D_average.size(), vector<double>(D_average[0].size(), 0));
      for (unsigned i = 0; i < D_average.size(); i++) {
        for (unsigned j = 0; j < D_average[i].size(); j++) {
          double nom = 0.0, denom = 0.0;
          for (unsigned k = 0; k < D_average[i][j].size(); k++) {
            const double av = D_average[i][j][k];
            const double exp = D_calc_exp.front()[i][j][k].first; // read experimental values from first frame
            nom += pow(av - exp, 2);
            denom += exp * exp;
            DEBUG(7, "av, exp: " << av << ", " << exp)
          }
          DEBUG(7, "n/d: " << nom << ", " << denom)
          Q_average[i][j] = sqrt(nom / denom);
          cout << "Q(" << i << "," << (backcalc ? j : i) << "): " << Q_average[i][j] << endl;
        }
      }
      DEBUG(5, "Q_average: " << Q_average)
    }

    delete boundary;
  } // try
  catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}
