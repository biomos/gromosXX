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
 * @file cg_potential.cc
 * calculates the potential of molecules on a grid
 */

/**
 * @page contrib Contrib Program Documentation
 *
 * @anchor cg_potential
 * @section cg_potential calculates coarse grain potentials from a fine grain trajectory
 * @author @ref ns
 * @date 10.10.2008
 *
 * In a first step program cg_potential groups solvent molecules in clusters of
 * a given size. Second, it calculates the energies between the
 * clusters and stores them as a function of the cluster's COG distances.
 *
 * The cluster size and shape of the clusters can be controlled using the 
 * \@cluster_size and \@max_cluster_rmsd arguments. The energy is calculated on
 * a grid with \@grid_size points from 0.0 nm up to a given cutoff (\@cutoff 
 * argument). Due to the repulsive nature of the interaction potentials
 * at short distances sampling close to 0.0 nm will be bad. The reaction field
 * parameters (RF cutoff, permittivity and inverse Debeye screeing length) can be
 * controlled by \@reaction_field.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;periodic boundary conditions&gt; </td></tr>
 * <tr><td> \@cutoff</td><td>&lt;cutoff for intercluster distance&gt;</td></tr>
 * <tr><td> \@traj</td><td>&lt;trajectory files&gt;</td></tr>
 * <tr><td>[\@cluster_size</td><td>&lt;number of solvent molecules that build a cluster&gt;]</td></tr>
 * <tr><td>[\@max_cluster_rmsd</td><td>&lt;maximum rmsd of a cluster (0.0 = disable)&gt;]</td></tr>
 * <tr><td>[\@grid_size</td><td>&lt;size of the grid for storage&gt;]</td></tr>
 * <tr><td>[\@reaction_field</td><td>&lt;cutoff, epsilon, kappa&gt;]</td></tr>
 * </table>
 *
 * Example:
 * @verbatim
   cg_potential
     @topo             ex.top
     @pbc              r
     @cutoff           1.4
     @traj             ex.tr
     @cluster_size     4
     @max_cluster_rmsd 0.1
     @grid_size        1000
     @reaction_field   1.4 62.0 0.0
   @endverbatim

 * <hr>
 */
#include <cassert>
#include <cstdlib>
#include <vector>
#include <vector>
#include <map>
#include <vector>
#include <string>
#include <iomanip>
#include <cmath>
#include <iostream>


#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/bound/Boundary.h"
#include "../src/gio/InG96.h"
#include "../src/gio/InTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/AtomPair.h"
#include "../src/gcore/LJType.h"
#include "../src/gcore/System.h"
#include "../src/gcore/GromosForceField.h"
#include "../src/gcore/Solvent.h"
#include "../src/gcore/SolventTopology.h"
#include "../src/gcore/Box.h"
#include "../src/gmath/Vec.h"
#include "../src/gmath/Physics.h"
#include "../src/utils/debug.h"
#include "../src/gromos/Exception.h"

#ifdef OMP
#include <omp.h>
#endif

using namespace gcore;
using namespace gmath;
using namespace gio;
using namespace bound;
using namespace args;
using namespace std;

struct cluster {
  /**
   * centre of geometry of the cluster
   */
  gmath::Vec cog;
  /**
   * rmsd from the cog
   */
  double rmsd;
  /**
   * atoms belonging to the cluster
   */
  std::vector<unsigned int> member_molecules;
};

class SolventEnergy {
public:
  /**
   * constructor
   */
  SolventEnergy(gcore::System &sys, unsigned int sol, gcore::GromosForceField &gff, 
           bound::Boundary &pbc);
  
  /**
   * destructor
   */
  ~SolventEnergy();
  
  /**
   * set the cutoff
   */
  void setCutoff(double c);
  /**
   * set reaction field parameters
   */
  void setRF(double eps, double kap);
  /** 
   * calculate the energy of a pair
   * @param mol_i first atom of molecule 1
   * @param mol_j first atom of molecule 2
   */
  void calcPair(unsigned int mol_i, unsigned int mol_j, double & e_lj, double & e_crf);
protected:
  /**
   * internal function to initialize the reaction field
   */
  void init();
  /**
   * the system on which we calculate everything
   */
  gcore::System & m_sys;
  /**
   * the solvent on which we calculate everything
   */
  unsigned int m_sol;
  /**
   * the force field parameters
   */
  gcore::GromosForceField & m_gff;
  /**
   * the periodic boundary conditions
   */
  bound::Boundary & m_pbc;
  
  /**
   * reaction field constants
   */
  double m_cutoff, m_cut3i, m_eps, m_kap, m_crf, m_crf_cut3i, m_crf_2cut3i, m_crf_cut;
  
  /**
   * @struct param
   * a struct that holds the precalculated interaction parameters
   */
  struct param {
    double q_eps, c6, c12;
  };
  
  /**
   * the precalculated interaction parameters
   */
  param * m_parameters;
  
};

SolventEnergy::SolventEnergy(gcore::System &sys, unsigned int sol,
           gcore::GromosForceField &gff,   bound::Boundary &pbc) :
           m_sys(sys), m_sol(sol), m_gff(gff), m_pbc(pbc), m_cutoff(1.4), m_eps(1.0), m_kap(0.0) {
  // precalculate the parameters
  const unsigned int num_atoms = m_sys.sol(m_sol).topology().numAtoms();
 
  // to parameters for all atom/atom combinations are stored in this array
  const unsigned int num_param = num_atoms * num_atoms;
  m_parameters = new param[num_param];
  unsigned int p = 0;
  // loop over first solvent molecule
  for(unsigned int i = 0; i < num_atoms; ++i) {
    // loop over second solvent molecule
    for(unsigned int j = 0; j < num_atoms; ++j, ++p) {
      // precalculate charge product
      m_parameters[p].q_eps = physConst.get_four_pi_eps_i() *
          m_sys.sol(m_sol).topology().atom(i).charge() *
          m_sys.sol(m_sol).topology().atom(j).charge();
      // cache lennard-jones parameters
      LJType lj(m_gff.ljType(
          AtomPair(m_sys.sol(m_sol).topology().atom(i).iac(),
                   m_sys.sol(m_sol).topology().atom(j).iac())));
      m_parameters[p].c6 = lj.c6();
      m_parameters[p].c12 = lj.c12();
      DEBUG(3, "P: " << p << " q: " << m_parameters[p].q_eps << " C6: " <<
              m_parameters[p].c6 << " C12: " << m_parameters[p].c12);
    }
  }
  init();
}

SolventEnergy::~SolventEnergy() {
  // free the memory
  delete[] m_parameters;
}

void SolventEnergy::setCutoff(double c) {
  m_cutoff = c;
  init();
}

void SolventEnergy::setRF(double eps, double kap) {
  m_eps = eps; m_kap = kap;
  init();
}

void SolventEnergy::init() {
  // prepare all the parameters needed for 
  // reaction field calculation
  m_cut3i = 1.0 / (m_cutoff * m_cutoff * m_cutoff);
  const double epsilon = 1.0;
  m_crf = 2.0 * (epsilon - m_eps) * (1.0 + m_kap * m_cutoff) -
      m_eps * (m_kap * m_cutoff * m_kap * m_cutoff);

  m_crf /= (epsilon + 2.0 * m_eps) *
      (1.0 + m_kap * m_cutoff) + m_eps * (m_kap * m_cutoff * m_kap * m_cutoff);
  m_crf_cut3i = m_crf * m_cut3i;
  m_crf_2cut3i = m_crf_cut3i / 2.0;
  m_crf_cut = (1 - m_crf / 2.0) / m_cutoff;
}

void SolventEnergy::calcPair(unsigned int mol_i, unsigned int mol_j, double & e_lj, double & e_crf) {
  const unsigned int num_atoms = m_sys.sol(m_sol).topology().numAtoms();
  
  e_crf = 0.0; e_lj = 0.0;
  // calculate the crf and lj energy between two solvent molecules
  
  // p is the index for the precalculated parameter table
  unsigned int p = 0;
  // loop over atoms in the first molecule
  for(unsigned int i = 0; i < num_atoms; ++i) {
    const gmath::Vec & r_i = m_sys.sol(m_sol).pos(mol_i + i);
    // loop over atoms in the second molecule
    for(unsigned int j = 0; j < num_atoms; ++j, ++p) {
      // calculate the distance
      const double dist2 = (r_i - m_pbc.nearestImage(r_i,
          m_sys.sol(m_sol).pos(mol_j + j), m_sys.box())).abs2();
      
      assert(dist2 != 0.0);
      // calculate the energy
      const double dist2i = 1.0 / dist2;
      const double dist6i = dist2i * dist2i * dist2i;
      const double disti = std::sqrt(dist2i);

      e_lj += (m_parameters[p].c12 * dist6i - m_parameters[p].c6) * dist6i;
      e_crf += m_parameters[p].q_eps * (disti - m_crf_2cut3i * dist2 - m_crf_cut);
    }
  }
}


/**
 * @struct grid_point_struct
 * grid point struct to hold the data needed to calculate the mean and rmsd
 */
struct grid_point_struct {
  /**
   * number of data points collected
   * as it can be huge it should be a big integer type
   */
  unsigned long long int n;
  /**
   * the sum of the data
   */
  double sum;
  /**
   * the sum of the squared data
   */
  double sum_squares;
  /**
   * constructor, sets sums and n to zero
   */
  grid_point_struct() : n(0), sum(0.0), sum_squares(0.0) {
  }
};

int main(int argc, char **argv) {

  Argument_List knowns;
  knowns << "topo" << "pbc" << "cutoff" << "traj" << "cluster_size" << "grid_size"
      << "reaction_field" << "max_cluster_rmsd";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo   <topology>\n";
  usage += "\t@pbc    <boundary type>\n";
  usage += "\t@cutoff <distance within which atoms are bound>\n";
  usage += "\t@traj   <trajectory files>\n";
  usage += "\t@cluster_size 4\n";
  usage += "\t@max_cluster_rmsd 0.42\n";
  usage += "\t@grid_size 2000\n";
  usage += "\t@reaction_field cutoff epsilon kappa\n";

  try {
    Arguments args(argc, argv, knowns, usage);

    // read topology
    InTopology it(args["topo"]);
    System sys(it.system());
    GromosForceField gff(it.forceField());

    // Parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(sys, args);

    // read in cutoff
    double cutoff = args.getValue<double>("cutoff", false, 1.4);
    const double cutoff2 = cutoff*cutoff;
    // read in cluster size
    unsigned int cluster_size = args.getValue<unsigned int>("cluster_size", false, 4);
    // read in maxium cluster rmsd
    double max_cluster_rmsd = args.getValue<double>("max_cluster_rmsd", false, 0.0);
    // read in grid size
    unsigned int grid_size = args.getValue<unsigned int>("grid_size", false, 1000);

    // this variable is used to map a distance r to a grid point by
    // unsigned int point_index = to_grid * r;
    const double to_grid = grid_size / cutoff;
    
    // create the grid as a struct holding the information needed to calculate
    // the mean and the rmsd
    vector<grid_point_struct> grid(grid_size);

    // read the reaction field parameters
    vector<double> rf_param = args.getValues<double>("reaction_field", 3, false,
            Arguments::Default<double>() << 1.4 << 1.0 << 0.0);

    // setup everything for energy calculation
    SolventEnergy ener(sys, 0, gff, *pbc);
    ener.setCutoff(rf_param[0]);
    ener.setRF(rf_param[1], rf_param[2]);

    // define input coordinate
    InG96 ic;

    // loop over all trajectories
    for (Arguments::const_iterator
      iter = args.lower_bound("traj"),
        to = args.upper_bound("traj");
        iter != to; ++iter) {

      // open file
      ic.open((iter->second).c_str());
      
      // select ALL atoms to get the solvent as well
      ic.select("ALL");

      // loop over single trajectory
      while (!ic.eof()) {
        ic >> sys;

        const unsigned int num_atoms_per_molecule = sys.sol(0).topology().numAtoms();
        const unsigned int num_solvent_atoms = sys.sol(0).numAtoms();
        
        DEBUG(1, "num solvent atoms: " << num_solvent_atoms);

        // first we have to determine the clusters and add them to a vector
        // for this we determine the cluster_size closest neighbor molecules
        // and store the first atom of each molecule
        vector<cluster> clusters;
        
        // let's do the whole fun in parallel if possible.
        // this will give different results because the summation
        // is numerically rather sensitive
        #ifdef OMP
        #pragma omp parallel
        #endif
        {   
#ifdef OMP
          const unsigned int tid = omp_get_thread_num();
	  const unsigned int size  = omp_get_num_threads();
#else
          const unsigned int tid = 0;
          const unsigned int size = 1;
#endif 
          // loop over solvent molecules
          // thread 0 starts with the first molecule
          // thread 1 with the second molecule and thus has to calculate the 
          // first atom of the second molecule.
          //
          // in general the starting point is given by
          // i = tid * num_atoms_per_molecule
          // of course we have to jump by the number of thread molecules at
          // every step and calculate the proper offset of the first atom
          for (unsigned int i = tid * num_atoms_per_molecule;
               i < num_solvent_atoms; i += num_atoms_per_molecule * size) {
            const gmath::Vec & pos_i = sys.sol(0).pos(i);

            // here we use a map because with the distance as keys because the
            // list is then sorted by key -> we get the sorted neigbours on fly
            // maybe there is a faster implementation with a std::pair and a
            // user definied compare operator. But this works fine.
            map<double, unsigned int> neighbours;

            // loop over solvent molecules again
            for (unsigned int j = 0; j < num_solvent_atoms; j += num_atoms_per_molecule) {
              // don't include i in the cluster calculation
              if (i == j)
                continue;

              // get the distances and save the cluster
              const gmath::Vec & pos_j = sys.sol(0).pos(j);
              const double r2 = (pos_i - pbc->nearestImage(pos_i, pos_j, sys.box())).abs2();
              // this is a bit dangerous as r2 could already exist. But this
              // is very unlikely: r2 needs to be exactly the same and in realistic
              // equilibrated boxes this is very unlikley.
              neighbours[r2] = j;
            }

            // only create a cluster if the water actually has at least 
            // cluster size neighbors
            if (neighbours.size() >= cluster_size) {
              // add the solvent molecule i to the cluster
              cluster c;
              c.cog = gmath::Vec(0.0, 0.0, 0.0);
              c.rmsd = 0.0;
              // this is incremented whenever an atom is added to the COG.
              unsigned int num_atoms = 0;
            
              // add the molecule itself to the cluster. It is the centre.
              c.member_molecules.push_back(i);
              for (unsigned int atom = 0; atom < num_atoms_per_molecule; ++atom, ++num_atoms) {
                // add positions to the cog
                c.cog += sys.sol(0).pos(i + atom);
              }
            
              // here we have to make sure that the gathering is correct. We 
              // gather the cluster to the first atom to calculate the cog
              const gmath::Vec & first_atom = sys.sol(0).pos(i);

              // loop over cluster_size closest neighbours.
              map<double, unsigned int>::const_iterator it = neighbours.begin();
              for (unsigned int mol = 0; mol < cluster_size; ++it, ++mol) {
                // add it to the cluster
                c.member_molecules.push_back(it->second);
                // loop over the molecule' atoms to calculate the cog
                for (unsigned int atom = 0; atom < num_atoms_per_molecule; ++atom, ++num_atoms) {
                  // gather to the first atom
                  gmath::Vec r = pbc->nearestImage(first_atom, sys.sol(0).pos(it->second + atom), sys.box());
                  // and add it to the cog
                  c.cog += r;
                }
              } // for neighbours

              // calculate the cog
              c.cog /= num_atoms;
              DEBUG(3, "mols: " << c.member_molecules[0] << ","
                  << c.member_molecules[1] << ","
                  << c.member_molecules[2] << ","
                  << c.member_molecules[3]);
              DEBUG(3, "cog:" << gmath::v2s(c.cog));

              if (max_cluster_rmsd) {
                // now we have the cog, let's calculate the rmsd from it
                // rmsd from molecule i
                for (unsigned int atom = 0; atom < num_atoms_per_molecule; ++atom) {
                  // calculate the squared distance from the COG.
                  c.rmsd += (c.cog - pbc->nearestImage(c.cog, sys.sol(0).pos(i + atom), sys.box())).abs2();
                }

                // rmsd from neighbours
                it = neighbours.begin();
                for (unsigned int mol = 0; mol < cluster_size; ++it, ++mol) {
                  for (unsigned int atom = 0; atom < num_atoms_per_molecule; ++atom) {
                    // calculate the squared distance from the COG.
                    c.rmsd += (c.cog - pbc->nearestImage(c.cog, sys.sol(0).pos(it->second + atom), sys.box())).abs2();
                  } // for atoms
                } // for neighbours

                // calculate the rmsd
                c.rmsd = sqrt(c.rmsd / num_atoms);
                DEBUG(3, "rmsd: " << c.rmsd);
               
                // decide whether the cluster fullfills the quality criterion
              if (c.rmsd <= max_cluster_rmsd) {
                  #ifdef OMP
                  #pragma omp critical
                  #endif
                  clusters.push_back(c);
                }
              } else {
                // add it to the clusters ignoring rmsd
                #ifdef OMP
                #pragma omp critical
                #endif
                clusters.push_back(c);
              }
            } // if has neigbours
          } // loop over solvent

          // here we have to synchronize the threads. First the clustering
          // needs to be done.
          #ifdef OMP
          #pragma omp barrier
          #endif
          
          // now we have the clusters and can calculate the energy between two
          // clusters
          DEBUG(1, "found " << clusters.size() << " clusters.");
        
          const unsigned int num_cluster = clusters.size();
          // loop over the clusters

          // here the thread stores it's grid. It is to slow to write to the
          // "grid" variable at every step because this is a critical
          // section. We just create a grid for every thread and copy it
          // to the "grid" variable at the end of the energy calculation
          
          vector<grid_point_struct> my_grid(grid_size);
          
          // loop over the clusters
          // this is again parallelized as we start at tid and the step
          // size is the number of threads.
          for (unsigned int i = tid; i < num_cluster; i += size) {
            // loop over the other clusters
            for (unsigned int j = i + 1; j < num_cluster; ++j) {
              // calculate the inter COG distance
              const double r2 = (clusters[i].cog - 
                  pbc->nearestImage(clusters[i].cog, clusters[j].cog, sys.box())).abs2();
              
              // check for range
              if (r2 <= cutoff2) {
                // here we have to guard for overlapping clusters
                bool overlap = false;
                for(unsigned int mol_i = 0; mol_i < cluster_size; ++mol_i) {
                  for (unsigned int mol_j = mol_i + 1; mol_j < cluster_size; ++mol_j) {
                    if (clusters[i].member_molecules[mol_i] ==
                        clusters[j].member_molecules[mol_j]) {
                      // clusters i and j share a molecule: overlap!
                      overlap = true;
                      break;
                    }
                  }
                }
              
                // don't calculate if they overlap
                if (overlap)
                  break;
              
                // calculate the distance and the grid point
                const double r = std::sqrt(r2);
                DEBUG(3, "clusters " << i << "," << j << " r :" << r);
                const unsigned int grid_point = int(to_grid * r);
                DEBUG(3, "grid point: " << grid_point);
                assert(grid_point < grid_size); // this should always be true!

                // calculate the energy for this grid point.
                double energy = 0.0;
                // loop over member molecules of first cluster
                for (unsigned int mol_i = 0; mol_i < cluster_size; ++mol_i) {
                  // loop over member molecules of the second cluster
                  for (unsigned int mol_j = mol_i + 1; mol_j < cluster_size; ++mol_j) {
                    // calculate the energy between the molecules
                    double e_lj, e_crf;
                    ener.calcPair(clusters[i].member_molecules[mol_i],
                                  clusters[j].member_molecules[mol_j], e_lj, e_crf);
                    energy += e_lj + e_crf;
                  }
                }
                // now we have the energy. Add it to the grid point!
                DEBUG(3, "energy: " << energy);
                ++my_grid[grid_point].n; // increment number of data collected
                my_grid[grid_point].sum += energy; // add the energy to the average
                my_grid[grid_point].sum_squares += energy * energy; // add the energy squared for the rmsd
              } // if cutoff
            } // cluster j
          } // cluster i
          
          
          // now we are done and have to copy back the grid. This is
          // a critical section (every thread executes it after the preceding
          // thread has finished.
          #ifdef OMP
          #pragma omp critical
          #endif
          {
            for (unsigned int i = 0; i < grid_size; i++) {
              grid[i].n += my_grid[i].n;
              grid[i].sum += my_grid[i].sum;
              grid[i].sum_squares += my_grid[i].sum_squares;
            }// copy grid
          } // critical
        } // omp parallel section
      } // loop over frames
      ic.close();

    } // loop over trajectory files

    // write the grid to the output
    cout.setf(std::ios::scientific, std::ios::floatfield);
    cout.precision(8);

    // title line
    cout << "#" << setw(19) << "r" << setw(20) << "energy" << setw(20)
        << "rmsd" << setw(20) << "error estimate" << endl
        << "#" << endl;

    // loop over the grid points
    for (unsigned int i = 0; i < grid_size; ++i) {
      // calculate  the grid point in nm
      const double r = double(i) / to_grid;
      
      // if data was collected
      if (grid[i].n) {
        // calculate average and rmsd
        const double ave = grid[i].sum / grid[i].n;
        const double ssum = grid[i].sum_squares / grid[i].n;
        const double rmsd = sqrt(ssum - ave*ave);
        
        // write them to the output
        cout << setw(20) << r << setw(20) << ave
             << setw(20) << rmsd << endl;
      } else {
        // no data collected: write a comment.
        cout << "#" << setw(19) << r << " no sampling" << endl;
      }
    }
  } catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}
