/*
 * This file is part of GROMOS.
 *
 * Copyright (c) 2011, 2012, 2016, 2018, 2021, 2023 Biomos b.v.
 * See <https://www.gromos.net> for details.
 */

/**
 * @file in_colvarres.cc
 * implements methods of In_Topology for COLVAR geometry specification files.
 */

#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../interaction/interaction_types.h"
#include "../../configuration/configuration.h"

#include "../../io/instream.h"
#include "../../io/blockinput.h"

#include "in_colvarres.h"

#undef MODULE
#undef SUBMODULE
#define MODULE io
#define SUBMODULE topology

#include "in_colvarres_helpers.h"

using namespace io::colvar_input;

/**
 * @section colvarresfile COLVAR restraint specification file
 *
 * The RESTRAINTS block stores the bias parameters. Geometry is read from the
 * matching *RESSPEC block. If a non-zero index is present, the geometry row is
 * matched to the RESTRAINTS row with the same index. Old files without indices
 * are still matched by order within each TYPE.
 *
 * TYPE values:
 * 0 DISTANCE, 1 ANGLE, 2 DIHEDRAL, 3 COORDNUM, 4 DISTANCE_DIFFERENCE
 *
 * DISTANCE_DIFFERENCE is |A-B| - |C-D| in distance units. It requires a
 * positive INDX and four virtual-atom endpoints in DISTANCEDIFFRESSPEC (or
 * PERTDISTANCEDIFFRESSPEC). RAH is restricted to -1/0/1. METHOD 0, 1, and 11
 * are supported; inverse-cubic averaging is not defined for a signed CV.
 *
 * ANGLE and DIHEDRAL TARGET/CVLIN values are read in degrees and converted
 * internally to radians. RAH -1/0/1 selects a lower wall/full harmonic/upper
 * wall for every CV. Only DISTANCE accepts dimensional projection codes
 * 10..60, with -1/0/1 added for wall behavior. ANGLE/DIHEDRAL walls require
 * complementary -1/+1 rows on the same ordered geometry with matching bias
 * settings. Dihedral pairs define a directed allowed arc and share an image
 * about its midpoint; they currently require METHOD 0.
 * METHOD encodes restraining/averaging behavior:
 * -1 colvar constraint (recognised but not implemented yet),
 * 0 no averaging, 1 exponential average, 2 inverse-cubic distance average,
 * 11/12 with relaxation force scaling, 22 with chain-rule force scaling.
 *
 * @verbatim
TITLE
collective variable restraints
END

RESTRAINTS
# NRES
4
# INDX TYPE METHOD TARGET CVK RAH TAU VIRIAL CVLIN
  1    0      0    0.24   10000   0   0.0   1   0.0
  2    1      0     104     490   0   0.0   0   0.0
  3    2      0      -5   15000   0   0.0   0   0.0
  4    3      0       1   10000   0   0.0   1   0.0
END

DISTANCERESSPEC
# DISH DISC
  0.1  0.153
# INDX i j k l TYPE   i j k l TYPE
  1    95 0 0 0 0     14 0 0 0 0
END

ANGRESSPEC
# INDX IPLR JPLR KPLR
  2    14   95   23
END

DIHEDRALRESSPEC
# INDX IPLR JPLR KPLR LPLR
  3    23   24   25   95
END

COORDNUMRESSPEC
# INDX RCUT N M
  4    0.3  20 40
# TYPE GROUP ATOM
  0    1     95
  0    2     96
  0    2     99
END
@endverbatim
 *
 * Perturbed restraints use PERTRESTRAINTS plus matching PERT*RESSPEC geometry
 * blocks. Each PERTRESTRAINTS row defines one perturbed term with A and B
 * state bias parameters. Geometry is not perturbed and can be matched by index
 * in the same way as normal restraints.
 *
 * @verbatim
PERTRESTRAINTS
# NRES
4
# INDX TYPE METHOD TARGETA CVKA TARGETB CVKB RAH TAU VIRIAL CVLIN [N M]
  5    0      0    0.24  10000   0.30  20000   0   0.0   1   0.0
  6    1      0     104    490    109    700   0   0.0   0   0.0
  7    2      0      -5  15000     10   9000   0   0.0   0   0.0
  8    3      0       1  10000      2  12000   0   0.0   1   0.0
END

PERTDISTANCERESSPEC
# DISH DISC
  0.1  0.153
# INDX i j k l TYPE   i j k l TYPE
  5    95 0 0 0 0     14 0 0 0 0
END

PERTANGRESSPEC
# INDX IPLR JPLR KPLR
  6    14   95   23
END

PERTDIHEDRALRESSPEC
# INDX IPLR JPLR KPLR LPLR
  7    23   24   25   95
END

PERTCOORDNUMRESSPEC
# INDX RCUT N M
  8    0.3  20 40
# TYPE GROUP ATOM
  0    1     95
  0    2     96
  0    2     99
END
@endverbatim
 */
void io::In_Colvarres::read(topology::Topology& topo,
                            simulation::Simulation & sim,
                            std::ostream & os) {

  DEBUG(7, "reading in a colvar restraints file");
  std::set<std::string> &block_read = m_colvar_blocks_read;
  block_read.clear();
  sim.param().colvarres.bias_specs.clear();
  sim.param().colvarres.pert_bias_specs.clear();
  sim.param().colvarres.eds_bias_specs.clear();
  topo.eds_distance_restraints().clear();

  if (!quiet)
    os << "COLVAR RESTRAINTS\n";


  read_restraints(topo, sim, os);
  read_pertrestraints(topo, sim, os);
  read_distanceresspec(topo, sim, os);
  read_angresspec(topo, sim, os);
  read_dihedralresspec(topo, sim, os);
  read_coordnumresspec(topo, sim, os);
  read_pertdistanceresspec(topo, sim, os);
  read_pertangresspec(topo, sim, os);
  read_pertdihedralresspec(topo, sim, os);
  read_pertcoordnumresspec(topo, sim, os);

  read_eds_colvar_blocks(m_block,topo,sim,block_read);
  read_distance_difference_geometries(m_block["DISTANCEDIFFRESSPEC"],
      "DISTANCEDIFFRESSPEC", sim.param().colvarres.bias_specs, false, topo, block_read);
  read_distance_difference_geometries(m_block["PERTDISTANCEDIFFRESSPEC"],
      "PERTDISTANCEDIFFRESSPEC", sim.param().colvarres.pert_bias_specs, true, topo, block_read);
  read_colvar_branches(m_block, sim, block_read, os, quiet);
  pair_angular_walls(sim.param().colvarres.bias_specs,
                    sim.param().colvarres.bias_specs, topo, false);
  pair_angular_walls(sim.param().colvarres.pert_bias_specs,
                    sim.param().colvarres.bias_specs, topo, true);

  for (std::map<std::string, std::vector<std::string> >::const_iterator
       it = m_block.begin(), to = m_block.end(); it != to; ++it) {
    if (block_read.count(it->first) == 0 && it->second.size()) {
      io::messages.add(" block " + it->first + " not supported!",
                       "In_Colvarres", io::message::warning);
    }
  }

  const unsigned int n_distance_specs = count_colvar_specs(sim, "DISTANCE");
  const unsigned int n_pertdistance_specs =
    count_colvar_specs(sim.param().colvarres.pert_bias_specs, "DISTANCE") / 2;
  const unsigned int n_distance_specs_total =
    n_distance_specs + n_pertdistance_specs;
  if (n_distance_specs_total > 0
      && n_distance_specs_total != topo.distance_restraints().size()) {
    std::ostringstream msg;
    msg << "RESTRAINTS+PERTRESTRAINTS/DISTANCERESSPEC+PERTDISTANCERESSPEC mismatch: "
        << n_distance_specs_total
        << " DISTANCE bias specification(s), but "
        << topo.distance_restraints().size()
        << " DISTANCE geometry specification(s).";
    io::messages.add(msg.str(), "In_Colvarres", io::message::warning);
  }

  const unsigned int n_coord_specs = count_colvar_specs(sim, "COORDNUM");
  const unsigned int n_pertcoord_specs =
    count_colvar_specs(sim.param().colvarres.pert_bias_specs, "COORDNUM") / 2;
  const unsigned int n_coord_specs_total =
    n_coord_specs + n_pertcoord_specs;
  if (n_coord_specs_total > 0
      && n_coord_specs_total != topo.coordnum_restraint().size()) {
    std::ostringstream msg;
    msg << "RESTRAINTS+PERTRESTRAINTS/COORDNUMRESSPEC+PERTCOORDNUMRESSPEC mismatch: "
        << n_coord_specs_total
        << " COORDNUM bias specification(s), but "
        << topo.coordnum_restraint().size()
        << " COORDNUM geometry specification(s).";
    io::messages.add(msg.str(), "In_Colvarres", io::message::warning);
  }

  const unsigned int n_angle_specs = count_colvar_specs(sim, "ANGLE");
  const unsigned int n_pertangle_specs =
    count_colvar_specs(sim.param().colvarres.pert_bias_specs, "ANGLE") / 2;
  const unsigned int n_angle_specs_total =
    n_angle_specs + n_pertangle_specs;
  if (n_angle_specs_total > 0
      && n_angle_specs_total != topo.colvar_angle_geometries().size()) {
    std::ostringstream msg;
    msg << "RESTRAINTS+PERTRESTRAINTS/ANGRESSPEC+PERTANGRESSPEC mismatch: "
        << n_angle_specs_total
        << " ANGLE bias specification(s), but "
        << topo.colvar_angle_geometries().size()
        << " ANGLE geometry specification(s).";
    io::messages.add(msg.str(), "In_Colvarres", io::message::warning);
  }

  const unsigned int n_dihedral_specs = count_colvar_specs(sim, "DIHEDRAL");
  const unsigned int n_pertdihedral_specs =
    count_colvar_specs(sim.param().colvarres.pert_bias_specs, "DIHEDRAL") / 2;
  const unsigned int n_dihedral_specs_total =
    n_dihedral_specs + n_pertdihedral_specs;
  if (n_dihedral_specs_total > 0
      && n_dihedral_specs_total != topo.colvar_dihedral_geometries().size()) {
    std::ostringstream msg;
    msg << "RESTRAINTS+PERTRESTRAINTS/DIHEDRALRESSPEC+PERTDIHEDRALRESSPEC mismatch: "
        << n_dihedral_specs_total
        << " DIHEDRAL bias specification(s), but "
        << topo.colvar_dihedral_geometries().size()
        << " DIHEDRAL geometry specification(s).";
    io::messages.add(msg.str(), "In_Colvarres", io::message::warning);
  }
}
