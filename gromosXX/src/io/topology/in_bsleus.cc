/**
 * @file bsleus.cc
 * implements methods of In_BSLEUS
 */

#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../simulation/parameter.h"
#include "../../interaction/interaction_types.h"
#include "../../configuration/configuration.h"

#include "../../io/instream.h"
#include "../../io/blockinput.h"
#include "../../io/configuration/in_configuration.h"

#include "../../util/bs_coordinate.h"
#include "../../util/bs_potentials.h"
#include "../../util/bs_subspace.h"
#include "../../util/bs_umbrella.h"
#include "../../util/bs_vector.h"
#include "in_bsleus.h"

#undef MODULE
#undef SUBMODULE
#define MODULE io
#define SUBMODULE topology

/**
 * @section BSLEUS input
 * 
 * The input is define in a separate file (@bsleus).
 * 
 * The BSLEUSSUB specifies the subspaces, especially the memory settings
 * @verbatim
BSLEUSSUB
#
# Define the subspaces here. Especially the memory updating
# NUMSPC:   The number of subspaces (currently only one allowed).
# SUBID:    The id of the subspace
# FRCINC:   Basis Force constant increment K_LE
# FRED:     force-constant reduction factor
# LVCO:     The local visiting cutoff
# GVCO:     The global visiting cutoff
#
# NUMSPC
  1
# SUBID FRCINC  FRED  LVCO    GVCO
  1     0.0002  0.5     1       2
END

 @endverbatim
 * 
 * The BSLEUSCOORD specifies the (internal) coordinates defining the active subspace
 * @verbatim
BSLEUSCOORD
#
# Define the possible coordinates for the B&S-LEUS algorithm
# NCDIM:    Number of coordinate dimensions
# CID:      ID of the coordinate dimension
# SUBSP:    Which subspace
# CTYPE:    Type of the coordinate
#   1:  Dihedral angles (i, j, k, l)
#   2:  Distance (i, j)
#   ...
# CREF:     The reference value of the coordinate (sigma_n)
#           Only important for the radius and the width,
#           since they should be normed and dimensionless.
# DEF:      The definitions of the coordinate
#
# NCDIM
  2
# CID   SUBSP   CTYPE   CREF    DEF
  1     1       1       1       2 4 6 8
  2     1       1       1       4 6 8 10
END
@endverbatim
 * 
 * The BSLEUSSPH block defines the spheres
 * @verbatim
BSLEUSSPH
#
# The Balls/Spheres of B&S-LEUS
# NSPH:     Number of Spheres
# SPHID:    Sphere ID
# SPHCEN:   Center of the Sphere (Dim = # of LEUS-Coordinates)
# SPHRAD:   Radius of the Sphere
# SPHCLE:   Force constant (c_LE) of the sphere)
# SPHNGP:   Number of grid points
#
# NSPH
  4
# SPHID SUBSP   SPHCEN      SPHRAD  SPHCLE  SPHNGP
  1     1       60   300    45      0.5     10
  2     1       300  300    45      0.5     10
  3     1       300  60     45      0.5     10
  4     1       60   60     45      0.5     10
END
@endverbatim
 * 
 * The BSLEUSSTK block defines the sticks
 * @verbatim
#
BSLEUSSTK
#
# The Stick of the B&S-LEUS algorithm
# NSTK:     The Number of Sticks
# STKID:    The ID of the Stick
# PRSTYP:   Wheter to use SIDs or coordinates to specify the sticks
#   0:          use SID
#   1:          use coordinates
# STKSTRT:  The SID of the sphere at the beginning
# STKEND:   The SID of the sphere at the end
# STKDW:    The half width of the stick
# STKCLE:   The force constant of the stick
# STKNGP:   The number of grid points of the stick
#
# NSTK
  2
# STKID SUBSP   PRSTYP  STKSTART    STKEND      STKWD   STKCLE  STKNGP
# or
# STKID SUBSP   PRSTYP  [SRTCOORD]  [ENDCOORD]  STKEND  STKWD   STKCLE  STKNGP
  1     1       0       1           3           10      0.5     20
  2     1       0       2           4           10      0.5     20
END
@endverbatim
 */
void io::In_BSLEUS::read(topology::Topology &topo,
          configuration::Configuration &conf,
	      simulation::Simulation & sim,
	      std::ostream & os)
{
  std::vector<std::string> buffer;
  std::vector<util::BS_Subspace *> bs_subspaces;
  int num_subspaces;
  
  { // BSLEUSSUB
    buffer = m_block["BSLEUSSUB"];
    DEBUG(10, "BSLEUSSUB block : " << buffer.size());

    if (!buffer.size()) {
      io::messages.add("no BSLEUSSUB block in B&S-LEUS definition file",
              "In_BSLEUS", io::message::error);
      return;
    }

    std::vector<std::string>::const_iterator it = buffer.begin() + 1,
            to = buffer.end() - 1;
    DEBUG(10, "reading in BSLEUSSUB data");
    _lineStream.clear();
    _lineStream.str(*it++);

    _lineStream >> num_subspaces;
    if (_lineStream.fail()) {
      io::messages.add("Couldn't get the number of Subspaces in BSLEUSSUB!",
              "In_BSLEUS", io::message::error);
      return;
    }
    if (num_subspaces != 1){
      io::messages.add("We currently don't support more than one subspace!",
              "In_BSLEUS", io::message::error);
      return;
    }
    
    double forceIncrement, reductionFactor;
    int id, localCutoff, globalCutoff;
    int last_id = 0;
    for (int i = 0; i < num_subspaces; i++) {
      _lineStream.clear();
      _lineStream.str(*it++);
      _lineStream >> id >> forceIncrement >> reductionFactor
                  >> localCutoff >> globalCutoff;
      if (_lineStream.fail()) {
        io::messages.add("Bad Block in BSLEUSSUB!",
                "In_BSLEUS", io::message::error);
        return;
      }
      util::BS_Subspace *bs_subspace = new util::BS_Subspace(id, forceIncrement, 
              reductionFactor, localCutoff, globalCutoff);
      bs_subspaces.push_back(bs_subspace);
      if (last_id != (id - 1)) {
        io::messages.add("The IDs of the subspaces are not given in a consecutive order!",
                "In_BSLEUS", io::message::error);
        io::messages.add("The Subspace has not been added!",
                "In_BSLEUS", io::message::error);
        return;
      }
    }  
  } // BSLEUSSUB
  // ==============================================================
  // Number of coordinates
  int numCoords;
  std::vector<double> references;

  { // BSLEUSCOORD
    buffer = m_block["BSLEUSCOORD"];
    DEBUG(10, "BSLEUSCOORD block : " << buffer.size());

    if (!buffer.size()) {
      io::messages.add("no BSLEUSCOORD block in B&S-LEUS definition file",
              "In_BSLEUS", io::message::error);
      return;
    }

    std::vector<std::string>::const_iterator it = buffer.begin() + 1,
            to = buffer.end() - 1;

    DEBUG(10, "reading in BSLEUSCOORD data");
    _lineStream.clear();
    _lineStream.str(*it++);

    _lineStream >> numCoords;
    if (_lineStream.fail()) {
      io::messages.add("Couldn't get the number of Coordinates in BSLEUSCOORD",
              "In_BSLEUS", io::message::error);
      return;
    }

    int id, subspace, type, i, j, k, l, numCoordRead = 0;
    double reference;


    // Loop over the coordinates
    for (; it != to; it++) {

      _lineStream.clear();
      _lineStream.str(*it);
      _lineStream >> id >> subspace >> type >> reference;
      if (_lineStream.fail()) {
        io::messages.add("bad line in BSLEUSCOORD block",
                "In_BSLEUS", io::message::error);
        return;
      }
      references.push_back(reference);

      switch (type) {
        case util::BS_Coordinate::dihedral: {
          _lineStream >> i >> j >> k >> l;
          if (_lineStream.fail()) {
            io::messages.add("bad line in BSLEUSCOORD block (Dihedral angle)",
                    "In_BSLEUS", io::message::error);
            return;
          }
          // Convert to gromos
          i--, j--, k--, l--;
          subspace--;
          if (i > topo.num_atoms() ||
                  j > topo.num_atoms() ||
                  k > topo.num_atoms() ||
                  l > topo.num_atoms()) {
            std::ostringstream msg;
            msg << "Dihedral (" << i + 1 << "-" << j + 1 << "-" << k + 1 << "-"
                    << l + 1 << ") atom indices out of range.";
            io::messages.add(msg.str(), "In_BSLEUS", io::message::error);
            return;
          }

          util::BS_Dihedral *dih = new util::BS_Dihedral(id, i, j, k, l, reference);
          DEBUG(10, dih->str());
          bs_subspaces[subspace]->addCoordinate(dih);
          numCoordRead++;
          break;
        }
        case util::BS_Coordinate::distance: {
          _lineStream >> i >> j;
          if (_lineStream.fail()) {
            io::messages.add("bad line in BSLEUSCOORD block (Distance)",
                    "In_BSLEUS", io::message::error);
            return;
          }
          // Convert to gromos
          i--, j--;
          subspace--;
          if (i > topo.num_atoms() || j > topo.num_atoms()) {
            std::ostringstream msg;
            msg << "Distance (" << i + 1 << "-" << j + 1 
                << ") atom indices out of range.";
            io::messages.add(msg.str(), "In_BSLEUS", io::message::error);
            return;
          }
          util::BS_Distance *dst = new util::BS_Distance(id, i, j, reference);
          DEBUG(10, dst->str());
          bs_subspaces[subspace]->addCoordinate(dst);
          numCoordRead++;
          break;
        }
        default:
          io::messages.add("Unknown Type in BSLEUSCOORD Block!",
                  "In_BSLEUS", io::message::error);
          return;
      }
    }

    if (numCoords != numCoordRead) {
      io::messages.add("The numbers of coordinates in BSLEUSCOORD seems wrong!",
              "In_BSLEUS", io::message::warning);
      return;
    }
  } // BSLEUSCOORD
  // ==============================================================
  { // BSLEUSSPH
    buffer = m_block["BSLEUSSPH"];
    DEBUG(10, "BSLEUSSPH block : " << buffer.size());

    if (!buffer.size()) {
      io::messages.add("no BSLEUSSPH block in B&S-LEUS definition file",
              "In_BSLEUS", io::message::error);
      return;
    }

    std::vector<std::string>::const_iterator it = buffer.begin() + 1,
            to = buffer.end() - 1;

    DEBUG(10, "reading in BSLEUSSPH data");
    _lineStream.clear();
    _lineStream.str(*it++);
    
    int subspace, numSpheres, numSpheresRead = 0;
    _lineStream >> numSpheres;
    if (_lineStream.fail()) {
      io::messages.add("Couldn't get the number of Coordinates in BSLEUSSPH",
              "In_BSLEUS", io::message::error);
      return;
    }
    
    int id = 0, last_id = 0, num_gp;
    double radius, forceConst, coord;
    std::vector<double> centerValues;
    util::BS_Vector center;
    for (; it != to; it++){
      _lineStream.clear();
      _lineStream.str(*it);
      _lineStream >> id >> subspace;
      for (int i = 0; i < numCoords; i++){
        _lineStream >> coord;
        if (_lineStream.fail()) {
          io::messages.add("bad center in BSLEUSSPH block",
                  "In_BSLEUS", io::message::error);
          return;
        }
        centerValues.push_back(coord / references[i]);
      }
      _lineStream >> radius >> forceConst >> num_gp;
      if (_lineStream.fail()) {
        io::messages.add("bad line in BSLEUSSPH block",
                "In_BSLEUS", io::message::error);
        return;
      }
      center.create(centerValues);
      DEBUG(10, "Center: " + center.str());
      // Convert to GROMOS
      subspace--;
      util::BS_Sphere *bs_sphere = new util::BS_Sphere(id, num_gp, forceConst, 
                                                        center, radius);
      bs_subspaces[subspace]->addPotential(bs_sphere);
      numSpheresRead++;
      centerValues.clear();
      if (id != (last_id + 1)) {
        io::messages.add("The IDs of the spheres are not given in a consecutive order!",
                "In_BSLEUS", io::message::error);
        io::messages.add("The Subspace has not been added!",
                "In_BSLEUS", io::message::error);
        return;
      }
      last_id = id;
    }
    DEBUG(10, "Finished Reading in Spheres");
    if (numSpheres != numSpheresRead) {
      io::messages.add("The numbers of spheres in BSLEUSSPH is wrong!",
              "In_BSLEUS", io::message::error);
      return;
    }    
    
  } // BSLEUSSPH
  // ==============================================================
  { // BSLEUSSTK
    buffer = m_block["BSLEUSSTK"];
    DEBUG(10, "BSLEUSSTK block : " << buffer.size());

    if (!buffer.size()) {
      io::messages.add("no BSLEUSSTK block in B&S-LEUS definition file",
              "In_BSLEUS", io::message::error);
      return;
    }

    std::vector<std::string>::const_iterator it = buffer.begin() + 1,
            to = buffer.end() - 1;

    DEBUG(10, "reading in BSLEUSSTK data");
    _lineStream.clear();
    _lineStream.str(*it++);
    
    int subspace, numSticks, numSticksRead = 0;
    _lineStream >> numSticks;
    if (_lineStream.fail()) {
      io::messages.add("Couldn't get the number of Coordinates in BSLEUSSTK",
              "In_BSLEUS", io::message::error);
      return;
    }
    
    int id = 0, last_id = 0, num_gp, defType, startSphere, endSphere;
    double width, forceConst;
    util::BS_Vector start, end;
    for (; it != to; it++){
      _lineStream.clear();
      _lineStream.str(*it);
      _lineStream >> id >> subspace >> defType;
      if (_lineStream.fail()) {
        io::messages.add("bad line in BSLEUSSTK block",
                "In_BSLEUS", io::message::error);
        return;
      }
      // Convert to GROMOS
      subspace--;
      switch (defType) {
        case 0: {// define start and end points in terms of sphere ids
          _lineStream >> startSphere >> endSphere;
          if (_lineStream.fail()) {
            io::messages.add("bad line in BSLEUSSTK block",
                    "In_BSLEUS", io::message::error);
            return;
          }
          start = bs_subspaces[subspace]->getCenter(startSphere);
          end = bs_subspaces[subspace]->getCenter(endSphere);
          break;
        }
        case 1: {// define end and start points with coordinates
          std::vector<double> coords;//, periodicities;
          for (int j = 0; j < 2; j++) { // loop over start and end coordinates
            coords.clear();
            for (int i = 0; i < numCoords; i++) {
              double value;
              _lineStream >> value;
              if (_lineStream.fail()) {
                io::messages.add("bad line in BSLEUSSTK block",
                        "In_BSLEUS", io::message::error);
                return;
              }
              coords.push_back(value / references[i]);
            }
            if (j == 0){
              start.create(coords);
            } else {
              end.create(coords);
            }
          }
          break;
        }
        default: {
          io::messages.add("Unkown specifier for coordinate definition of sticks!",
                "In_BSLEUS", io::message::error);
          return;
        }
      }
      DEBUG(10, "Start and end");
      DEBUG(10, start.str());
      DEBUG(10,end.str());
      
      _lineStream >> width >> forceConst >> num_gp;
      if (_lineStream.fail()) {
        io::messages.add("bad line in BSLEUSSTK block",
                "In_BSLEUS", io::message::error);
        return;
      }
      
      util::BS_Stick *bs_stick = new util::BS_Stick(id, num_gp, forceConst,
                                                        start, end, width);
      bs_subspaces[subspace]->addPotential(bs_stick);
      numSticksRead++;
      if (last_id != (id - 1)) {
        io::messages.add("The IDs of the sticks are not given in a consecutive order!",
                "In_BSLEUS", io::message::warning);
        io::messages.add("The Subspace has not been added!",
                "In_BSLEUS", io::message::warning);
        return;
      }
      last_id = id;
    }
    DEBUG(5, "The number of sticks according to file: " << numSticks << "; actually read: " << numSticksRead);
    if (numSticks != numSticksRead) {
      io::messages.add("The numbers of Sticks in BSLEUSSTK seems wrong!",
              "In_BSLEUS", io::message::warning);
      return;
    }
    
  } // BSLEUSSTK
  conf.special().bs_umbrella.addSubspace(bs_subspaces);
}
