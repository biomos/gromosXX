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
 * The BSLEUSCOORD specifies the (internal) coordinates defining the active subspace
 * @verbatim
BSLEUSCOORD
#
# Define the possible coordinates for the B&S-LEUS algorithm
# NCDIM:    Number of coordinate dimensions
# CID:      ID of the coordinate dimension
# CTYPE:    Type of the coordinate
#   1:  Dihedral angles (i, j, k, l)
#   ...
# CREF:     The reference value of the coordinate (sigma_n)
#           Only important for the radius and the width,
#           since they should be normed and dimensionless.
# DEF:      The definitions of the coordinate
#
# NCDIM
  2
# CID   CTYPE   CREF    DEF
  1     1       1       2 4 6 8
  2     1       1       4 6 8 10
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
# SPHID SPHCEN      SPHRAD  SPHCLE  SPHNGP
  1     60   300    45      0.5     10
  2     300  300    45      0.5     10
  3     300  60     45      0.5     10
  4     60   60     45      0.5     10
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
# STKID PRSTYP  STKSTART    STKEND      STKWD   STKCLE  STKNGP
# or
# STKID PRSTYP  [SRTCOORD]  [ENDCOORD]  STKEND  STKWD   STKCLE  STKNGP
  1     0       1           3           10      0.5     20
  2     0       2           4           10      0.5     20
END
@endverbatim
 */
void io::In_BSLEUS::read(topology::Topology &topo,
          configuration::Configuration &conf,
	      simulation::Simulation & sim,
	      std::ostream & os)
{
  std::vector<std::string> buffer;
  util::BS_Subspace *bs_subspace = new util::BS_Subspace();
  std::vector<double> periodicities;
  // Number of coordinates
  int numCoords;

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

    int id, type, i, j, k, l, numCoordRead = 0;
    double reference;


    // Loop over the coordinates
    for (; it != to; it++) {

      _lineStream.clear();
      _lineStream.str(*it);
      _lineStream >> id >> type >> reference;
      if (_lineStream.fail()) {
        io::messages.add("bad line in BSLEUSCOORD block",
                "In_BSLEUS", io::message::error);
        return;
      }

      switch (type) {
        case util::BS_Coordinate::dihedral: {
          _lineStream >> i >> j >> k >> l;
          if (_lineStream.fail()) {
            io::messages.add("bad line in BSLEUSCOORD block",
                    "In_BSLEUS", io::message::error);
            return;
          }
          // Convert to gromos
          i--, j--, k--, l--;
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
          periodicities.push_back(dih->getPeriodicity());
          DEBUG(10, dih->str());
          bs_subspace->addCoordinate(dih);
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
    
    int numSpheres, numSpheresRead = 0;
    _lineStream >> numSpheres;
    if (_lineStream.fail()) {
      io::messages.add("Couldn't get the number of Coordinates in BSLEUSSPH",
              "In_BSLEUS", io::message::error);
      return;
    }
    
    int id = 0, sumID = 0, num_gp;
    double radius, forceConst, coord;
    std::vector<double> centerValues;
    util::BS_Vector center;
    for (; it != to; it++){
      _lineStream.clear();
      _lineStream.str(*it);
      _lineStream >> id;
      for (int i = 0; i < numCoords; i++){
        _lineStream >> coord;
        centerValues.push_back(coord);
      }
      _lineStream >> radius >> forceConst >> num_gp;
      if (_lineStream.fail()) {
        io::messages.add("bad line in BSLEUSSPH block",
                "In_BSLEUS", io::message::error);
        return;
      }
      center.create(centerValues, periodicities);
      DEBUG(10, "Center: " + center.str())
      util::BS_Sphere *bs_sphere = new util::BS_Sphere(id, num_gp, forceConst, 
                                    sim.param().bsleus.forceConstantIncrement,
                                                        center, radius);
      bs_subspace->addPotential(bs_sphere);
      numSpheresRead++;
      sumID += id;
      centerValues.clear();
    }
    DEBUG(10, "Finished Reading in Spheres");
    if (numSpheres != numSpheresRead) {
      io::messages.add("The numbers of spheres in BSLEUSSPH is wrong!",
              "In_BSLEUS", io::message::error);
      return;
    }    
    if (sumID != (id * (id + 1) / 2)){
      io::messages.add("The IDs of the spheres are not given in a consecutive order!",
              "In_BSLEUS", io::message::error);
      io::messages.add("The Subspace has not been added!",
              "In_BSLEUS", io::message::error);
      return;
    }
  } // BSLEUSSPH
  
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
    
    int numSticks, numSticksRead = 0;
    _lineStream >> numSticks;
    if (_lineStream.fail()) {
      io::messages.add("Couldn't get the number of Coordinates in BSLEUSSTK",
              "In_BSLEUS", io::message::error);
      return;
    }
    
    int id = 0, sumID = 0, num_gp, defType, startSphere, endSphere;
    double width, forceConst;
    util::BS_Vector start, end;
    for (; it != to; it++){
      _lineStream.clear();
      _lineStream.str(*it);
      _lineStream >> id >> defType;
      if (_lineStream.fail()) {
        io::messages.add("bad line in BSLEUSSTK block",
                "In_BSLEUS", io::message::error);
        return;
      }
      switch (defType) {
        case 0: {// define start and end points in terms of sphere ids
          _lineStream >> startSphere >> endSphere;
          if (_lineStream.fail()) {
            io::messages.add("bad line in BSLEUSSTK block",
                    "In_BSLEUS", io::message::error);
            return;
          }
          start = bs_subspace->getCenter(startSphere);
          end = bs_subspace->getCenter(endSphere);
          break;
        }
        case 1: {// define end and start points with coordinates
          std::vector<double> coords, periodicities;
          bs_subspace->getPeriodicities(periodicities);
          for (int j = 0; j < 2; j++) {
            coords.clear();
            for (int i = 0; i < numCoords; i++) {
              double value;
              _lineStream >> value;
              if (_lineStream.fail()) {
                io::messages.add("bad line in BSLEUSSTK block",
                        "In_BSLEUS", io::message::error);
                return;
              }
              coords.push_back(value);
            }
            if (j == 0){
              start.create(coords, periodicities);
            } else {
              end.create(coords, periodicities);
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
                                    sim.param().bsleus.forceConstantIncrement,
                                                        start, end, width);
      bs_subspace->addPotential(bs_stick);
      numSticksRead++;
      sumID += id;
    }
    DEBUG(5, "The number of sticks according to file: " << numSticks << "; actually read: " << numSticksRead);
    if (numSticks != numSticksRead) {
      io::messages.add("The numbers of Sticks in BSLEUSSTK seems wrong!",
              "In_BSLEUS", io::message::warning);
      return;
    }
    if (sumID != (id * (id + 1) / 2)){
      io::messages.add("The IDs of the sticks are not given in a consecutive order!",
              "In_BSLEUS", io::message::warning);
      io::messages.add("The Subspace has not been added!",
              "In_BSLEUS", io::message::warning);
      return;
    }
  } // BSLEUSSTK
  conf.special().bs_umbrella.addSubspace(bs_subspace);
}
