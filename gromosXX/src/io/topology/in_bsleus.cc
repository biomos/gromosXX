/**
 * @file in_bsleus.cc
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
#include "../../util/template_split.h"
#include "../../math/periodicity.h"

#include "in_bsleus.h"

#undef MODULE
#undef SUBMODULE
#define MODULE io
#define SUBMODULE topology

/**
 * @section bsleussub BSLEUSSUB block
 * 
 * The parmater input is define elsewhere (@ref bsleus).
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
 * @section bsleuscoord BSLEUSCOORD block
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
#   3:  Cartesian Coordinate (numAtoms, Atoms[1..numAtoms])
#   4:  Sum of two dihedral angles ({i,j,k,l}_1, {i,j,k,l}_2)
#   5:  Lambda Variable
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
#  3     1       2       0.1     1 10
#  4     1       3       1       4  2 4 6 8
#  5     1       4       1       4 6 8 10 6 8 10 12
#  6     1       5       1       
END
@endverbatim
 * 
 * @section bsleussph BSLEUSSPH block
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
 * @section bsleusstk BSLEUSSTK block
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
 * @section bsleussnake BSLEUSSNAKE block
 * 
 * The BSLEUSSNAKE block defines the snake-like potentials
 * @verbatim
BSLEUSSNAKE
#
# A snake-like potential
# NSNK:     The number of snake potentials
# SNKID:    The ID
# SUBSP:    The subspace
# NPNTS:    The number of points
# POINTS:   The Points
# STKHWD:   The halfwidth
# SNKCLE:   The halfharmonic force constant
# SNKNGP:   The number of grid points between two points (every point is one)
#
# NSNK
  1
# SNKID SUBSP   NPNTS   POINTS[1..NPNTS]                    STKHWD  SNKCLE  SNKNGP
  1     1       3       {tp_1.cnf} {tp_2.cnf} {tp_3.cnf}    10      0.5     2
END
@endverbatim
 * 
 * @section bsleuspipe BSLEUSPIPE block
 * 
 * The BSLEUSSNAKE block defines the snake-like potentials
 * @verbatim
BSLEUSPIPE
# A pipe like potential (active subspace only in the shell of a cylinder)
# NPIPES:   The number of pipe-like potentials
# ID:       The id of the potential
# SUBSP:    The subspace
# SPECTYP:  Use either spheres (0) or coordinates (1) to define the end points.
# START:    All the values for the start point
# END:      All the values for the end point
#   The included values are:
#   POINT:      The point of the start/end
#   IHW:        The inner half width
#   OHW:        The outer half width
# CLE:      The half-harmonic restraint potential
# NGPL:     The number of grid points in the longitudinal direction
# NGPP:     The number of grid points in the perpendicular direction
#
# NPIPES
  1
#                           START               END
# ID    SUBSP   SPECTYPE    POINT   IHW   OHW   POINT   IHW   OHW   CLE     NGPL NGPP
  2     1       1           1 1     0.5   1.0   2 2     0.6   1.1   0.05    10   2
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
  int num_subspaces = 0;
  
  { // BSLEUSSUB
    buffer = m_block["BSLEUSSUB"];
    DEBUG(10, "BSLEUSSUB block : " << buffer.size());

    if (!buffer.size()) {
      io::messages.add("no BSLEUSSUB block in B&S-LEUS definition file",
              "In_BSLEUS", io::message::error);
      return;
    }

    std::vector<std::string>::const_iterator it = buffer.begin() + 1;
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
    
    double forceIncrement = 0.0, reductionFactor = 0.0;
    int id = 0, localCutoff = 0, globalCutoff = 0;
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
  int numCoords = 0;
  std::vector<unsigned int> numCoordsPerSubspace(num_subspaces, 0);
  std::vector<std::vector<double> > references(num_subspaces, std::vector<double>());
  std::map<unsigned int, std::string> refFiles;
  std::vector<unsigned int> cartAtoms;

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

    int id = 0, subspace = 0, type = 0, numCoordRead = 0;
    double reference = 0.0;


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
      DEBUG(10, "ID: " << id << " SUBSP: " << subspace << " TYPE: " << type << " REF: " << reference);
      // Convert to GROMOS
      subspace--;
      util::BS_Coordinate *coord = nullptr;
      
      // The Type of the coordinate
      switch (type) {
        case util::BS_Coordinate::dihedral: {
          unsigned int i = 0, j = 0, k = 0, l = 0;
          _lineStream >> i >> j >> k >> l;
          if (_lineStream.fail()) {
            io::messages.add("bad line in BSLEUSCOORD block (Dihedral angle)",
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

          coord = new util::BS_Dihedral(id, i, j, k, l, reference);
          break;
        }
        case util::BS_Coordinate::distance: {
          unsigned int i = 0, j = 0;
          _lineStream >> i >> j;
          if (_lineStream.fail()) {
            io::messages.add("bad line in BSLEUSCOORD block (Distance)",
                    "In_BSLEUS", io::message::error);
            return;
          }
          // Convert to gromos
          i--, j--;
          if (i > topo.num_atoms() || j > topo.num_atoms()) {
            std::ostringstream msg;
            msg << "Distance (" << i + 1 << "-" << j + 1 
                << ") atom indices out of range.";
            io::messages.add(msg.str(), "In_BSLEUS", io::message::error);
            return;
          }
          coord = new util::BS_Distance(id, i, j, reference);
          break;
        }
        case util::BS_Coordinate::cartesian: {
          int num_atoms_int = 0;
          unsigned int num_atoms = 0;
          _lineStream >> num_atoms_int;
          if (_lineStream.fail()) {
            io::messages.add("bad line in BSLEUSCOORD block (Cartesian)",
                    "In_BSLEUS", io::message::error);
            return;
          }
          // workaround to check for negative input, but also make sure 
          // num_atoms is an unsigned int
          if (num_atoms_int < 0){
            io::messages.add("Cartesian: Number of Atoms is negative!",
                    "In_BSLEUS", io::message::error);
            return;
          } else{
            num_atoms = (unsigned int) num_atoms_int;
          }
          bool allAtoms = true; // num_atoms == 0
          std::vector<unsigned int> atom_list;
          if (num_atoms){
            allAtoms = false;
            unsigned int atom = 0;
            DEBUG(10, "We will add in total " << num_atoms << " atoms.")
            for (unsigned int i = 0; i < num_atoms; i++){
              _lineStream >> atom;
              if (_lineStream.fail()) {
                io::messages.add("bad atoms in BSLEUSCOORD block (Cartesian)",
                        "In_BSLEUS", io::message::error);
                return;
              }
              atom--;
              if (atom >= topo.num_atoms()) {
                std::ostringstream msg;
                msg << "Cartesian: atom indices (" << atom + 1
                    << ") out of range.";
                io::messages.add(msg.str(), "In_BSLEUS", io::message::error);
                return;
              }
              DEBUG(8, "Cartesian: Adding atom " << atom + 1);
              atom_list.push_back(atom);
              cartAtoms.push_back(atom);
            }
          }
          coord = new util::BS_Cartesian(id,atom_list, allAtoms, reference);
          break;
        }
        case util::BS_Coordinate::dihedralSum: {
          unsigned int i = 0, j = 0, k = 0, l = 0, ii = 0, jj = 0, kk = 0, ll = 0;
          _lineStream >> i >> j >> k >> l >> ii >> jj >> kk >> ll;
          if (_lineStream.fail()) {
            io::messages.add("bad line in BSLEUSCOORD block (Dihedral angle)",
                    "In_BSLEUS", io::message::error);
            return;
          }
          // Convert to gromos
          i--, j--, k--, l--, ii--, jj--, kk--, ll--;
          if (i > topo.num_atoms() ||
                  j > topo.num_atoms() ||
                  k > topo.num_atoms() ||
                  l > topo.num_atoms() || 
                  ii > topo.num_atoms() ||
                  jj > topo.num_atoms() ||
                  kk > topo.num_atoms() ||
                  ll > topo.num_atoms()) {
            std::ostringstream msg;
            msg << "Dihedral (" << i + 1 << "-" << j + 1 << "-" << k + 1 << "-"
                    << l + 1 << ") or (" 
                    << ii + 1 << "-" << jj + 1 << "-" << kk + 1 << "-"
                    << ll + 1 << ") atom indices out of range.";
            io::messages.add(msg.str(), "In_BSLEUS", io::message::error);
            return;
          }

          coord = new util::BS_DihedralSum(id, i, j, k, l, ii, jj, kk, ll, reference);
          break;
        }
        case util::BS_Coordinate::lambda: {
          coord = new util::BS_Lambda(id, reference);
          break;
        }
        default: {
          io::messages.add("Unknown Type in BSLEUSCOORD Block!",
                  "In_BSLEUS", io::message::error);
          return;
        }
      } // end switch
      bs_subspaces[subspace]->addCoordinate(coord);
      DEBUG(10, coord->str() << " added to subspace " << subspace + 1);
      numCoordRead++;
      numCoordsPerSubspace[subspace]++;
      references[subspace].push_back(reference);
    }

    if (numCoords != numCoordRead) {
      std::ostringstream os;
      os << "The numbers of coordinates in BSLEUSCOORD seems wrong! "
         << "Expected " << numCoords <<", but read " << numCoordRead << " in!";
      io::messages.add(os.str(), "In_BSLEUS", io::message::error);
      return;
    }
    DEBUG(5, "Read in " << numCoordRead << " different definitions of coordinates.");
  } // BSLEUSCOORD
  
  // ==============================================================  
  std::vector<std::vector<int> > num_dimensions; // [subspace][coordinate]
  for (int i = 0; i < num_subspaces; i++){
    num_dimensions.push_back(bs_subspaces[i]->getDimensionality());
    //DEBUG(5, "BS_Subspace(" << i+1 << ") has " << num_dimensions.back() << " dimensions");
  }
  
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
    
    int subspace = 0, numSpheres = 0, numSpheresRead = 0;
    _lineStream >> numSpheres;
    if (_lineStream.fail()) {
      io::messages.add("Couldn't get the number of Spheres in BSLEUSSPH",
              "In_BSLEUS", io::message::error);
      return;
    }
    
    int id = 0, num_gp = 0;
    double radius = 0.0, forceConst = 0.0;
    std::vector<double> centerValues;
    util::BS_Vector center;
    for (; it != to; it++){
      _lineStream.clear();
      _lineStream.str(*it);
      _lineStream >> id >> subspace;
      // Convert to GROMOS
      subspace--;
      
      // Get the center
      for (unsigned int i = 0; i < numCoordsPerSubspace[subspace]; i++) {
        for (int j = 0; j < num_dimensions[subspace][i];) {
          std::string coordStr;
          DEBUG(10, "Reading in Dimension no. " << j + 1 << " of " << num_dimensions[subspace][i]);
          _lineStream >> coordStr;
          if (_lineStream.fail()) {
            io::messages.add("bad center in BSLEUSSPH block",
                    "In_BSLEUS", io::message::error);
            return;
          }
          std::vector<double> coords;
          parseSpecifier(topo, sim, conf, 
                         coordStr, refFiles, cartAtoms, 
                         coords, os);
          for (unsigned int k = 0; k < coords.size(); k++, j++)
            centerValues.push_back(coords[k] / references[subspace][i]);
        }
      }
      
      // The rest of the constants
      _lineStream >> radius >> forceConst >> num_gp;
      if (_lineStream.fail()) {
        io::messages.add("bad line in BSLEUSSPH block",
                "In_BSLEUS", io::message::error);
        return;
      }
      if (num_gp < 1){
        io::messages.add("BSLEUSSPH: number of grid points must be at least 1!",
                "In_BSLEUS", io::message::error);
        return;
      }
      center.create(centerValues);
      DEBUG(10, "Center: " + center.str());
      util::BS_Sphere *bs_sphere = new util::BS_Sphere(id, num_gp, forceConst, 
                                                        center, radius);
      bs_subspaces[subspace]->addPotential(bs_sphere);
      numSpheresRead++;
      centerValues.clear();
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
    
    int subspace = 0, numSticks = 0, numSticksRead = 0;
    _lineStream >> numSticks;
    if (_lineStream.fail()) {
      io::messages.add("Couldn't get the number of Sticks in BSLEUSSTK",
              "In_BSLEUS", io::message::error);
      return;
    }
    
    int id = 0, num_gp = 0, defType = 0, startSphere = 0, endSphere = 0;
    double width = 0.0, forceConst = 0.0;
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
      
      // Start and End Points of the Sticks
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
        } // case 0: start & end spheres
        case 1: {// define end and start points with coordinates
          std::vector<double> coords;//, periodicities;
          for (int m = 0; m < 2; m++) { // loop over start and end coordinates
            coords.clear();
            for (unsigned int i = 0; i < numCoordsPerSubspace[subspace]; i++) {
              for (int j = 0; j < num_dimensions[subspace][i];) {
                std::string coordStr;
                _lineStream >> coordStr;
                if (_lineStream.fail()) {
                  io::messages.add("bad line in BSLEUSSTK block",
                          "In_BSLEUS", io::message::error);
                  return;
                }
                std::vector<double> new_coords;
                parseSpecifier(topo, sim, conf, 
                               coordStr, refFiles, cartAtoms, 
                               new_coords, os);
                for (unsigned int k = 0; k < new_coords.size(); k++, j++) {
                  DEBUG(10, "(" << i << ", " << j << ", " << k << ") Push in " << new_coords[k] << " / " << references[subspace][i])
                  coords.push_back(new_coords[k] / references[subspace][i]);
                }
              }
            }
            if (m == 0){
              start.create(coords);
            } else {
              end.create(coords);
            }
          }
          break;
        } // case 1: coordinates
        default: {
          io::messages.add("Unknown specifier for coordinate definition of sticks!",
                "In_BSLEUS", io::message::error);
          return;
        }
      }
      DEBUG(10, "Start and end");
      DEBUG(10, start.str());
      DEBUG(10,end.str());
      
      // Specifications of the stick
      _lineStream >> width >> forceConst >> num_gp;
      if (_lineStream.fail()) {
        io::messages.add("bad line in BSLEUSSTK block",
                "In_BSLEUS", io::message::error);
        return;
      }
      if (num_gp < 1){
        io::messages.add("BSLEUSSTK: number of grid points must be at least 1!",
                "In_BSLEUS", io::message::error);
        return;
      }
      
      util::BS_Stick *bs_stick = new util::BS_Stick(id, num_gp, forceConst,
                                                        start, end, width);
      bs_subspaces[subspace]->addPotential(bs_stick);
      numSticksRead++;
    }
    DEBUG(5, "The number of sticks according to file: " << numSticks << "; actually read: " << numSticksRead);
    if (numSticks != numSticksRead) {
      io::messages.add("The numbers of Sticks in BSLEUSSTK seems wrong!",
              "In_BSLEUS", io::message::warning);
      return;
    }
    
  } // BSLEUSSTK

  // ==============================================================
  // BSSNAKE
  buffer = m_block["BSLEUSSNAKE"];
  DEBUG(10, "BSSNAKE block : " << buffer.size());

  if (!buffer.size()) {
    io::messages.add("no BSLEUSSNAKE block in B&S-LEUS definition file",
            "In_BSLEUS", io::message::notice);
  } else {

    std::vector<std::string>::const_iterator it = buffer.begin() + 1,
            to = buffer.end() - 1;

    DEBUG(10, "reading in BSLEUSSNAKE data");
    _lineStream.clear();
    _lineStream.str(*it++);

    int subspace = 0, numSnakes = 0, numSnakesRead = 0;
    _lineStream >> numSnakes;
    if (_lineStream.fail()) {
      io::messages.add("Couldn't get the number of Snakes in BSLEUSSNAKE",
              "In_BSLEUS", io::message::error);
      return;
    }

    int id = 0, num_gp = 0, numPoints = 0;
    double half_width = 0.0, forceConst = 0.0;
    for (; it != to; it++) {
      _lineStream.clear();
      _lineStream.str(*it);
      _lineStream >> id >> subspace >> numPoints;
      if (_lineStream.fail()) {
        io::messages.add("bad line in BSLEUSSNAKE block",
                "In_BSLEUS", io::message::error);
        return;
      }
      // Convert to GROMOS
      subspace--;

      std::vector<double> coords;
      std::vector<util::BS_Vector> points;
      util::BS_Vector point;

      for (int m = 0; m < numPoints; m++) {
        coords.clear();
        for (unsigned int i = 0; i < numCoordsPerSubspace[subspace]; i++) {
          for (int j = 0; j < num_dimensions[subspace][i];) {
            std::string coordStr;
            _lineStream >> coordStr;
            if (_lineStream.fail()) {
              io::messages.add("bad line in BSLEUSSNAKE block",
                      "In_BSLEUS", io::message::error);
              return;
            }
            std::vector<double> new_coords;
            parseSpecifier(topo, sim, conf,
                    coordStr, refFiles, cartAtoms,
                    new_coords, os);
            for (unsigned int k = 0; k < new_coords.size(); k++, j++) {
              DEBUG(10, "(" << i << ", " << j << ", " << k << ") Push in " << new_coords[k] << " / " << references[subspace][i])
              coords.push_back(new_coords[k] / references[subspace][i]);
            }
          }
        }

        point.create(coords);
        points.push_back(point);
        DEBUG(10, "Point [" << m << "]: " << point.str());
      }

      // Specifications of the stick
      _lineStream >> half_width >> forceConst >> num_gp;
      if (_lineStream.fail()) {
        io::messages.add("bad line in BSLEUSSNAKE block",
                "In_BSLEUS", io::message::error);
        return;
      }
      if (num_gp < 0) {
        io::messages.add("BSLEUSSNAKE: number of grid points must be at least 0!",
                "In_BSLEUS", io::message::error);
        return;
      }

      util::BS_Snake *bs_snake = new util::BS_Snake(id, num_gp, forceConst,
              points, half_width);
      bs_subspaces[subspace]->addPotential(bs_snake);
      numSnakesRead++;
    }
    DEBUG(5, "The number of snake according to file: " << numSnakes << "; actually read: " << numSnakesRead);
    if (numSnakes != numSnakesRead) {
      io::messages.add("The numbers of Snakes in BSLEUSSNAKE seems wrong!",
              "In_BSLEUS", io::message::warning);
      return;
    }

  } // BSSNAKE
  // ==============================================================
  // BSLEUSPIPE
  buffer = m_block["BSLEUSPIPE"];
  DEBUG(10, "BSPIPE block : " << buffer.size());

  if (!buffer.size()) {
    io::messages.add("no BSLEUSPIPE block in B&S-LEUS definition file",
            "In_BSLEUS", io::message::notice);
  } else {

    std::vector<std::string>::const_iterator it = buffer.begin() + 1,
            to = buffer.end() - 1;

    DEBUG(10, "reading in BSLEUSPIPE data");
    _lineStream.clear();
    _lineStream.str(*it++);

    int subspace = 0, numPipes = 0, numPipesRead = 0;
    _lineStream >> numPipes;
    if (_lineStream.fail()) {
      io::messages.add("Couldn't get the number of Pipes in BSLEUSPIPE",
              "In_BSLEUS", io::message::error);
      return;
    }

    int id = 0, num_gp_l = 0, num_gp_p = 0, spec_type = 0;
    double forceConst = 0.0;
    for (; it != to; it++) {
      _lineStream.clear();
      _lineStream.str(*it);
      _lineStream >> id >> subspace >> spec_type;
      if (_lineStream.fail()) {
        io::messages.add("bad line in BSLEUSPIPE block",
                "In_BSLEUS", io::message::error);
        return;
      }
      // Convert to GROMOS
      subspace--;

      util::BS_Vector point;
      util::BS_Pipe_Param params[2];

      for (int m = 0; m < 2; m++) { // start, end point
        if (spec_type == 0) {
          int sph_id = 0;
          _lineStream >> sph_id;
          if (_lineStream.fail()) {
            io::messages.add("Could not read the sphere id in BSLEUSPIPE block",
                    "In_BSLEUS", io::message::error);
            return;
          }
          point = bs_subspaces[subspace]->getCenter(sph_id);
        } else {
          std::vector<double> coords;
          for (unsigned int i = 0; i < numCoordsPerSubspace[subspace]; i++) {
            for (int j = 0; j < num_dimensions[subspace][i];) {
              std::string coordStr;
              _lineStream >> coordStr;
              if (_lineStream.fail()) {
                io::messages.add("bad line in BSLEUSPIPE block",
                        "In_BSLEUS", io::message::error);
                return;
              }
              std::vector<double> new_coords;
              parseSpecifier(topo, sim, conf,
                      coordStr, refFiles, cartAtoms,
                      new_coords, os);
              for (unsigned int k = 0; k < new_coords.size(); k++, j++) {
                DEBUG(10, "(" << i << ", " << j << ", " << k << ") Push in " << new_coords[k] << " / " << references[subspace][i])
                coords.push_back(new_coords[k] / references[subspace][i]);
              }
            }
          }
          point.create(coords);
        }
        DEBUG(10, "Point [" << m << "]: " << point.str());
        params[m].point = point;

        _lineStream >> params[m].inner_width >> params[m].outer_width;
        if (_lineStream.fail()) {
          io::messages.add("Could not read the half widths in BSLEUSPIPE block",
                  "In_BSLEUS", io::message::error);
          return;
        }
      }

      _lineStream >> forceConst >> num_gp_l >> num_gp_p;
      if (_lineStream.fail()) {
        io::messages.add("bad line in BSLEUSPIPE block",
                "In_BSLEUS", io::message::error);
        return;
      }
      if (num_gp_l < 0 || num_gp_p < 0) {
        io::messages.add("BSLEUSPIPE: number of grid points must be at least 0!",
                "In_BSLEUS", io::message::error);
        return;
      }

      util::BS_Pipe *bs_pipe = new util::BS_Pipe(id, num_gp_l, num_gp_p,
              forceConst, params[0], params[1]);
      bs_subspaces[subspace]->addPotential(bs_pipe);
      numPipesRead++;
    }
    DEBUG(5, "The number of snake according to file: " << numPipes << "; actually read: " << numPipesRead);
    if (numPipes != numPipesRead) {
      io::messages.add("The numbers of Pipes in BSLEUSPIPE seems wrong!",
              "In_BSLEUS", io::message::warning);
      return;
    }
  } // BSPIPE
  
  conf.special().bs_umbrella.addSubspaces(bs_subspaces);
}

void
io::In_BSLEUS::parseSpecifier(topology::Topology &topo,
                              simulation::Simulation &sim,
                              configuration::Configuration& conf, 
                              std::string& coordStr, 
                              std::map<unsigned int,std::string>& refFiles, 
                              std::vector<unsigned int> &cartAtoms,
                              std::vector<double>& coords,
                              std::ostream & os)
{
  size_t found = 0;;
  if ((found = coordStr.find("{") )!= std::string::npos) {
    size_t to = coordStr.find("}");
    std::string refFileName(coordStr, found + 1, --to);
    io::igzstream refFile;
    refFile.open(refFileName.c_str());
    if (!refFile.is_open()) {
      io::messages.add("opening reference structure file failed!\n",
              "In_BSLEUS", io::message::error);
      return;
    }
    configuration::Configuration myConf;
    simulation::Simulation mySim = sim;
    mySim.param().boundary.boundary = sim.param().boundary.boundary;
    io::In_Configuration in_conf(refFile);
    in_conf.read_position_plain(topo, myConf, mySim, os);
    math::VArray &refpos = myConf.current().pos;
    
    put_into_box(myConf, refpos);
    
    if (cartAtoms.size() == 0) { // all atoms
      for (unsigned int i = 0; i < refpos.size(); i++){
        for (int j = 0; j < 3; j++){
          coords.push_back(refpos[i][j]);
        }
        DEBUG(8, "Using atom " << i << " with the position " << v2s(refpos[i]));
      }          
    }
    else {
      for (unsigned int i = 0; i < cartAtoms.size(); i++){
        for (int j = 0; j < 3; j++){
          coords.push_back(refpos[cartAtoms[i]][j]);
        }
        //DEBUG(8, "Get Position of atom " << cartAtoms[i] + 1)
        DEBUG(8, "Using atom " << i << " with the position " << v2s(refpos[i]));
      }
    }
  }
  else {
    coords.push_back(atof(coordStr.c_str()));
  }
}
void io::In_BSLEUS::
put_into_box(configuration::Configuration &conf, math::VArray& pos){
  SPLIT_BOUNDARY(_put_into_box, conf, pos);
}

template<math::boundary_enum B>
void io::In_BSLEUS::
_put_into_box(configuration::Configuration& conf, math::VArray& pos){
  DEBUG(4, "Put the coordinates into the GROMOS box");
  math::Periodicity<B> periodicity(conf.current().box);
  
  math::VArray::iterator it = pos.begin(),
          to = pos.end();
  for (; it != to; it++){
    periodicity.put_into_positive_box(*it);
    DEBUG(5, v2s(*it));
  }
}
