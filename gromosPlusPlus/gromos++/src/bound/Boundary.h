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

#ifndef INCLUDED_BOUND_BOUNDARY
#define INCLUDED_BOUND_BOUNDARY

#include <string>
#include <vector>

namespace gmath{
  class Vec;
}

namespace gcore{
  class System;
  class Box;
}

using gmath::Vec;

namespace bound{

  class Boundary_i;
  /**
   * Class Boundary
   * defines some basic function for the periodic boundary 
   * conditions. It provides a virtual function for the nearest image
   * calculation and virtual functions to the specific @ref gathermethods "gather methods".
   * Usually it is constructed by the @ref args::BoundaryParser class.
   *
   * @section gathermethods Gather Methods
   *
   * The following methods are recognized by the @ref args::GatherParser "GatherParser",
   * usually given as the second argument to \@pbc after the boundary condition:
   * 
   * - @ref bound::Boundary::nogather() "nog or 0"
   *    do not gather 
   * - @ref bound::Boundary::gatherlist() "glist or 1"
   *   gathering based on a list of atom pairs
   *   - The atom pair should be in the sequence: A B, where A is an atom of 
   *   the molecule to be gathered, and B is an atom of the reference molecule.
   * - @ref bound::Boundary::gathertime()  "gtime or 2"
   *   gathering based on previous frame
   * - @ref bound::Boundary::gatherref()  "gref or 3"
   *   gathering based on a reference structure
   * - @ref bound::Boundary::gatherltime() "gltime or 4"
   *   gather first frame based on a list, next frames based on previous frame
   * - @ref bound::Boundary::gatherrtime()  "grtime or 5"
   *   gather first frame based on a reference structure, next frames based on previous frame
   * - @ref bound::Boundary::gatherbond()  "gbond or 6"
   *   gather based on bond connectivity
   * - @ref bound::Boundary::coggather()  "cog or 7"
   *   gather with respect to the centre of geometry of all atoms of the first molecule in the system
   * - @ref bound::Boundary::gfitgather()   "gfit or 8"
   *   gather based on a reference structure which has been superimposed on the 
   *   first frame of the trajectory
   *   - This method should work with multimers for which otherwise often only the
   *     glist method works, but it has the advantage that no manual setting up 
   *     of an atom pair list is necessary. If a molecule selection is given with
   *     secondary argument "molecules", the remaining molecules will not be 
   *     gathered to their reference positions but to the overall cog of the 
   *     selected molecules, which is usually desired for ions or lipids.
   *   - The method depends on an existing correctly-gathered reference frame.
   *
   * Some of the methods take compulsory or optional secondary arguments:
   * - refg &lt;reference coordinates&gt; (gref, grtime, gtime, gfit)
   * - list &lt;atom pair list&gt; (glist, gltime, gtime)
   * - molecules &lt;molecule numbers&gt; (gfit)
   *
   * Usage examples:
   * - @verbatim @pbc r gref refg coord.cnf @endverbatim
   * - @verbatim @pbc r gltime list 2:res(15:CA) 1:res(43:CA) 3:134 1:res(43:CA) @endverbatim
   * - @verbatim @pbc r gfit refg coord.cnf molecules 1-5 @endverbatim
   *
   * @class Boundary
   * @author R. Buergi, M.K. Kastenholz
   * @ingroup bound
   */
  class Boundary {

    Boundary_i *d_this;
    std::vector<int > d_refmol;
    bool firstframe;

    // not implemented
    Boundary& operator=(const Boundary& rhs);
    Boundary(const Boundary& b);
    Boundary();

  public:
    /** Constructor */
    Boundary (gcore::System *sys);
    /** Destructor */
    virtual ~Boundary();

    // Methods

    /**
     * sets reference position for molecule i to v.
     */
    void setReference(int i, const gmath::Vec &v);
    /**
     * sets reference position for all molecules to the first
     * atom of the corresponding molecule in sys
     */
    void setReference(gcore::System const & sys);
   
    /**
     * Given the reference position r1, we give r2 back so that 
     * r2 is the nearest image to r1. Used to reconnect molecules.
     * Note that solvent molecules do never need to be reconnected
     */
    virtual gmath::Vec nearestImage(const gmath::Vec &r1,
				    const  gmath::Vec &r2, 
				    const gcore::Box &box) const = 0;
   
    /**
     * Using nearestImage, check if the v is in the same box as the ref
     */
    virtual bool inBox(const gmath::Vec &ref,
				       const gmath::Vec &v,
				       const gcore::Box &box) const;
    
    /**
     * No gathering
     */
    virtual void nogather(); 
    /**
     * gather based on a list of atom pairs. The atom pair should be in the
     * sequence: A B, where A is an atom of the molecule to be gathered, and B
     * is an atom of the reference molecule.
     */
    virtual void gatherlist();
    /**
     * gather in term of time
     */
    virtual void gathertime();
    /**
     * gather based on a reference structure
     */
    virtual void gatherref();
    /**
     * gather the first frame based on an atom list, then the rest in term of time
     */
    virtual void gatherltime();
    /**
     * gather the first frame based on a reference, then the rest in term of time
     */
    virtual void gatherrtime();
    /**
     * gather based on bond connection
     */
    virtual void gatherbond();
    /**
     * gather solute and solvent with respect to the cog of mol(0)
     */
    virtual void coggather();
    /**
     * Method to work with multimers or protein/membrane systems.
     * Selected molecules: gather first frame w.r.t. to a reference which was 
     * superimposed on the first molecule of the first frame, gather following 
     * frames to previous frame.
     * Always gather non-selected molecules (ions, lipids) with respect to cog 
     * of selected molecules.
     */
    virtual void gfitgather();
    /**
     * reference vector (set to pos(0) of mol(i)) of each molecule upon 
     * creation of object boundary.
     * if the system does not have any coordinates yet, they are initialized
     * with 0.0.
     */
    const gmath::Vec &reference(int i)const;
    /**
     * system accessor.
     */
    gcore::System &sys();
    /**
     * reference system accessor.
     */
    gcore::System &refSys();
    /**
     * the boundary type.
     */
    char type() const;
    /**
     * set the boundary type.
     */
    void setType(char t);
    /**
     * set to the path of a reference frame that is used in the
     * reference gathering method
     */
    void setReferenceFrame(std::string file);
    /**
     * set the reference system to sys
     */
    void setReferenceSystem(gcore::System system);
    /**
     * add molecules to be gathered by reference
     */
    void addRefMol(int molnum);
    /**
     * member pointer to gather function
     */
    typedef void (Boundary::*MemPtr)();
    
    // old/deprecated gathering methods
    /**
     * @deprecated old gathering methods which should not be used anymore
     * gather the whole System in gromos style (per molecule). 
     * Solvent is not regarded.
     */
    virtual void gathergr();
    /**
     * @deprecated old gathering methods which should not be used anymore
     * gather each molecule in gromos style (per molecule), then translate
     * the molecules so that their cogs are inside the box. 
     * Solvent is not regarded.
     */
    virtual void gathermgr();
    /**
     * @deprecated old gathering methods which should not be used anymore
     * gather the whole system in gromos++ style (per first molecule).
     */
    virtual void gather();
    /**
     * @deprecated old gathering methods which should not be used anymore.
     * gathering of e.g amyloid crystals.
     * 1) Gather each molecule by sequence to its first atom.
     * 2) Gather cogs w.r.t. cog of the previous molecule.
     * 3) Gather molecules w.r.t. to the new cog positions.
     * 4) Gather solvent to overall cog.
     */
    virtual void crsgather(); 
    /**
     * @deprecated old gathering methods which should not be used anymore
     * same as crsgather, but in 2) gather to overall cog of previous molecules
     */
    virtual void seqgather();
    /**
     * @deprecated old gathering methods which should not be used anymore
     * attempt for a generalized gathering method (A. Choutko / D. Geerke / A.-P. Kunz).
     * as crsgather, but in 2) gather cogs in order of closeness, the closest one first
     */ 
    virtual void gengather();
  };

}
#endif


