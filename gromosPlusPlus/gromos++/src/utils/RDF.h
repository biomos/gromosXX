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

#ifndef INCLUDED_RDF
#define INCLUDED_RDF

#include "Noe.h"
#include "../args/Arguments.h"

namespace utils {

  // the implementation class, just to let the compiler know that it exists
  class iRDF;

   /**
   * Class RDF:
   * a class to calculate radial distribution functions of a given system
   * described by a topology file and the corresponding trajectory files.
   *
   * The radial distribution function, g(r), is defined as
   *
   * @f[ g(r) = \frac{N_J(k)}{4\pi r^2 dr \rho_J} @f]
   *
   * which is the probability of finding a particle of type J at distance r
   * from a central particle of type I, relative to the probability for a homogenous
   * distribution of particle of type J around particle of type I.
   *
   * @class RDF
   * @ingroup utils
   * @author A. Eichenberger
   * */
  class RDF {

  private:
    
    // THE CLASSES DATA (or a pointer to the implementation class)
    // ================
    //
    /**
     * pointer to the data of the PDB class defined by the implementation
     * class iPDB
     */
    iRDF *d_this;

  public:

    // CONSTRUCTORS
    // ============
    //
    /**
     * Standard constructor, don't forget to set the system and trajectories
     * by hand!
     */
    RDF();
    /**
     * Constructor to initialize the class.
     * @param sys The system
     */
    RDF(gcore::System *sys, const args::Arguments *args);
    /**
     * Constructor
     */
    RDF(const RDF &rdf);


    // DESTRUCTORS
    // ===========
    /**
     * Destructor
     */
    ~RDF(void);


    // MEMBER FUNCTIONS
    // ================
    //
    /**
     * Sets the atoms to be considered as centre atoms
     */
     int addCenters(std::string s);
     /**
      * Adds an atom to the centre atoms
      */
     void addCentersAtom(int m, int a);
     /**
     * Removes all centre atoms
     */
     void clearCenters(void);
     /**
     * Sets the atoms to be considered as with atoms
     */
     int addWiths(std::string s);
     /**
      * Adds an atom to the with atoms
      */
     void addWithAtom(int m, int a);
    /**
     * Removes all with atoms
     */
     void clearWiths(void);
     /**
      * Calculate the rdf (all, intra- and intermolecular)
      */
     void calculateAll(void);
     /**
      * Calculate the rdf (intermolecular only)
      */
     void calculateInter(void);
     /**
      * Calculate the local number of particle j around particle i
      */
     void calculateLocal(void);
     /**
      * Calculate the rdf (intermolecular only) multiplied by the particle
      * density (needed for the calculation of neutron scattering intensities)
      *
      * The anom density is defined as numatoms / box_volume.
      *
      * @param numatoms the total number of atoms
      */
     double calculateInterPDens(unsigned int numatoms);
     /**
      * Sets all data points of the d_rdf vector to zero
      */
     void clearRDF(void);
     /**
      * Sets all data points of the d_local vector to zero
      */
     void clearLocal(void);
     /**
      * Sets the grid number for the rdf calculation to the number grid
      */
     void setGrid(unsigned int grid);
     /**
      * Sets the cutoff for the rdf calculation to cut
      */
     void setCut(double cut);
     /**
      * Sets a system to the class
      */
     void setSystem(gcore::System *sys);
     /**
      * Set the flag to compute dipole moment correlations between molecules
      */
     void setDCF(bool dcf);
     /**
      * Normalize the dipole moment correlations between molecules ?
      */
     void DCFnorm(bool dcfnorm);
     /**
      * Prints the contents of the d_rdf vector
      */
     void print(std::ostream &os);
     /**
      * Prints the contents of the d_local_mix and d_local_self vector
      */
     void printLocal(std::ostream &os);
     /**
      * Prints the contents of the d_dcf vector
      */
     void print_DCF(std::ostream &os);
     /**
      * Returns the rdf at position r
      */
     double rdf(unsigned int i);
     /**
      * Returns the dcf at position r
      */
     double dcf(unsigned int i);

  };

} /* end of namespace utils */
#endif
