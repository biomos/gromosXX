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

#ifndef INCLUDED_GCORE_PTTOPOLOGY
#define INCLUDED_GCORE_PTTOPOLOGY

#include <set>
#include <string>
#include <vector>

#include "CrossDihedral.h"
#include "AtomPair.h"
#include "Angle.h"
#include "Bond.h"
#include "Angle.h"
#include "Dihedral.h"
#include "Improper.h"
#include "System.h"

namespace gcore
{
  /**
   * @class AtomPairParam
   * holds an AtomPair and a parameter
   */
  class AtomPairParam : public AtomPair {
    int d_param;
  public:
    /**
     * constructor
     * @param a the first atom
     * @param b the second atom
     * @param p the connected parameter
     */
    AtomPairParam(int a, int b, int p) : gcore::AtomPair(a, b), d_param(p) {}
    /**
     * accessor to the parameter
     */
    int & param() { return d_param; }
    /**
     * const accessor to the parameter
     */
    int param() const { return d_param; }
  };
  
  /**
   * Class PtTopology
   * contains one or several perturbation topologies for 
   * nonbonded interactions
   *
   * @class PtTopology
   * @author C. Oostenbrink, N. Schmid
   * @ingroup gcore
   */
  class PtTopology
  {
    std::vector<int> d_atomnum;
    std::vector<int> d_residuenum;
    std::vector<std::string> d_atomnames;
    std::vector<std::string> d_residuenames;
    std::vector<std::string> d_pertnames;
    std::vector< std::vector <int> > d_iac;
    std::vector< std::vector <double> > d_charge;
    std::vector< std::vector <double> > d_mass;
    std::vector< std::vector <double> > d_polarisability;
    std::vector< std::vector <double> > d_dampingLevel;
    std::vector<double> d_alphaLJ;
    std::vector<double> d_alphaCRF;
    bool d_hasPolarisationParams;
    bool d_multipt;
    std::vector<std::set<AtomPairParam> > d_atompairs;
    std::vector<std::set<Bond> > d_bonds;
    std::vector<std::set<Angle> > d_angles;
    std::vector<std::set<Improper> > d_impropers;
    std::vector<std::set<Dihedral> > d_dihedrals;
    std::vector<std::set<CrossDihedral> > d_crossdihedrals;
  
  public:
    /**
     * constructor
     */
    PtTopology() : d_hasPolarisationParams(false), d_multipt(false) {};
    /**
     * copy constructor
     */
    PtTopology(const PtTopology & pt);
    /**
     * construct the perturbation topology from
     * two states.
     * 
     * @param sysA the first system
     * @param sysB the second system
     * @param quiet warnings omitted if true
     */
    PtTopology(System & sysA, System & sysB, bool quiet = false);
    /**
     * construct the multiple perturbation topology from
     * many systems (bonded and exclusions ignored)
     *
     * @param sys a vector containing pointers to the states.
     * @param quiet warnings omitted if true
     */
    PtTopology(std::vector<System*> & sys, bool quiet = false);
    /**
     * deconstructor
     */
    ~PtTopology(){};
    /**
     * set the dimensions of the perturbation topology
     * @param a number of atoms that are to be perturbed
     * @param p number of perturbation topologies that 
     *          should be stored
     */
    void setSize(int a, int p);
    /**
     * function to set the iac for atom a in perturbation p
     */
    void setIac(int a, int p, int iac);
    /**
     * function to set the mass for atom a in perturbation p
     */
    void setMass(int a, int p, double m);
    /**
     * function to set the charge for atom a in perturbation p
     */
    void setCharge(int a, int p, double q);
    /**
     * function to set the polarisability for atom a in perturbation p
     */
    void setPolarisability(int a, int p, double al);
    /**
     * function to set the damping level for atom a in perturbation p
     */
    void setDampingLevel(int a, int p, double l);
    /**
     * function to set the atom name of atom a
     */
    void setAtomName(int a, std::string name);
    /**
     * function to set the residue name of atom a
     */
    void setResidueName(int a, std::string name);
    /**
     * function to set the name of perturbation p
     */
    void setPertName(int p, std::string name);
    /**
     * function to set the atom number of atom a
     */
    void setAtomNum(int a, int num);
    /**
     * function to set the residue number of atom a
     */
    void setResidueNum(int a, int num);
    /**
     * function to set the softness for LJ for of atom a
     */
    void setAlphaLJ(int a, double al);
    /**
     * function to set the softness for CRF for of atom a
     */
    void setAlphaCRF(int a, double al);
    /**
     * function to set whether the topology contained polarisation
     * parameters
     */
    void setHasPolarisationParameters(bool pol);
    /**
     * function to set whether the pt topology is in a format
     * which can contain more than two perturbations (then
     * we have only IAC and charges)
     */
    void setMultiPt(bool multipt);
    /**
     * function to apply a given perturbation to the system
     *
     * The strategy is the following: The topology is converted to a linear 
     * version and the changes are applied to it. The linear version is then
     * used to create a new topology. Finally the topology of the molecules
     * is exchanged by the new version. Positions, velocities etc. are not 
     * affected.
     *
     * @param sys the system to which the perturbation is applied
     * @param iipt the perturbation which is applied (default 1, state B)
     * @param first The perturbation topology is shifted by this value
     */
    void apply(gcore::System &sys, int iipt=1, int first = 0)const;
    /**
     * accessor to a vector of the charges in perturbation p
     */
    std::vector<double> charge(int p=0)const{return d_charge[p];}
    /**
     * accessor to the charge of atom a in perturbation p
     */
    double charge(int a, int p)const{return d_charge[p][a];}
    /**
     * accessor to a vector of the iac in perturbation p
     */
    std::vector<int> iac(int p=0)const{return d_iac[p];}
    /**
     * accessor to the iac of atom a in perturbation p
     */
    int iac(int a, int p)const{return d_iac[p][a];}
    /**
     * accessor to a vector of the masses in perturbation p
     */
    std::vector<double> mass(int p=0)const{return d_charge[p];}
    /**
     * accessor to the mass of atom a in perturbation p
     */
    double mass(int a, int p)const{return d_mass[p][a];}
    /**
     * accessor to a vector of the polaisabilities in perturbation p
     */
    std::vector<double> polarisability(int p=0)const{return d_polarisability[p];}
    /**
     * accessor to the polarisability of atom a in perturbation p
     */
    double polarisability(int a, int p)const{return d_polarisability[p][a];}
    /**
     * accessor to a vector of the damping levels in perturbation p
     */
    std::vector<double> dampingLevel(int p=0)const{return d_dampingLevel[p];}
    /**
     * accessor to the damping level of atom a in perturbation p
     */
    double dampingLevel(int a, int p)const{return d_dampingLevel[p][a];}
    /**
     * accessor to the the atom pairs in perturbation p
     */
    std::set<gcore::AtomPairParam> & atompairs(int p){return d_atompairs[p];}
    /**
     * accessor to the the atom pairs in perturbation p (const)
     */
    const std::set<gcore::AtomPairParam> & atompairs(int p)const{return d_atompairs[p];}
    /**
     * accessor to the the bonds in perturbation p
     */
    std::set<gcore::Bond> & bonds(int p){return d_bonds[p];}
    /**
     * accessor to the the bonds in perturbation p (const)
     */
    const std::set<gcore::Bond> & bonds(int p)const{return d_bonds[p];}
    /**
     * accessor to the the angles in perturbation p
     */
    std::set<gcore::Angle> & angles(int p){return d_angles[p];}
    /**
     * accessor to the the angles in perturbation p
     */
    const std::set<gcore::Angle> & angles(int p)const{return d_angles[p];}
    /**
     * accessor to the the improper dihedrals in perturbation p
     */
    std::set<gcore::Improper> & impropers(int p){return d_impropers[p];}
    /**
     * accessor to the the improper dihedrals in perturbation p (const)
     */
    const std::set<gcore::Improper> & impropers(int p)const{return d_impropers[p];}
    /**
     * accessor to the the dihedrals in perturbation p
     */
    std::set<gcore::Dihedral> & dihedrals(int p){return d_dihedrals[p];}
    /**
     * accessor to the the dihedrals in perturbation p (const)
     */
    const std::set<gcore::Dihedral> & dihedrals(int p)const{return d_dihedrals[p];}
    /**
     * accessor to the the cross dihedrals in perturbation p
     */
    std::set<gcore::CrossDihedral> & crossdihedrals(int p){return d_crossdihedrals[p];}
    /**
     * accessor to the the crossdihedrals in perturbation p (const)
     */
    const std::set<gcore::CrossDihedral> & crossdihedrals(int p)const{return d_crossdihedrals[p];}
    
    /**
     * accessor to the name of atom a
     */
    std::string atomName(int a)const{return d_atomnames[a];}
    /**
     * accessor to the reside name of atom a
     */
    std::string residueName(int a)const{return d_residuenames[a];}
    /**
     * accessor to the name of perturbation p
     */
    std::string pertName(int p)const{return d_pertnames[p];}
    /**
     * accessor to the atom number of atom a
     */
    int atomNum(int a)const{return d_atomnum[a];}
    /**
     * accessor to the residue number of atom a
     */
    int residueNum(int a)const{return d_residuenum[a];}
    /**
     * accessor to the softness in LJ of atom a
     */
    double alphaLJ(int a)const{return d_alphaLJ[a];}
    /**
     * accessor to the softness in CRF of atom a
     */
    double alphaCRF(int a)const{return d_alphaCRF[a];}
    /**
     * accessor to the number of perturbations in the topology
     */
    int numPt()const{return d_iac.size();}
    /**
     * accessor to the number of atoms in the perturbation topology
     */
    int numAtoms()const{return d_atomnames.size();}
    /**
     * accessor to polarisation
     */
    bool hasPolarisationParameters()const{return d_hasPolarisationParams;}
    /**
     * accessor to whether we have a multi-perturbation topology
     * (which has only IAC and charges)
     */
    bool multiPt()const{return d_multipt;}
  
  };
}

#endif
