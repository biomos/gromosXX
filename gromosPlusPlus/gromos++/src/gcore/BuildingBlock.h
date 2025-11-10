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

// gcore_BuildingBlock.h

#ifndef INCLUDED_GCORE_BUILDINGBLOCK
#define INCLUDED_GCORE_BUILDINGBLOCK

#include <cassert>
#include <vector>
#include <string>
#include <set>
#include <sstream>

namespace gcore{

  class BbSolute;
  class SolventTopology;
  /**
   * Class BuildingBlock
   * Purpose: contains all different kinds of building blocks
   *
   * Description:
   * The BuildingBlock class contains all building blocks and additional
   * information from the gromos96 mtb-file. Three different kinds of
   * building blocks are known: MTBUILDBLSOLUTE (class BbSolute), 
   * MTBUILDBLSOLVENT (class SolventTopology) and MTBUILDBLEND (class 
   * BbSolute)
   *
   * @class BuildingBlock
   * @author C. Oostenbrink
   * @ingroup gcore
   * @version $Date: Wed Jul 10 14:00:00 MEST 2002
   * @sa gcore::BbSolute
   * @sa gcore::SolventTopology
   * @sa gcore::BBEnd
   */
  class BuildingBlock{
    std::vector<BbSolute*> d_bb;
    std::vector<BbSolute*> d_be;
    std::vector<SolventTopology*> d_bs;
    double d_fpepsi;
    double d_hbar;
    double d_spdl;
    double d_boltz;
    bool d_physConstRead;
    std::set<std::string> d_ffcode;
    int d_linkExclusions;
    bool d_empty;

  public:
    //Constructors
    /**
     * BuildingBlock constructor
     */
    BuildingBlock();
    /**
     * BuildingBlock copy constructor
     * @param bld BuildingBlock to be copied
     */
    BuildingBlock(const BuildingBlock &bld);
    /**
     * BuildingBlock deconstructor
     */
    ~BuildingBlock();

    // Methods
    /**
     * Member function addBuildingBlock adds another set of building blocks
     * to this one
     * They should have the same force field code, fpepsi, hbar and kB
     */
    void addBuildingBlock(const BuildingBlock &bld);
    /**
     * Member function operator = copies one set of building blocks
     * into the other
     */
    BuildingBlock &operator=(const BuildingBlock &bld);
    /**
     * Method to add a solute building block to the building blocks
     * @param mol a BbSolute (corresponds to a MTBUILDBLSOLUTE block)
     */
    void addBbSolute(const BbSolute &mol);
    /**
     * Method to add a solvent building block to the building blocks
     * @param sol a SolventTopology (corresponds to a MTBUILDBLSOLVENT 
     *        block
     */
    void addBbSolvent(const SolventTopology &sol);
    /**
     * Method to add an end group building block to the building blocks
     * @param mol a BbSolute that is an end group
     *            (corresponds to a MTBUILDBLEND block)
     */
    void addBbEnd(const BbSolute &mol);
    /**
     * Set the value of Fpepsi. It would probably make more sense to store
     * Fpepsi in the GromosForceField, but this is the gromos96 structure
     */
    void setFpepsi(const double a){d_fpepsi=a;};
    /**
     * Set the value of Hbar. It would probably make more sense to store 
     * Hbar in the GromosForceField, but this is the gromos96 structure
     */
    void setHbar(const double a){d_hbar=a;};
    /** Set the value of the speed of light. It would probably make more sense to store
     * it in the GromosForceField, but this is the gromos96 structure
     */
    void setSpdl(const double a){d_spdl=a;};
    /**
     * Set the value of kB. It would probably make more sense to store
     * kB in the GromosForceField, but this is the gromos96 structure
     */
    void setBoltz(const double a){d_boltz=a;};
    /**
     * To know if the physical constants block was read
     */
    void setPhysConstRead(bool b){d_physConstRead=b;};
    /**
     * Set the number of exclusions when linking (= number of trailing
     * atoms)
     */
    void setLinkExclusions(const int i){d_linkExclusions=i;};
    /**
     * Set the force field code
     */
    void setForceField(const std::string s){
      std::istringstream is(s);
      std::string temp_ffcode;
      while(is >> temp_ffcode){
        d_ffcode.insert(temp_ffcode);
      }
    };
    
    
    
    // Accessors
    /**
     * Accessor, returns the i-th BbSolute as a const
     */
    const BbSolute &bb(int i)const;
    /**
     * Accessor, returns the i-th BbSolute
     */
    BbSolute &bb(int i);
    /**
     * Accessor, returns the i-th end-group BbSolute as a const
     */
    const BbSolute &be(int i)const;
    /**
     * Accessor, returns the i-th BbEnd
     */
    BbSolute &be(int i);
    /**
     * Accessor, returns the i-th SolventTopology as a const
     */
    const SolventTopology &bs(int i)const;
    /**
     * Accessor, returns the i-th SolventTopology
     */
    SolventTopology &bs(int i);

    /**
     * Accessor, returns the number of BbSolutes in the BuildingBlock
     */
    int numBbSolutes()const;
    /**
     * Accessor, returns the number of BbSolvents (SolventTopologies) in 
     * the BuildingBlock
     */
    int numBbSolvents()const;
    /**
     * Accessor, returns the number of BbEnds in the BuildingBlock
     */
    int numBbEnds()const;
    /**
     * Accessor, returns the value of Fpepsi
     */
    double Fpepsi()const;
    /**
     * Accessor, returns the value of Hbar
     */
    double Hbar()const;
    /**
     *Accessor, returns the value of the speed of light
     */
    double Spdl()const;
    /**
     * Accessor, returns the value of kB
     */
    double Boltz()const;
    /**
     * Accessor, returns true if the physical constants have been read.
     */
    bool physConstRead()const;
    /**
     * Accessor, returns the number of exclusions for linking
     */
    int LinkExclusions()const;
    /**
     * Accessor, returns the set of force field codes
     */
    std::set<std::string> ForceField()const;
    /**
     * Method, returns an index for the first building block that is 
     * has the name s. This method searches through both the solute
     * and the end-group building blocks
     * @param s String with the building block name to search for
     * @return integer i with value<br>
     *         0  if s is not found<br>
     *         >0 s is found as the (i-1)-th solute building block.
     *            <i>Don't forget to substract the 1!</i><br>
     *         <0 s is found as the (abs(i)-1)-th end-group building 
     *            block. <i>Don't forget to convert i to the index!</i>
     */ 
    int findBb(std::string s);
    /**
     * Method, returns an index for the first building block that 
     * has the name s. This methos searches through both the solute
     * and the end-group building blocks and also gives back the total
     * of building blocks with name s that were found
     * @param s String with the building block name to search for
     * @param n integer that will be returned with the number of blocks
     * @return integer i with value<br>
     *         0  if s is not found<br>
     *         >0 s is found as the (i-1)-th solute building block.
     *            <i>Don't forget to substract the 1!</i><br>
     *         <0 s is found as the (abs(i)-1)-th end-group building 
     *            block. <i>Don't forget to convert i to the index!</i>
     */
    int findBb(std::string s, int &n);
    /** 
     * Method, returns an index for the first solvent building block that
     * has the name s.
     * @param s String with the solvent building block name to search for
     * @return integer i with value<br>
     *         0  if s is not found<br>
     *         >0 s is found as the (i-1)-th solvent building block.
     *            <i>Don't forget to substract the 1!</i>
     */
    int findBs(std::string s);
    /** 
     * Method, returns an index for the first solvent building block that
     * has the name s and also gives back the total number of building blocks
     * s that were found.
     * @param s String with the solvent building block name to search for
     * @return integer i with value<br>
     *         0  if s is not found<br>
     *         >0 s is found as the (i-1)-th solvent building block.
     *            <i>Don't forget to substract the 1!</i>
     */
    int findBs(std::string s, int &n);
    
    
};

  inline const BbSolute &BuildingBlock::bb(int i)const{
    assert (i < this->numBbSolutes());
    return *d_bb[i];
  }

  inline BbSolute &BuildingBlock::bb(int i){
    assert (i < this->numBbSolutes());
    return *d_bb[i];
  }

  inline const BbSolute &BuildingBlock::be(int i)const{
      assert (i < this->numBbEnds());
      return *d_be[i];
  }

  inline BbSolute &BuildingBlock::be(int i){
      assert (i < this->numBbEnds());
      return *d_be[i];
  }

  inline const SolventTopology &BuildingBlock::bs(int i)const{
      assert (i< this->numBbSolvents());
      return *d_bs[i];
  }

  inline SolventTopology &BuildingBlock::bs(int i){
      assert (i < this->numBbSolvents());
      return *d_bs[i];
  }

  inline int BuildingBlock::numBbSolutes()const{
    return d_bb.size();
  }

  inline int BuildingBlock::numBbEnds()const{
      return d_be.size();
  }

  inline int BuildingBlock::numBbSolvents()const{
    return d_bs.size();
  }

  inline double BuildingBlock::Fpepsi()const{
    return d_fpepsi;
  }
  
  inline double BuildingBlock::Hbar()const{
    return d_hbar;
  }

  inline double BuildingBlock::Spdl() const{
    return d_spdl;
  }

  inline double BuildingBlock::Boltz()const{
    return d_boltz;
  }
  
  inline bool BuildingBlock::physConstRead()const{
    return d_physConstRead;
  }
  
  inline int BuildingBlock::LinkExclusions()const{
    return d_linkExclusions;
  }

  inline std::set<std::string> BuildingBlock::ForceField()const{
    return d_ffcode;
  }
  
  
}
#endif

