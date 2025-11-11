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

// TrajArray.cc

#include "TrajArray.h"

#include "../gromos/Exception.h"
#include "../gcore/Molecule.h"
#include "../gcore/System.h"


using gcore::System;

// Constructor
namespace utils{

TrajArray::TrajArray(const System &sys) {

  nAtoms = 0;

  // count the atoms in the system to find frame size
  for (int molIndex = 0; molIndex < sys.numMolecules(); molIndex++) {
    nAtoms += sys.mol(molIndex).numAtoms();
  } 
}

// Destructor
TrajArray::~TrajArray() {
  for(unsigned int i = 0; i < trajectoryData.size(); i++)
    delete [] trajectoryData[i];
}

// Function to store a frame
void TrajArray::store(const gcore::System &sys, 
  const unsigned int frameIndex){

  int i;
  int nAtomsMol, molAtomIndex;
  int molIndex;
  unsigned int nAtomsSystem = 0;
  double *framePointer;

  // count the atoms in the system to find frame size
  for (molIndex = 0; molIndex < sys.numMolecules(); molIndex++)
    nAtomsSystem += sys.mol(molIndex).numAtoms();

  if (frameIndex < 0 || nAtomsSystem != nAtoms )
    throw gromos::Exception("TrajArray", "Unable to store frame.\n");

  if(frameIndex + 1 >= trajectoryData.size()){
  // guard againts off-by-one
    trajectoryData.resize(frameIndex + 1);
  }
  else
    delete [] trajectoryData[frameIndex];

  trajectoryData[frameIndex] = new double[3 * nAtoms];
  framePointer = trajectoryData[frameIndex];

  // read all the coords from the sys and store them in the array
  for (molIndex = 0; molIndex < sys.numMolecules(); molIndex++) {
    nAtomsMol = sys.mol(molIndex).numAtoms();
    for (molAtomIndex = 0; molAtomIndex < nAtomsMol; molAtomIndex++) {
      for(i = 0; i < 3; i++){
        *framePointer = sys.mol(molIndex).pos(molAtomIndex)[i];
        framePointer++;
      }
    }
  }
}

// Function to extract a frame
void TrajArray::extract( gcore::System &sys,
  const unsigned int frameIndex ) const {

  int i, molIndex, nAtomsMol, molAtomIndex;
  unsigned int nAtomsSystem = 0;
  double *framePointer;

  // count the atoms in the system to find frame size
  for (molIndex = 0; molIndex < sys.numMolecules(); molIndex++)
    nAtomsSystem += sys.mol(molIndex).numAtoms();
   
  if (frameIndex < 0 || frameIndex >= trajectoryData.size())
    throw gromos::Exception("TrajArray", 
      "Can't extract. Invalid frame index.\n");
  if(nAtomsSystem != nAtoms)
    throw gromos::Exception("TrajArray", 
      "Can't extract. Invalid system.\n");
  if(!(framePointer = trajectoryData[frameIndex]))
    throw gromos::Exception("TrajArray", 
      "Can't extract. Empty frame.\n");

  for (molIndex = 0; molIndex < sys.numMolecules(); molIndex++) {
    nAtomsMol = sys.mol(molIndex).numAtoms();
    for (molAtomIndex = 0; molAtomIndex < nAtomsMol; molAtomIndex++) {
      for(i = 0; i < 3; i++){
        sys.mol(molIndex).pos(molAtomIndex)[i] = *framePointer;
        framePointer++;
      }
    }
  }
}

unsigned int TrajArray::numAtoms() {
    return nAtoms;
}
}
