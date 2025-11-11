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

// pb_FFTVacuumField.cc
#include "FFTVacuumField.h"

#include <cstdlib>
#include <cassert>

#include "FFTBoundaryCondition.h"
#include "FFTGridType.h"
#include "PB_Parameters.h"
#include "../fit/PositionUtils.h"
#include "../utils/AtomSpecifier.h"


using pb::FFTVacuumField;

FFTVacuumField::FFTVacuumField(utils::AtomSpecifier atoms, FFTGridType gt, FFTBoundaryCondition bci, ofstream &os):bc(os), ppp(os), gt(os){

  //FFTBoundaryCondition bc(os);
  //PB_Parameters ppp(os);

  this->atoms = atoms;
  this->gt = gt;
  this->bc = bc;
  bc = bci;
  this->tinynum=ppp.getTiny_real();
  this->csfunc=ppp.get_default_chshape();
  this->pi2=2*ppp.getPI();
  this->eps0=1.0/(ppp.getFPEPSI()*4*ppp.getPI());
  this->fpepsi=ppp.getFPEPSI();
}
	
	/* FFTVacuumField FFTVacuumField::getVF(
			int type,
			AtomSpecifier atoms,
			FFTGridType gt, FFTBoundaryCondition bc){

               try{
		if (FFTInteractionTypeCodes.lsType == type)
			return new FFTVacuumField_LS(atoms,  gt, bc);

                else if (FFTInteractionTypeCodes.rfType == type)
			return new FFTVacuumField_RF(atoms,  gt, bc);

		else{
		throw gromos::Exception("FFTVacuumField","Invalid interaction type code. Exiting.");
	} // endif
    } //end of try

                                catch (const gromos::Exception &e){
                                         cerr << e.what() << endl;
                                         exit(1);
                                }

}*/
