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

/*
Highly Optimized Object-Oriented Molecular Dynamics (HOOMD) Open
Source Software License
Copyright (c) 2008 Ames Laboratory Iowa State University
All rights reserved.

Redistribution and use of HOOMD, in source and binary forms, with or
without modification, are permitted, provided that the following
conditions are met:

* Redistributions of source code must retain the above copyright notice,
this list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.

* Neither the name of the copyright holder nor the names HOOMD's
contributors may be used to endorse or promote products derived from this
software without specific prior written permission.

Disclaimer

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDER AND
CONTRIBUTORS ``AS IS''  AND ANY EXPRESS OR IMPLIED WARRANTIES,
INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY
AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. 

IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS  BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
THE POSSIBILITY OF SUCH DAMAGE.
*/

// $Id: WallData.h 1826 2009-04-27 22:18:31Z joaander $
// $URL: https://codeblue.umich.edu/hoomd-blue/svn/tags/hoomd-0.8.2/src/data_structures/WallData.h $
// Maintainer: joaander

/*! \file WallData.h
 	\brief Contains declarations for WallData.
 */

#include <math.h>
#include "ParticleData.h"

#include <boost/utility.hpp>

#ifndef __WALLDATA_H__
#define __WALLDATA_H__

//! Simple structure representing a single wall
/*! Walls are represented by an origin and a unit length normal.
	\ingroup data_structs
*/
struct Wall
	{
	//! Constructor
	/*!	\param ox Origin x-component
		\param oy Origin y-component
		\param oz Origin z-component
		\param nx Origin x-component
		\param ny Normal y-component
		\param nz Normal z-component
	*/
	Wall(Scalar ox=0.0, Scalar oy=0.0, Scalar oz=0.0, Scalar nx=1.0, Scalar ny=0.0, Scalar nz=0.0)
		: origin_x(ox), origin_y(oy), origin_z(oz)
		{
		// normalize nx,ny,nz
		Scalar len = sqrt(nx*nx + ny*ny + nz*nz);
		normal_x = nx / len;
		normal_y = ny / len;
		normal_z = nz / len;
		}

	Scalar origin_x;	//!< x-component of the origin
	Scalar origin_y;	//!< y-component of the origin
	Scalar origin_z;	//!< z-component of the origin

	Scalar normal_x;	//!< x-component of the normal
	Scalar normal_y;	//!< y-component of the normal
	Scalar normal_z;	//!< z-component of the normal
	};

//! Stores information about all the walls defined in the simulation
/*!	WallData is responsible for storing all of the walls in the simulation. 
	Walls are specified by the Wall struct and any number can be added. 

	On the CPU, walls can be accessed with getWall()

	An optimized data structure for the GPU will be written later.
	It will most likely take the form of a 2D texture.
	\ingroup data_structs
*/
class WallData : boost::noncopyable
	{
	public:
		//! Creates an empty structure with no walls
		WallData() : m_walls() {}
		
		//! Get the number of walls in the data
		/*! \return Number of walls
		*/
		unsigned int getNumWalls() const
			{
			return (unsigned int)m_walls.size();
			}
		
		//! Access a specific wall
		/*! \param idx Index of the wall to retrieve
			\return Wall stored at index \a idx
		*/
		const Wall& getWall(unsigned int idx) const
			{
			assert(idx < m_walls.size());
			return m_walls[idx];
			}
		
		//! Adds a wall to the data structure
		void addWall(const Wall& wall);

	private:
		//! Storage for the walls
		std::vector<Wall> m_walls;
	};

#endif
