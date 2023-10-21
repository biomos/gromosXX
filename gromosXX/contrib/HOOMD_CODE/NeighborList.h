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

// $Id: NeighborList.h 1983 2009-07-24 18:04:33Z akohlmey $
// $URL: https://codeblue.umich.edu/hoomd-blue/svn/tags/hoomd-0.8.2/src/computes/NeighborList.h $
// Maintainer: joaander

#ifdef WIN32
#pragma warning( push )
#pragma warning( disable : 4103 )
#endif

#include <boost/shared_ptr.hpp>
#include <boost/signals.hpp>
#include <vector>

#include "Compute.h"

//#incldue <string.h>
using std::string;

#ifdef ENABLE_CUDA
#include "NeighborList.cuh"
#endif

/*! \file NeighborList.h
	\brief Declares the NeighborList class
*/

#ifndef __NEIGHBORLIST_H__
#define __NEIGHBORLIST_H__

//! Computes a Neibhorlist from the particles
/*! Specification: A particle \c i is a neighbor of particle \c j if the distance between
	particle them is less than or equal to \c r_cut. The neighborlist for a given particle
	\c i includes all of these neighbors at a minimum. Other particles particles are included
	in the list: those up to \c r_max which includes a buffer distance so that the neighbor list
	doesn't need to be updated every step.
	
	There are two ways of storing this information. One is to store only half of the 
	neighbors (only those with i < j), and the other is to store all neighbors. There are
	potential tradeoffs between number of computations and memory access complexity for 
	each method. Since the goal of HOOMD is to explore all these possibilities, the 
	NeighborLists will support both of these modes via a switch: setStorageMode(); 
	
	In particular, classes such as LJForceCompute can work with
	either setting, full or half, but they are faster with the half setting. However,
	LJForceComputeGPU requires that the neighbor list storage mode is set to full.
	
	This base class supplys a dumb O(N^2) algorithm for generating this list. It is very 
	slow, but functional. For a more efficient implementation on the CPU, create a 
	BinnedNeighborList, which is a subclass of this class and is used through the 
	same interface documented here. Similarly, BinnedNeighborlistGPU will compute neighbor lists
	on the GPU.
	
	When compute() is called, the neighbor list is updated, but only if it needs to be. Checks
	are performed to see if any particle has moved more than half of the buffer distance, and
	only then is the list actually updated. This check can even be avoided for a number of time
	steps by calling setEvery(). If the caller wants to forces a full update, forceUpdate() 
	can be called before compute() to do so. Note that if the particle data is resorted,
	an update is automatically forced.
	
	Particles pairs can be excluded from the list by calling addExclusion(). Each particle can 
	exclude up to 4 others from appearing in its list. copyExclusionsFromBonds() adds an 
	exclusion for every bond specified in the ParticleData. \b NOTE: this is a one time thing.
	Any dynamic bonds added after calling copyExclusionsFromBonds() will not be added to the 
	neighborlist exclusions: this must be done by hand.
	
	The calculated neighbor list can be read using getList().
	
	\ingroup computes
*/
class NeighborList : public Compute
	{
	public:
		//! Simple enum for the storage modes
		enum storageMode
			{
			half,	//!< Only neighbors i,j are stored where i < j
			full	//!< All neighbors are stored
			};
#if !defined(LARGE_EXCLUSION_LIST)
		static const unsigned int MAX_NUM_EXCLUDED = 4;  //!< Maximum number of allowed exclusions per atom
#else
		static const unsigned int MAX_NUM_EXCLUDED = 16; //!< Maximum number of allowed exclusions per atom
#endif
		static const unsigned int EXCLUDE_EMPTY = 0xffffffff; //!< Signifies this element of the exclude list is empty
		//! Simple data structure for storing excluded tags
		/*! We only need to exclude up to 4 beads for typical molecular systems with bonds only, or bonds 
			and angles, if the molecules are linear. For regular systems with angles, dihedrals and/or impropers,
			HOOMD needs to be recompiled with -DLARGE_EXCLUSION_LIST. NOTE: this is a temporary hack until
			neighbor lists and exclusions are re-implemented to handle this in a better way.
			Initialized to empty for all elements. Empty is defined to be the maximum sized unsigned int.
			Since it isn't likely we will be performing 4 billion pariticle sims anytime soon, this is OK.
			
			\note The elements of this list are particle TAGS, not indices into the particle array. This is done
				to support particle list reordering.
		*/
		struct ExcludeList
			{
			//! Constructs an empty ExcludeList
			ExcludeList() : e1(EXCLUDE_EMPTY), e2(EXCLUDE_EMPTY), e3(EXCLUDE_EMPTY), e4(EXCLUDE_EMPTY) {}
			unsigned int e1; //!< Exclusion tag 1
			unsigned int e2; //!< Exclusion tag 2
			unsigned int e3; //!< Exclusion tag 3
			unsigned int e4; //!< Exclusion tag 4
			};
			
		//! Constructs the compute
		NeighborList(boost::shared_ptr<ParticleData> pdata, Scalar r_cut, Scalar r_buff);

		//! Destructor
		virtual ~NeighborList();

		//! Print statistics on the neighborlist
		virtual void printStats();

		//! Computes the NeighborList if it needs updating
		void compute(unsigned int timestep);
		
		//! Benchmark the neighbor list
		virtual double benchmark(unsigned int num_iters);

		//! Change the cuttoff radius
		void setRCut(Scalar r_cut, Scalar r_buff);

		//! Change how many timesteps before checking to see if the list should be rebuilt
		void setEvery(unsigned int every);

		//! Access the calculated neighbor list on the CPU
		virtual const std::vector< std::vector<unsigned int> >& getList();
		
		#ifdef ENABLE_CUDA
		//! Acquire the list on the GPU
		vector<gpu_nlist_array>& getListGPU();
		#endif
		
		//! Set the storage mode
		void setStorageMode(storageMode mode);
		
		//! Get the storage mode
		storageMode getStorageMode() { return m_storage_mode; }

		//! Exclude a pair of particles from being added to the neighbor list
		void addExclusion(unsigned int tag1, unsigned int tag2);
		
		//! Clear all existing exclusions
		void clearExclusions();
		
		//! Collect some statistics on exclusions.
		void countExclusions();
		
		//! Add an exclusion for every bond in the ParticleData
		void addExclusionsFromBonds();

		//! Add exclusions from angles
		void addExclusionsFromAngles();
		
		//! Add exclusions from dihedrals
		void addExclusionsFromDihedrals();
		
		//! Test if an exclusion has been made
		bool isExcluded(unsigned int tag1, unsigned int tag2);
		
		//! Add an exclusion for every 1,3 pair 
		void addOneThreeExclusionsFromTopology();
		
		//! Add an exclusion for every 1,4 pair
		void addOneFourExclusionsFromTopology();

		//! Forces a full update of the list on the next call to compute()
		void forceUpdate() { m_force_update = true; }
		
		//! Gives an estimate of the number of nearest neighbors per particle
		virtual Scalar estimateNNeigh();

	protected:
		Scalar m_r_cut;	//!< The cuttoff radius
		Scalar m_r_buff; //!< The buffer around the cuttoff

		std::vector< std::vector<unsigned int> > m_list; //!< The neighbor list itself
		storageMode m_storage_mode;	//!< The storage mode
		
		boost::signals::connection m_sort_connection;	//!< Connection to the ParticleData sort signal
		
		std::vector< ExcludeList > m_exclusions; //!< Stores particle exclusion lists BY TAG
#if defined(LARGE_EXCLUSION_LIST)
		std::vector< ExcludeList > m_exclusions2; //!< Stores particle exclusion lists BY TAG
		std::vector< ExcludeList > m_exclusions3; //!< Stores particle exclusion lists BY TAG
		std::vector< ExcludeList > m_exclusions4; //!< Stores particle exclusion lists BY TAG
#endif
		#ifdef ENABLE_CUDA
		//! Simple type for identifying where the most up to date particle data is
		enum DataLocation
			{
			cpu,	//!< Particle data was last modified on the CPU
			cpugpu,	//!< CPU and GPU contain identical data
			gpu		//!< Particle data was last modified on the GPU
			};
		
		DataLocation m_data_location;			//!< Where the neighborlist data currently lives
		vector<gpu_nlist_array> m_gpu_nlist;	//!< Stores pointers and dimensions of GPU data structures
		unsigned int *m_host_nlist;				//!< Stores a temporary copy of the neighbor list on the host
		unsigned int *m_host_n_neigh;			//!< Stores a temporary count of the number of neighbors on the host
		uint4 *m_host_exclusions;				//!< Stores a translated copy of the excluded particles on the host
#if defined(LARGE_EXCLUSION_LIST)
		uint4 *m_host_exclusions2;				//!< Stores a translated copy of the excluded particles on the host
		uint4 *m_host_exclusions3;				//!< Stores a translated copy of the excluded particles on the host
		uint4 *m_host_exclusions4;				//!< Stores a translated copy of the excluded particles on the host
#endif
		
		//! Helper function to move data from the host to the device
		void hostToDeviceCopy();
		//! Helper function to move data from the device to the host
		void deviceToHostCopy();
		//! Helper function to update the exclusion data on the device
		void updateExclusionData();
		//! Helper function to allocate data
		void allocateGPUData(int height);
		//! Helper function to free data
		void freeGPUData();
		#endif
	
		//! Performs the distance check
		virtual bool distanceCheck();
		
		//! Updates the previous position table for use in the next distance check
		virtual void setLastUpdatedPos();

		//! Builds the neighbor list
		virtual void buildNlist();
		
	private:
		int64_t m_updates;			//!< Number of times the neighbor list has been updated
		int64_t m_forced_updates;	//!< Number of times the neighbor list has been foribly updated
		int64_t m_dangerous_updates;	//!< Number of dangerous builds counted

		bool m_force_update;	//!< Flag to handle the forcing of neighborlist updates

		unsigned int m_last_updated_tstep; //!< Track the last time step we were updated
		unsigned int m_every; //!< No update checks will be performed until m_every steps after the last one

		Scalar *m_last_x;		//!< x coordinates of last updated particle positions
		Scalar *m_last_y;		//!< y coordinates of last updated particle positions
		Scalar *m_last_z;		//!< z coordinates of last updated particle positions
		
		//! Test if the list needs updating
		bool needsUpdating(unsigned int timestep);		

	};
	

#endif

#ifdef WIN32
#pragma warning( pop )
#endif
