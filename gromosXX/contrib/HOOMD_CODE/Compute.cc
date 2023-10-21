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

// $Id: Compute.cc 1826 2009-04-27 22:18:31Z joaander $
// $URL: https://codeblue.umich.edu/hoomd-blue/svn/tags/hoomd-0.8.2/src/computes/Compute.cc $
// Maintainer: joaander

#ifdef WIN32
#pragma warning( push )
#pragma warning( disable : 4244 )
#endif

#include <iostream>
#include <stdexcept>
using namespace std;

#include "Compute.h"

/*! \file Compute.cc
	\brief Contains code for the Compute class
*/

/*! \param pdata Particle data this compute will act on. Must not be NULL.
	\post The Compute is constructed with the given particle data and a NULL profiler.
*/
Compute::Compute(boost::shared_ptr<ParticleData> pdata) : m_pdata(pdata), exec_conf(pdata->getExecConf()), m_last_computed(0), m_first_compute(true)
	{
	// sanity check
	assert(pdata);
	assert(pdata->getN() > 0);
	}
	
/*! \param num_iters Number of iterations to average for the benchmark
	\returns Milliseconds of execution time per calculation
	Derived classes can optionally implement this method. */
double Compute::benchmark(unsigned int num_iters)
	{
	cerr << endl << "***Error! This compute doesn't support benchmarking" << endl << endl;
	throw runtime_error("Error benchmarking compute");
	return 0.0;
	}
		
/*! It is useful for the user to know where computation time is spent, so all Computes
	should profile themselves. This method sets the profiler for them to use.
	This method does not need to be called, as Computes will not profile themselves
	on a NULL profiler
	\param prof Pointer to a profiler for the compute to use. Set to NULL 
		(boost::shared_ptr<Profiler>()) to stop the 
		analyzer from profiling itself.
	\note Derived classes MUST check if m_prof is set before calling any profiler methods.
*/
void Compute::setProfiler(boost::shared_ptr<Profiler> prof)
	{
	m_prof = prof;
	}

/*! \param timestep Current time step
 	\returns true if computations should be performed, false if they have already been done
 		at this \a timestep.
 	\note This method is designed to only be called once per call to compute() like so:
\code
void SomeClass::compute(unsgned int timestep)
	{
	if (!shouldCompute(timestep))
		return;
	... do compute tasks
	}
\endcode
*/
bool Compute::shouldCompute(unsigned int timestep)
	{
	// handle case where no computation has been performed yet
	if (m_first_compute)
		{
		m_first_compute = false;
		m_last_computed = timestep;
		return true;
		}

	// otherwise, we update if the last computed timestep is less than the current
	if (m_last_computed != timestep)
		{	
		m_last_computed = timestep;
		return true;
		}
	
	// failing the above, we perform no computation
	return false;
	}
	


#ifdef WIN32
#pragma warning( pop )
#endif
