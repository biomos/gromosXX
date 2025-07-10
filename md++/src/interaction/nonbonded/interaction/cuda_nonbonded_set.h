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

/**
 * @file cuda_nonbonded_set.h
 * the non bonded interactions for a set of atoms:
 * Lennard-Jones and Coulomb interactions
 */

#ifndef INCLUDED_CUDA_NONBONDED_SET_H
#define INCLUDED_CUDA_NONBONDED_SET_H

#include "pairlist.h"
#include "storage.h"
#include "nonbonded_outerloop.h"
#include "nonbonded_set_interface.h"
#include "nonbonded_set.h"
#include <util/cycle_thread.h>
#ifndef HAVE_LIBCUDART
#define gpu_status void
#endif


namespace interaction {

  /**
   * @class CUDA_Nonbonded_Set
   * calculates the nonbonded interactions.
   */
  class CUDA_Nonbonded_Set : public Nonbonded_Set/*, public util::CycleThread*/ {
  public:

    /**
     * Constructor.
     */
    CUDA_Nonbonded_Set(Pairlist_Algorithm & pairlist_alg, Nonbonded_Parameter & param,
            int rank, int num_threads);

    /**
     * Destructor
     */
    virtual ~CUDA_Nonbonded_Set();

    /**
     * initialize some things
     */
    virtual int init(topology::Topology const & topo,
            configuration::Configuration const & conf,
            simulation::Simulation const & sim,
            std::ostream & os = std::cout,
            bool quiet = false);

    /**
     * calculate the interactions.
     */
    virtual int calculate_interactions(topology::Topology & topo,
            configuration::Configuration & conf,
            simulation::Simulation & sim);

    Storage & shortrange_storage() {
      return m_storage;
    }

    Storage & longrange_storage() {
      return m_longrange_storage;
    }
    /**
     * maximal number of neighbors per atom in shortrange pairlist
     */
    unsigned int estNeigh_short;
    /**
     * maximal number of neighbors per atom in longrange pairlist
     */
    unsigned int estNeigh_long;

  private:
    /**
     * Pointer to the gpu status
     */
    // gpu_status * gpu_stat;
    /**
     * Pointer to the topology
     */
    topology::Topology * mytopo;
    /**
     * Pointer to the configuration
     */
    configuration::Configuration * myconf;
    /**
     * Pointer to the simulation
     */
    simulation::Simulation * mysim;
    /**
     * Nonbonded parameters
     */
    Nonbonded_Parameter * m_parameter;
    /**
     * Error integer to communicate problems between threads
     */
    int error;

    /**
     * Calculates the contstants needed for
     * further calculation
     */
    //virtual inline void init_run();
    /**
     * contains the calculations, executed at every step
     */
    //virtual void cycle();
    virtual inline void calculate();
    /**
     * Clean up
     */
    virtual inline void end_run();
    

  protected:
    Storage m_longrange_storage;
  };

} // interaction

#endif
