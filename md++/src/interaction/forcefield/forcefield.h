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
 * @file forcefield.h
 * the Forcefield class.
 */

#ifndef INCLUDED_FORCEFIELD_H
#define INCLUDED_FORCEFIELD_H

namespace configuration{
	class Configuration;
}
namespace topology{
	class Topology;
}
namespace simulation{
	class Simulation;
}

/**
 * @namespace interaction
 * namespace that contains the classes to
 * handle the interactions between the particles.
 * (energies, forces).
 */
namespace interaction
{
  /**
   * @class Forcefield
   * contains the specific interactions.
   * (clear does not call them (i guess) -- sorry, don't know what this means anymore)
   * Strategy Pattern.
   */
  class Forcefield : public std::vector<Interaction *>,
		     public algorithm::Algorithm
  {
  public:
    /**
     * Constructor
     */
    Forcefield() 
      : std::vector<Interaction *>(),
	algorithm::Algorithm("Forcefield") {}
    /**
     * Destructor
     */
    ~Forcefield();
    /**
     * initialise
     */
    virtual int init(topology::Topology & topo,
		     configuration::Configuration & conf,
		     simulation::Simulation & sim,
		     std::ostream & os = std::cout,
		     bool quiet = false);

    /**
     * calculate all interactions.
     */
    int calculate_interactions(topology::Topology & topo,
			       configuration::Configuration & conf,
			       simulation::Simulation & sim);


    /**
     * let the forcefield be used as an algorithm
     */
    int apply(topology::Topology & topo,
	      configuration::Configuration & conf,
	      simulation::Simulation & sim)
    {
      return calculate_interactions(topo, conf, sim);
    }
    
    virtual void print_timing(std::ostream & os);

    /**
     * const interaction accessor.
     */
    Interaction const * interaction(std::string name)const;

    /**
     * interaction accessor.
     */
    Interaction * interaction(std::string name);

    /**
     * GPU-native Interactions (NonBonded, bonded terms) atomicAdd into
     * the mirror's shared force/virial buffers and leave them resident
     * (calculate_interactions() zeros them once via sim.cuda().
     * zero_mirror_force(), not per-Interaction) -- exempting FORCE/
     * VIRIAL here would have kept Algorithm_Sequence::run()'s default
     * post-apply() invalidation from immediately erasing that
     * freshness. Narrowed further, to 0, once the real dependency this
     * mask was accidentally covering got fixed at its actual source
     * (see below) -- Forcefield never writes POS/VEL/BOX/
     * CONSTRAINT_FORCE/VIRIAL itself, so it has nothing to publish or
     * invalidate on its own behalf.
     *
     * History: narrowing this mask (PERFORMANCE.md's residual-resync-
     * cost investigation -- Leap_Frog_Velocity<gpuBackend> paid for a
     * POS/VEL resync every step the old mask forced, ~9.5s -> ~24s
     * across a 10000-step benchmark) used to break CUDA_Shake
     * (ubiquitin_gpu: "SHAKE error, vectors orthogonal" by step 3).
     * Root-caused via GROMOS_DEBUG_MIRROR-instrumented tracing:
     * Lattice_Shift_Tracker<gpuBackend> writes the periodic-image-
     * corrected position into the GPU mirror's current.pos and leaves
     * it GPU-dirty every step (its own gpu_mirror_touches() == 0). The
     * old (POS-including) Forcefield mask happened to flush that dirty
     * POS to the CPU-authoritative conf right before Leap_Frog_
     * Velocity<gpuBackend>'s conf.exchange_state() (a CPU-side pointer
     * swap of current<->old) -- without that flush, the swap relocated
     * the *stale* pre-shift CPU value into conf.old() while the
     * mirror's own swap (sim.cuda().exchange_mirror_state(), right
     * after) relocated the correct GPU value into its own old() half,
     * silently diverging conf.old() between CPU and GPU for any atom
     * that wrapped that step. CUDA_Shake reads conf.old().pos()
     * directly from the CPU array (bypassing the mirror), so it was
     * the one algorithm directly exposed to that divergence -- a
     * handful of wrapped atoms was enough to corrupt its reference
     * geometry into "vectors orthogonal" within a few steps.
     *
     * Fixed at the actual source instead of leaning on this mask:
     * Leap_Frog_Velocity<gpuBackend>::apply() (leap_frog_gpu.cc) now
     * calls sim.cuda().flush_gpu_dirty(conf, MIRROR_POS | MIRROR_VEL)
     * itself, immediately before conf.exchange_state() -- the
     * invariant ("no dirty GPU POS/VEL survives across an
     * exchange_state() swap") is now guaranteed locally, not as a
     * side effect of an unrelated algorithm's mask. Verified: full
     * ctest (31/32, only the known/documented aladip_cuda perturbation
     * gate failing), compute-sanitizer --tool memcheck clean on
     * ubiquitin_gpu (100 steps). See PERFORMANCE.md.
     */
    virtual unsigned gpu_mirror_touches() const override {
      return 0u;
    }

  protected:

  };
  
} // interaction

#endif
