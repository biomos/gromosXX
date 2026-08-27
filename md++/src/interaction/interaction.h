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
 * @file src/interaction/interaction.h
 * the interaction interface.
 */

#ifndef INCLUDED_INTERACTION_H
#define INCLUDED_INTERACTION_H

namespace configuration{
  class Configuration;
}
namespace topology{
  class Topology;
}
namespace simulation{
  class Simulation;
}
namespace util {
  class Algorithm_Timer;
}

namespace interaction
{
  /**
   * @class Interaction
   * @interface Interaction
   * declares the interaction interface.
   */
  class Interaction
  {
  public:
    /**
     * Constructor.
     */
    Interaction(std::string name) : name(name), m_timer(name) {};
    /**
     * Destructor.
     */
    virtual ~Interaction(){};
    /**
     * the name of the interaction.
     * can be used to identify a special class.
     */
    std::string name;
    /**
     * initialise
     */
    virtual int init(topology::Topology & topo,
		     configuration::Configuration & conf,
		     simulation::Simulation & sim,
		     std::ostream & os = std::cout,
		     bool quiet = false) = 0;
    // { return 0; }
    /**
     * calculate the interactions.
     */
    virtual int calculate_interactions(topology::Topology & topo,
				       configuration::Configuration & conf,
				       simulation::Simulation & sim) = 0;

    /**
     * @brief Does this Interaction touch conf.current().force through
     * the plain CPU array (conf.current().force(i) += ..., a read-
     * modify-write, whether or not it explicitly reads the total first)
     * rather than writing into the GPU-resident mirror itself? Default
     * true -- the overwhelming majority of Interactions (every special/
     * restraint force, QMMM, every non-GPU-native bonded term, and
     * Molecular_Virial_Interaction's correction) fall into this
     * category. Forcefield::calculate_interactions() checks this flag
     * and publishes the GPU mirror's force before calling such an
     * Interaction -- without it, a CPU-side += would silently operate
     * on a stale/zero array (GPU-native terms leave force GPU-resident,
     * gpu::MIRROR_FORCE marked dirty via mark_gpu_dirty(), not synced
     * back every call), and worse, that addition would later be wiped
     * out entirely when the mirror's force is eventually published
     * (flush_gpu_dirty() *overwrites* conf.current().force from the
     * mirror, it doesn't merge). Found the hard way: POSITIONRES active
     * alongside GPU-native bonded/nonbonded silently dropped the
     * position-restraint contribution every step, producing a slow,
     * escalating LINCS-rotation divergence over thousands of steps that
     * a short/step-0-only test never caught.
     *
     * Only the known GPU-native force writers override this to false:
     * CUDA_Angle_Interaction, CUDA_Dihedral_Interaction, CUDA_Improper_
     * Dihedral_Interaction, CUDA_Quartic_Bond_Interaction, CUDA_
     * Nonbonded_Interaction, CUDA_Position_Restraint_Interaction --
     * each atomicAdd's into the mirror directly and must NOT trigger a
     * flush before its own call (that would force a premature,
     * wasteful publish of whatever an earlier GPU-native term already
     * wrote this step, defeating GPU residency). Deliberately the
     * minority list: new CPU-only Interactions need no changes to be
     * safe by default; only a new GPU-native force writer needs to
     * remember to opt out.
     */
    virtual bool needs_fresh_cpu_force() const { return true; }

    /**
     * @brief Same contract as needs_fresh_cpu_force(), for virial_tensor.
     * Default false: unlike force, most CPU-only Interactions never read
     * virial_tensor mid-step, so an unconditional flush here would be
     * pure waste for the common case. Molecular_Virial_Interaction is the
     * one real consumer -- it reads the atomic virial accumulated by
     * GPU-native bonded/nonbonded terms (their mark_gpu_dirty(MIRROR_
     * VIRIAL) calls) and corrects it to a molecular virial in place, so
     * it must see the true accumulated total, not the zeroed-this-step
     * mirror value that's still sitting on the CPU side otherwise. Found
     * the same way as the POSITIONRES force-drop bug: GPU virial ended
     * up wildly wrong (raw atomic virial, symmetric) vs. CPU's correct
     * molecular virial (asymmetric, smaller) despite matching kinetic
     * energy -- because nothing flushed MIRROR_VIRIAL before this
     * Interaction's turn (Pressure_Calculation's own MIRROR_VIRIAL flush
     * happens later, outside Forcefield entirely).
     */
    virtual bool needs_fresh_cpu_virial() const { return false; }

    /**
     * @brief Does this Interaction have a real GPU-native implementation
     * (writes force directly into the GPU-resident mirror, no per-call
     * CPU round trip)? Default false. Used only for diagnostics --
     * Forcefield::init() warns once, when the run's accelerator is
     * gpu_cuda, about every child Interaction that returns false here,
     * so a mixed CPU/GPU configuration (e.g. a restraint type with no
     * CUDA port yet) is visible up front instead of only showing up as
     * an unexplained slowdown (or, before needs_fresh_cpu_force()
     * above existed, silently wrong dynamics). Override to true in the
     * same classes that override needs_fresh_cpu_force() to false.
     */
    virtual bool is_gpu_native() const { return false; }

    /**
     * timing information.
     */
    virtual void print_timing(std::ostream & os)
    {
      m_timer.print(os);
    }

  protected:
    /**
     * store time used in algorithm.
     */
    util::Algorithm_Timer m_timer;

  };  
  
} // interaction

#endif
