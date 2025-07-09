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
 * @file remove_com_motion_gpu.h
 * remove center of mass translational and angular momentum (GPU variant)
 */

#pragma once


namespace algorithm
{
  /**
   * @class Remove_COM_Motion_GPU
   * implements GPU variant of centre of mass motion removal
   * input switches determine whether to remove
   * translational or angular centre of mass
   * momentum (or both).
   * For periodic boundary conditions, only
   * translational centre of mass motion is
   * removed.
   * It is recommended to remove it at every step.
   */
  class Remove_COM_Motion_GPU : public Remove_COM_Motion
  {
  public:
    /**
     * Constructor.
     */
    Remove_COM_Motion_GPU(gpu::CudaManager& cuda_manager,
                          std::ostream & os = std::cout,
                          const std::string& name = "RemoveCOMMotion") : Remove_COM_Motion(os, "RemoveCOMMotionGPU"),
                                                            os(os),
                                                            cuda_manager_(cuda_manager) {}

    /**
     * Destructor.
     */
    virtual ~Remove_COM_Motion_GPU() {}
    
    /**
     * apply COM removal.
     */
    virtual int apply(topology::Topology & topo,
		      configuration::Configuration & conf,
		      simulation::Simulation & sim);

    /**
     * init
     */
    virtual int init(topology::Topology &topo, 
		     configuration::Configuration &conf,
		     simulation::Simulation &sim,
		     std::ostream &os = std::cout,
		     bool quiet = false);

    /**
     * calculate and remove translational centre of mass motion
     */
    double remove_com_translation(topology::Topology & topo,
				  configuration::Configuration & conf,
				  simulation::Simulation & sim,
				  bool remove_trans = true);
    
    /**
     * calculate and remove angular centre of mass motion
     */
    double remove_com_rotation(topology::Topology & topo,
			       configuration::Configuration & conf,
			       simulation::Simulation & sim,
			       bool remove_rot = true);

    /**
     * add centre of mass rotation
     */
    double add_com_rotation(topology::Topology & topo,
			    configuration::Configuration & conf,
			    simulation::Simulation & sim,
			    math::Vec com_L);
    

  protected:
    std::ostream & os;
    gpu::CudaManager& cuda_manager_;
  };
  
} //algorithm
