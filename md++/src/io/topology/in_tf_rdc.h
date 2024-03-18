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
 * @file in_tf_rdc.h
 * read in a tensor-free RDC restraints file.
 */
/**
 * @page tfrdc tensor-free RDC restraints format
 * @date 25-04-2019
 *
 * A tensor-free RDC restraints specification file may contain the following
 * blocks:
 * - @ref title
 * - @ref tfrdcconversion
 * - @ref tfrdcresspec
 * - @ref tfrdcmolaxis
 */


#ifndef IN_TF_RDC_H
#define	IN_TF_RDC_H

#include "../instream.h"

namespace io {

  /**
   * @class In_Tfrdcresspec
   * reads in a tensor-free RDC restraints file
   */
  class In_Tfrdcresspec : public GInStream {

  private:
    double dish, disc;

  public:
    /**
     * Default constructor.
     */
    In_Tfrdcresspec() {}
    /**
     * Constructor.
     */
    In_Tfrdcresspec(std::istream& is) : GInStream(is) { readStream(); };
    /**
     * Read in a tensor-free RDC restraints file.
     */
    void read(topology::Topology &topo,
	      simulation::Simulation & sim,
	      std::ostream & os = std::cout);

    void read_TFRDCCONV(topology::Topology &topo,
                          simulation::Simulation &sim,
                          std::ostream & os);
    void read_TFRDCRESSPEC(topology::Topology &topo,
                          simulation::Simulation &sim,
                          std::ostream & os);
    void read_TFRDCMOLAXIS(topology::Topology &topo,
                          simulation::Simulation &sim,
                          std::ostream & os);
    /**
     * Maximum number of atoms that can be specified to define a virtual atom
     */
    static const unsigned int MAX_ATOMS = 4;
  };

} // io

#endif	/* IN_TF_RDC_H */
