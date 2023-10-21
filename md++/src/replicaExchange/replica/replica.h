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
 * @file replica.h
 * Holds all the information for a single replica
 */
#include <stdheader.h>
#include "_replica_Interface.h"

#ifndef REPLICA_H
#define REPLICA_H


namespace re {

  /**
   * replica Class
   * All the data members are held public, so Master/Slave classes have full 
   * control over it.
   */
  class replica: public replica_Interface {
  private:
    /**
     * copy constructor
     * don't allow copies because there's a bug in the copy constructor of configuration::configurtiaon
     */
    replica(const replica & r);
  public:
    /**
     * Constructor
     * Reads in all input parameters and assigns rank and ID of the replica
     * @param _args io::Argument, copy of the one initialized in repex_mpi.cc 
     * @param _ID integer, unique identifier of the replica
     * @param _rank integer, for knowing where to send the data
     */
    replica(io::Argument _args, int cont, int globalThreadID, simulation::MpiControl &replica_mpi_control);

    /**
     * runs an MD simulation
     */
    void run_MD();
    /**
     * init an MD simulation
     */
    void init();
    /**
     *
     */
    double calculateEnergies() override;

  };
}

#endif  /* REPLICA_H */
