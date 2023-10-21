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
 * @file repex_mpi.h
 * some definitions for MPI communication - types etc.
 * Modified June 18, 2021 - bschroed, srieder
 */

#ifndef REPEX_MPI_H
#define	REPEX_MPI_H

#ifdef XXMPI
#include <mpi.h>
#endif

#ifdef REPEX_MPI
#define EXTERN
#else
#define EXTERN extern
#endif

/*
 * define some MPI datatypes
 */
//  EXTERN MPI_Datatype MPI_INITVAL;
#ifdef XXMPI
    EXTERN MPI_Datatype MPI_VARRAY;
    EXTERN MPI_Datatype MPI_CONFINFO;
    EXTERN MPI_Datatype MPI_BOX;
    EXTERN MPI_Datatype MPI_REPINFO;
    EXTERN MPI_Datatype MPI_EDSINFO;
#endif


enum {
  INIT, CONF, POS, POSV, VEL, LATTSHIFTS, STOCHINT, BOX, ANGLES, DF, ENERGIES, SWITCHENERGIES, SENDCOORDS, REPINFO, EDSINFO,
};

namespace re {
  /**
   * @struct repInfo
   * holds the information that is sent to master after MD run
   */
  struct repInfo {
    int switched;
    int run;
    int partner;
    double epot;
    double epot_partner;
    double probability;
  };

}

/**
 * @typedef ID_t
 * ID type
 */
typedef unsigned int ID_t;
/**
 * @typedef rank_t
 * rank_type
 */
typedef unsigned int rank_t;

#endif	/* REPEX_MPI_H */
