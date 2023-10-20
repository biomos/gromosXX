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
 * @file template_split.h
 * split function call based on template paramters
 */

////////////////////////////////////////////////////////////////////////////////
/**
 * call a function with the appropriate boundary as template parameter
 */
#define SPLIT_BOUNDARY(f, ...) \
switch(conf.boundary_type){ \
  case math::vacuum : f<math::vacuum>(__VA_ARGS__); break; \
  case math::rectangular : f<math::rectangular>(__VA_ARGS__); break; \
  case math::truncoct : \
  case math::triclinic : f<math::triclinic>(__VA_ARGS__); break; \
  default: io::messages.add("wrong boundary type", "template_split", io::message::error); \
} \

////////////////////////////////////////////////////////////////////////////////
/**
 * call a function with the appropriate boundary as template parameter
 */
#define SPLIT_MY_BOUNDARY(bound, f, ...) \
switch(bound){ \
  case math::vacuum : f<math::vacuum>(__VA_ARGS__); break; \
  case math::rectangular : f<math::rectangular>(__VA_ARGS__); break; \
  case math::truncoct : \
  case math::triclinic : f<math::triclinic>(__VA_ARGS__); break; \
  default: io::messages.add("wrong boundary type", "template_split", io::message::error); \
} \

////////////////////////////////////////////////////////////////////////////////

/**
 * call a function with the appropriate values for virial and boundary
 * as template parameters.
 */
#define SPLIT_VIRIAL_BOUNDARY(f, ...) \
switch(conf.boundary_type){ \
  case math::vacuum : f<math::vacuum, math::no_virial>(__VA_ARGS__); break; \
  case math::rectangular : \
    switch(sim.param().pcouple.virial){ \
      case math::no_virial : f<math::rectangular, math::no_virial>(__VA_ARGS__); break; \
      case math::molecular_virial : f<math::rectangular, math::molecular_virial>(__VA_ARGS__); break; \
      case math::atomic_virial : f<math::rectangular, math::atomic_virial>(__VA_ARGS__); break; \
      default: io::messages.add("wrong virial type", "template_split", io::message::error); \
    } \
    break; \
  case math::truncoct : \
  case math::triclinic : \
    switch(sim.param().pcouple.virial){ \
      case math::no_virial : f<math::triclinic, math::no_virial>(__VA_ARGS__); break; \
      case math::molecular_virial : f<math::triclinic, math::molecular_virial>(__VA_ARGS__); break; \
      case math::atomic_virial : f<math::triclinic, math::atomic_virial>(__VA_ARGS__); break; \
      default: io::messages.add("wrong virial type", "template_split", io::message::error); \
    } \
    break; \
  default: io::messages.add("wrong boundary type", "template_split", io::message::error); \
} \

////////////////////////////////////////////////////////////////////////////////


