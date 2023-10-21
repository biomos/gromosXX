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
 * @file check_parameter.h
 * check parameters
 */

#ifndef INCLUDED_CHECK_PARAMETER_H
#define INCLUDED_CHECK_PARAMETER_H

namespace io
{
  /**
   * cross checks on the parameters
   */
  int check_parameter(simulation::Simulation &sim);
  
  /**
   * does basic cross checks on parameters from different blocks
   * cross-checks of parameters within one block should be done in the 
   * read_XXX functions in @ref In_Parameter
   */
  int simple_crosschecks(simulation::Simulation &sim);
  
  /**
   * does basic cross checks on the parameters
   * checks each feature against all others
   * every allowed combination has to be unlocked
   */
  int check_features(simulation::Simulation &sim);
}

#endif
