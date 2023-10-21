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
 * @file error.h
 * define error values.
 */

#ifndef INCLUDED_ERROR_H
#define INCLUDED_ERROR_H

// errors
#define E_UNSPECIFIED 1
#define E_NAN 2
#define E_NOT_IMPLEMENTED 3
#define E_ILLEGAL 4
#define E_SHAKE_FAILURE 10
#define E_SHAKE_FAILURE_SOLUTE 11
#define E_SHAKE_FAILURE_SOLVENT 12
#define E_TYPE_NOT_IMPLEMENTED 13
#define E_INPUT_ERROR 14
#define E_USAGE 15
#define E_BOUNDARY_ERROR 16

// termination
#define E_MINIMUM_REACHED 1000

#endif
