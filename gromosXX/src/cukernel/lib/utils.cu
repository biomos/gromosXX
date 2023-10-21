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
 * @file utils.cu
 * implementation of utils
 */

#include <cuda.h>
#include <iostream>
#include "utils.h"

int cudakernel::checkError(const char * err_msg) {
  cudaError_t error = cudaGetLastError();
  if (error != cudaSuccess)
    std::cout << "CUDA-ERROR " << err_msg << ": " << cudaGetErrorString(error) << std::endl;

  return (int) error;
}

