/*
 * File:   pairlist.h
 * Author: Axel && Mel
 *
 * Created on October 20, 2011, 11:34 AM
 */
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <map>
#include <set>
#include <stdexcept>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <typeinfo>
#include <list>
#include </usr/local/cuda/include/cuda.h>
#include </usr/local/cuda/include/driver_types.h>
#include "../math/gmath.h"
#include "gromos_cuda.h"
#ifndef GROMOS_CUDA_PAIRLIST_H
#define	GROMOS_CUDA_PAIRLIST_H

namespace gcuda {

struct devices {
    int number_devices;
    std::vector<cudaDeviceProp> property;
    bool error;
};

extern "C" int has_gpu ();

extern "C" devices get_dev_properties (int);

//update atomic function

extern "C" void update_atomic(gcuda::memmap &gpumemory, double* h_positions, math::Box box, unsigned int nratoms, unsigned int nr_solute_atoms,
                          gcuda::TwoDArray<unsigned int> *h_results_solute_short,
                          gcuda::TwoDArray<unsigned int> *h_results_solute_long,
                          gcuda::TwoDArray<unsigned int> *h_results_solvent_short,
                          gcuda::TwoDArray<unsigned int> *h_results_solvent_long);

//update chargegroup function
extern "C" void update_cg(gcuda::memmap &gpumemory, double* h_positions, math::Box box, unsigned int nratoms, unsigned int nr_solute_atoms, unsigned int nrcg,
                          gcuda::TwoDArray<unsigned int> *h_results_solute_short,
                          gcuda::TwoDArray<unsigned int> *h_results_solute_long,
                          gcuda::TwoDArray<unsigned int> *h_results_solvent_short, 
                          gcuda::TwoDArray<unsigned int> *h_results_solvent_long);

// initialise the GPU memory
extern "C" void initcuda_atomic (gcuda::memmap &gpumemory, unsigned int nratoms, unsigned int nr_solute_atoms, double *cutoff, size_t host_results_pitch);
extern "C" void initcuda_cg (gcuda::memmap &gpumemory, unsigned int nratoms, unsigned int nr_solute_atoms, unsigned int nrcg, unsigned int nrsolutecg, double *cutoff, size_t host_results_pitch);

//free the GPU memory
extern "C" void freecuda (gcuda::memmap &gpumemory, bool atomic);

} // namespace

#endif	/* GROMOS_CUDA_PAIRLIST_H */
