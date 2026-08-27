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
 * @file configuration_struct.cu
 * Implementation of the light-weight configuration struct for GPU
 */

#include "stdheader.h"

#include "configuration/configuration_global.h"

#include "algorithm/algorithm.h"
#include "topology/topology.h"
#include "configuration/configuration.h"
#include "configuration/mesh.h"
#include "configuration/influence_function.h"
#include "simulation/simulation.h"
#include "simulation/multibath.h"
#include "simulation/parameter.h"

#include "math/periodicity.h"
#include "math/boundary_checks.h"
#include "util/template_split.h"

#include "configuration_struct.h"

#include "gpu/cuda/utils.h"

void gpu::Configuration::copy_pos_vel_to_device(const configuration::Configuration& conf) {
    const size_t n = conf.current().pos.size();
    current.pos.resize(n);
    current.vel.resize(n);
    old.pos.resize(n);
    old.vel.resize(n);
    for (size_t i = 0; i < n; ++i) {
        current.pos[i] = static_cast<FPL3_TYPE>(conf.current().pos(i));
        current.vel[i] = static_cast<FPL3_TYPE>(conf.current().vel(i));
        old.pos[i]     = static_cast<FPL3_TYPE>(conf.old().pos(i));
        old.vel[i]     = static_cast<FPL3_TYPE>(conf.old().vel(i));
    }
}

void gpu::Configuration::copy_forces_from_device(configuration::Configuration& conf) {
    CUDA_CHECK(cudaDeviceSynchronize());
    const size_t n = current.force.size();
    for (size_t i = 0; i < n; ++i) {
        const FPH3_TYPE& f = current.force[i];
        conf.current().force(i) = math::Vec(f.x, f.y, f.z);
    }
}

void gpu::Configuration::copy_pos_vel_from_device(configuration::Configuration& conf) {
    CUDA_CHECK(cudaDeviceSynchronize());
    const size_t n = current.pos.size();
    for (size_t i = 0; i < n; ++i) {
        const FPL3_TYPE& p = current.pos[i];
        const FPL3_TYPE& v = current.vel[i];
        conf.current().pos(i) = math::Vec(p.x, p.y, p.z);
        conf.current().vel(i) = math::Vec(v.x, v.y, v.z);
    }
}

namespace {
    // math::Matrix (9 doubles) <-> FPH9_TYPE element-wise conversion --
    // a raw memcpy is wrong whenever FPH_TYPE != double (sizes differ,
    // e.g. FP_PRECISION==1 where high==float too), which is exactly the
    // bug this replaces below. virial_tensor/kinetic_energy_tensor/
    // pressure_tensor are FPH9_TYPE (accumulator targets, see
    // configuration_struct.h) -- under FP_PRECISION==2/3 this is a
    // genuine double, matching math::Matrix exactly (no precision loss
    // in this specific conversion, unlike the FPL9 case).
    FPH9_TYPE matrix_to_fph9(const math::Matrix & m) {
        FPH9_TYPE r;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                r(i, j) = static_cast<FPH_TYPE>(m(i, j));
        return r;
    }
    math::Matrix fph9_to_matrix(const FPH9_TYPE & r) {
        math::Matrix m;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                m(i, j) = static_cast<double>(r(i, j));
        return m;
    }
}

void gpu::Configuration::copy_constraint_data_from_device(configuration::Configuration& conf) {
    CUDA_CHECK(cudaDeviceSynchronize());

    const size_t n = current.constraint_force.size();
    for (size_t i = 0; i < n; ++i) {
        const FPH3_TYPE& cf_c = current.constraint_force[i];
        const FPH3_TYPE& cf_o = old.constraint_force[i];
        conf.current().constraint_force(i) = math::Vec(cf_c.x, cf_c.y, cf_c.z);
        conf.old().constraint_force(i) = math::Vec(cf_o.x, cf_o.y, cf_o.z);
    }

    FPH9_TYPE vt_c, vt_o;
    CUDA_CHECK(cudaMemcpy(&vt_c, current.virial_tensor, sizeof(FPH9_TYPE), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(&vt_o, old.virial_tensor, sizeof(FPH9_TYPE), cudaMemcpyDeviceToHost));
    conf.current().virial_tensor = fph9_to_matrix(vt_c);
    conf.old().virial_tensor = fph9_to_matrix(vt_o);
}

void gpu::Configuration::copy_virial_to_device(const configuration::Configuration& conf) {
    const FPH9_TYPE vt_c = matrix_to_fph9(conf.current().virial_tensor);
    const FPH9_TYPE vt_o = matrix_to_fph9(conf.old().virial_tensor);
    CUDA_CHECK(cudaMemcpy(current.virial_tensor, &vt_c, sizeof(FPH9_TYPE), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(old.virial_tensor, &vt_o, sizeof(FPH9_TYPE), cudaMemcpyHostToDevice));
}

void gpu::Configuration::copy_lattice_shifts_to_device(const configuration::Configuration& conf) {
    const size_t n = conf.special().lattice_shifts.size();
    lattice_shifts.resize(n);
    for (size_t i = 0; i < n; ++i)
        lattice_shifts[i] = static_cast<FPL3_TYPE>(conf.special().lattice_shifts(i));
}

void gpu::Configuration::copy_lattice_shifts_from_device(configuration::Configuration& conf) {
    CUDA_CHECK(cudaDeviceSynchronize());
    const size_t n = lattice_shifts.size();
    for (size_t i = 0; i < n; ++i) {
        const FPL3_TYPE& s = lattice_shifts[i];
        conf.special().lattice_shifts(i) = math::Vec(s.x, s.y, s.z);
    }
}

void gpu::Configuration::copy_to_device(configuration::Configuration& conf) {
    CUDA_CHECK_ERROR("At gpu::Configuration::copy_to_device");
    const size_t num_atoms = conf.current().pos.size();
    using Vec = typename decltype(conf.current().pos)::value_type;

    static_assert(std::is_convertible<Vec, FPL3_TYPE>::value,
                  "Vec must be convertible to FPL3_TYPE");
    static_assert(std::is_convertible<Vec, FPH3_TYPE>::value,
                  "Vec must be convertible to FPH3_TYPE");

    // Two casts, not one: pos/vel/lattice_shifts are FPL (state, low
    // precision is fine); force/constraint_force are FPH (accumulator
    // targets -- see configuration_struct.h's doc comment on why).
    auto convert_and_copy_low = [num_atoms](const auto& src, auto& dst) {
        CUDA_CHECK_ERROR("Before resize");
        dst.resize(num_atoms);
        CUDA_CHECK_ERROR("After resize");
        for (size_t i = 0; i < num_atoms; ++i)
            dst[i] = static_cast<FPL3_TYPE>(src[i]);
    };
    auto convert_and_copy_high = [num_atoms](const auto& src, auto& dst) {
        CUDA_CHECK_ERROR("Before resize");
        dst.resize(num_atoms);
        CUDA_CHECK_ERROR("After resize");
        for (size_t i = 0; i < num_atoms; ++i)
            dst[i] = static_cast<FPH3_TYPE>(src[i]);
    };

    // Current state
    convert_and_copy_low(conf.current().pos, current.pos);
    convert_and_copy_low(conf.current().vel, current.vel);
    convert_and_copy_high(conf.current().force, current.force);
    convert_and_copy_high(conf.current().constraint_force, current.constraint_force);

    // Old state
    convert_and_copy_low(conf.old().pos, old.pos);
    convert_and_copy_low(conf.old().vel, old.vel);
    convert_and_copy_high(conf.old().force, old.force);
    convert_and_copy_high(conf.old().constraint_force, old.constraint_force);

    // Persistent, non-cycling (see configuration_struct.h's doc comment)
    convert_and_copy_low(conf.special().lattice_shifts, lattice_shifts);


    CUDA_CHECK_ERROR("Before tensors");
    // Box is left as a raw memcpy (pre-existing, separate issue --
    // gpu::Box is also FPL_TYPE-based, so this has the same size-
    // mismatch bug the tensor copies below used to have; not fixed
    // here since nothing currently reads the mirror's .box at all --
    // every kernel takes math::Box by value at launch instead. See
    // KNOWN_ISSUES.md.
    CUDA_CHECK(cudaMemcpy(current.box, &conf.current().box, sizeof(Box), cudaMemcpyHostToDevice));
    // Element-wise convert, not a raw memcpy: math::Matrix is 9
    // doubles (72 bytes); FPH9_TYPE matches exactly under
    // FP_PRECISION 2/3 (high==double) but would still mismatch under
    // FP_PRECISION==1 (high==float, 36 bytes) -- a straight memcpy
    // would silently truncate/misinterpret in that case. Found while
    // wiring the constraint algorithms onto this mirror's virial_
    // tensor field, which made this bug live for the first time
    // (previously nothing read these three tensor fields at all).
    const FPH9_TYPE vt_c = matrix_to_fph9(conf.current().virial_tensor);
    const FPH9_TYPE ket_c = matrix_to_fph9(conf.current().kinetic_energy_tensor);
    const FPH9_TYPE pt_c = matrix_to_fph9(conf.current().pressure_tensor);
    CUDA_CHECK(cudaMemcpy(current.virial_tensor, &vt_c, sizeof(FPH9_TYPE), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(current.kinetic_energy_tensor, &ket_c, sizeof(FPH9_TYPE), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(current.pressure_tensor, &pt_c, sizeof(FPH9_TYPE), cudaMemcpyHostToDevice));

    CUDA_CHECK(cudaMemcpy(old.box, &conf.old().box, sizeof(Box), cudaMemcpyHostToDevice));
    const FPH9_TYPE vt_o = matrix_to_fph9(conf.old().virial_tensor);
    const FPH9_TYPE ket_o = matrix_to_fph9(conf.old().kinetic_energy_tensor);
    const FPH9_TYPE pt_o = matrix_to_fph9(conf.old().pressure_tensor);
    CUDA_CHECK(cudaMemcpy(old.virial_tensor, &vt_o, sizeof(FPH9_TYPE), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(old.kinetic_energy_tensor, &ket_o, sizeof(FPH9_TYPE), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(old.pressure_tensor, &pt_o, sizeof(FPH9_TYPE), cudaMemcpyHostToDevice));
    CUDA_CHECK_ERROR("After tensors");
}