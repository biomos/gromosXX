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
 * @file backend.h
 * tags and traits for CPU and GPU variants of algorithms
 */
#pragma once

namespace util
{
    #ifdef USE_CUDA
        inline constexpr bool cuda_enabled = true;
    #else
        inline constexpr bool cuda_enabled = false;
    #endif
    /**
     * @brief 
     *
     */
    struct cpuBackend {};
    /**
     * @brief
     *
     */
    struct gpuBackend {};

    /**
     * @brief Checks if the given backend is allowed (if USE_CUDA).
     *
     * @tparam Backend The backend to check.
     * @return true If the backend is a valid algorithm backend.
     * @return false Otherwise.
     */
    template <typename Backend>
    struct backend_is_available
    {
        static constexpr bool value =
            std::is_same_v<Backend, cpuBackend> ||
#ifdef USE_CUDA
            std::is_same_v<Backend, gpuBackend>;
#else
            false;
#endif
    };

    // Traits definition

    /**
     * @brief Checks if the given algorithm template has a specific backend implementation.
     * 
     * @tparam AlgT The algorithm template to check.
     * @tparam Backend The backend to check.
     * @return true_type If the algorithm template has a specific backend implementation.
     * @return false_type Otherwise.
     */
    template <template <typename> class AlgT, typename Backend, typename = void>
    struct has_backend : std::false_type {};
    template <template <typename> class AlgT, typename Backend>
    struct has_backend<
        AlgT,
        Backend,
        std::enable_if_t<
            AlgT<Backend>::template is_supported_backend<Backend> &&
            (!std::is_same_v<Backend, util::gpuBackend> || cuda_enabled)
        >
    > : std::true_type {};

    /**
     * @brief Checks if the given algorithm template has a CPU backend implementation.
     * 
     * @tparam AlgT The algorithm template to check.
     * @return true_type If the algorithm template has a CPU backend implementation.
     * @return false_type Otherwise.
     */
    template <template <typename> class AlgT>
    struct has_cpu_backend : has_backend<AlgT, util::cpuBackend> {};

    /**
     * @brief Checks if the given algorithm template has a CPU backend implementation.
     * 
     * @tparam AlgT 
     * @return true If the algorithm template has a CPU backend implementation.
     * @return false Otherwise.
     */
    template <template <typename> class AlgT>
    constexpr bool has_cpu_backend_v =
        has_cpu_backend<AlgT>::value;


    /**
     * @brief Checks if the given algorithm template has a GPU backend implementation.
     * 
     * @tparam AlgT The algorithm template to check.
     * @return true_type If the algorithm template has a GPU backend implementation and CUDA is activated
     * @return false_type Otherwise.
     */
    template <template <typename> class AlgT>
    struct has_gpu_backend : has_backend<AlgT, util::gpuBackend> {};

    /**
     * @brief Checks if the given algorithm template has a GPU backend implementation.
     * 
     * @tparam AlgT 
     * @return true If the algorithm template has a Gpu backend implementation.
     * @return false Otherwise.
     */
    template <template <typename> class AlgT>
    constexpr bool has_gpu_backend_v =
        has_gpu_backend<AlgT>::value;
}