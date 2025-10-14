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
 * @file cumath.h 
 * math operations for host and device
 */

#pragma once

/**
 * calculates the center of geometry
 * @param pos array of positions
 * @param begin start index
 * @param end end index
 * @return center of geometry
 */
template <typename ArrT, typename SizeType,
          typename = std::enable_if_t<
              std::is_integral<SizeType>::value &&
              std::is_convertible<
                  decltype(std::declval<ArrT&>()(std::declval<SizeType>())),
                  typename ArrT::value_type
              >::value
          >
>
HOSTDEVICE typename ArrT::value_type
calculate_centre_of_geometry(ArrT& pos, SizeType begin, SizeType end) {
    using value_type = typename ArrT::value_type;
    value_type cog{0., 0., 0.};
    for (SizeType i = begin; i < end; ++i)
        cog += pos(i);
    return cog / (end - begin);
}

template <typename IType>
HOSTDEVICE auto calculate_centre_of_geometry(IType pos_begin, IType pos_end) {
    using value_type = std::remove_reference_t<decltype(*pos_begin)>;
    value_type cog{0., 0., 0.};
    auto n = 0;
    for (auto it = pos_begin; it != pos_end; ++it, ++n)
        cog += *it;
    return cog / n;
}