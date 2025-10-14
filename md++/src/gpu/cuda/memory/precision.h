
#pragma once

#include "gpu/cuda/cuhostdevice.h"

// Choose one by defining FP_PRECISION as 1, 2, or 3
// 1 = float only
// 2 = double only
// 3 = mixed precision
#ifndef FP_PRECISION
  #define FP_PRECISION 1  // default to float only
#endif

#if FP_PRECISION == 1
    using FPL_TYPE = float;
    using FPH_TYPE = float;
#elif FP_PRECISION == 2
    using FPL_TYPE = float;
    using FPH_TYPE = double;
#elif FP_PRECISION == 3
    using FPL_TYPE = double;
    using FPH_TYPE = double;
#else
    #error "FP_PRECISION must be 1, 2, or 3"
#endif

template<typename FP, int N>
struct FPVec;

template<> struct FPVec<float,2>  { using type = float2; };
template<> struct FPVec<float,3>  { using type = float3; };
template<> struct FPVec<float,4>  { using type = float4; };

template<> struct FPVec<double,2> { using type = double2; };
template<> struct FPVec<double,3> { using type = double3; };
template<> struct FPVec<double,4> { using type = double4; };

// convenience aliases
template<typename FP, int N>
using FPVec_t = typename FPVec<FP,N>::type;
template<typename FP>
using FP2 = typename FPVec<FP,2>::type;
using FPL2_TYPE = FP2<FPL_TYPE>;
using FPH2_TYPE = FP2<FPH_TYPE>;

template<typename FP>
using FP3 = typename FPVec<FP,3>::type;
using FPL3_TYPE = FP3<FPL_TYPE>;
using FPH3_TYPE = FP3<FPH_TYPE>;

template<typename FP>
using FP4 = typename FPVec<FP,4>::type;
using FPL4_TYPE = FP4<FPL_TYPE>;
using FPH4_TYPE = FP4<FPH_TYPE>;


template <typename FP>
struct FPMat;

template<> struct FPMat<float>  { using type = float9; };
template<> struct FPMat<double> { using type = double9; };

template<typename FP>
using FPMat_t = typename FPMat<FP>::type;
using FPL9_TYPE = FPMat_t<FPL_TYPE>;
using FPH9_TYPE = FPMat_t<FPH_TYPE>;


/*
use:
    FPL_TYPE x;             // low precision type
    FPH_TYPE y;       // high precision type

compilation:
    # Float only
    cmake -DFP_PRECISION=1 ...

    # Mixed precision
    cmake -DFP_PRECISION=2 ...

    # Double only
    cmake -DFP_PRECISION=3 ...

*/

/**
 * @brief 
 * 
 * @tparam FP 
 * @param x 
 * @param y 
 * @param z 
 * @return HOSTDEVICE 
 */
template<typename FP>
HOSTDEVICE FPVec_t<FP,3> make_FP3(FP x, FP y, FP z) {
    if constexpr (std::is_same<FP,float>::value) {
        return make_float3(x, y, z);
    } else if constexpr (std::is_same<FP,double>::value) {
        return make_double3(x, y, z);
    } else {
        static_assert(!sizeof(FP*), "FP must be float or double");
    }
}

// Forwarding template: cast any convertible types to FP
template<typename FP, typename T,
         typename = std::enable_if_t<
             std::is_convertible_v<T, FP> &&
             !std::is_same_v<T, FP>
         >>
HOSTDEVICE FPVec_t<FP,3> make_FP3(T x, T y, T z) {
    return make_FP3(FP(x), FP(y), FP(z)); // calls terminal overload
}

// Single-value overload: sets all components to the same value
template<typename FP, typename T>
HOSTDEVICE FPVec_t<FP,3> make_FP3(T val) {
    return make_FP3(val, val, val);
}

HOSTDEVICE FPVec_t<FPL_TYPE,3> make_FPL3(FPL_TYPE x, FPL_TYPE y, FPL_TYPE z) {
    return make_FP3<FPL_TYPE>(x,y,z);
}

HOSTDEVICE FPVec_t<FPH_TYPE,3> make_FPH3(FPH_TYPE x, FPH_TYPE y, FPH_TYPE z) {
    return make_FP3<FPH_TYPE>(x,y,z);
}

HOSTDEVICE FPVec_t<FPL_TYPE,3> make_FPL3(FPL_TYPE val) {
    return make_FP3<FPL_TYPE>(val);
}

HOSTDEVICE FPVec_t<FPH_TYPE,3> make_FPH3(FPH_TYPE val) {
    return make_FP3<FPH_TYPE>(val);
}