
#pragma once


// Choose one by defining FP_PRECISION as 1, 2, or 3
// 1 = float only
// 2 = double only
// 3 = mixed precision
#ifndef FP_PRECISION
#define FP_PRECISION 1  // default to double only
#endif

#if FP_PRECISION == 1  // float only
typedef float FP_TYPE;
typedef float FP_MIXED_TYPE;  // not used here, but keep consistent

#elif FP_PRECISION == 2  // double only
typedef double FP_TYPE;
typedef double FP_MIXED_TYPE;

#elif FP_PRECISION == 3  // mixed precision
typedef double FP_TYPE;        // main working type double
typedef float FP_MIXED_TYPE;   // "mixed" type is float here

#else
#error "Unsupported FP_PRECISION"
#endif

/*
use:
    FP_TYPE x;             // main variables use FP_TYPE
    FP_MIXED_TYPE y;       // mixed precision vars use FP_MIXED_TYPE

compilation:
    # Float only
    g++ -DFP_PRECISION=1 ...

    # Double only
    g++ -DFP_PRECISION=2 ...

    # Mixed precision
    g++ -DFP_PRECISION=3 ...

*/

// Helper macros for literals and casting
#if FP_PRECISION == 1
  #define FP_LITERAL(x) x##f
#elif FP_PRECISION == 2
  #define FP_LITERAL(x) x
#else
  #define FP_LITERAL(x) x  // use double literal as default
#endif
