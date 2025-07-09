/**
 * @file constants.h
 * device constant memory
 */

#pragma once

namespace device {
    /**
     * @struct cutoff_struct
     * holds cutoffs in device constant memory
     */
    template<typename R>
    struct cutoff_struct {
        /**
         * the long-range cutoff
         */
        R clong;
        /**
         * the squared long-range cutoff
         */
        R clong_2;
        /**
         * the short-range cutoff
         */
        R cshort;
        /**
         * the squared short-range cutoff
         */
        R cshort_2;
    };

    __constant__ __device__ cutoff_struct<float> cutoff;
    //__constant__ __device__ float cutoff_short;
    //__constant__ __device__ float cutoff_long;
    texture<float, 1, cudaReadModeElementType> tex;
}
