/*
 *================================================================================
 *
 *    Copyright (c) 2016 Vortex Co.,Ltd.
 *    Unpublished - All rights reserved
 *
 *================================================================================
 *    File description:
 *    Utility functions for using OpenMP.
 *
 *================================================================================
 *    Date            Name                    Description of Change
 *    07-Oct-2016     Jiamin Xu               Creation
 *================================================================================
 */
#ifndef ARIES_ARIESOMP_HPP
#define ARIES_ARIESOMP_HPP

/*
 *  @brief Macros defined for using OpenMP, with sensible definitions
 *  when not using it.
 *
 *  OpenMP symbols beginning with omp_ is prepended with ARIES_ to
 *  indicate it is a ARIES macro.  Macros for OpenMP functions
 *  must have argument list, even if it is empty.
 */

#ifdef _OPENMP

#include <omp.h>

#define ARIES_omp_version _OPENMP
#define ARIES_omp_lock_t omp_lock_t

#define ARIES_omp_init_lock(LOCK_PTR) omp_init_lock(LOCK_PTR)
#define ARIES_omp_destroy_lock(LOCK_PTR) omp_destroy_lock(LOCK_PTR)

#define ARIES_omp_set_lock(LOCK_PTR) omp_set_lock(LOCK_PTR)
#define ARIES_omp_unset_lock(LOCK_PTR) omp_unset_lock(LOCK_PTR)

#define ARIES_omp_get_num_threads() omp_get_num_threads()
#define ARIES_omp_get_max_threads() omp_get_max_threads()

#define ARIES_IF_SINGLE_THREAD(CODE)            \
    {                                           \
        if (omp_get_num_threads() == 1) {       \
            CODE                                \
        }                                       \
    }

#define ARIES_IF_MULTI_THREAD(CODE)             \
    {                                           \
        if (omp_get_num_threads() > 1) {        \
            CODE                                \
        }                                       \
    }

#define ARIES_IF_IN_PARALLEL_REGION(CODE)       \
    {                                           \
        if (omp_in_parallel()) {                \
            CODE                                \
        }                                       \
    }

#define ARIES_IF_NOT_IN_PARALLEL_REGION(CODE)   \
    {                                           \
        if (!omp_in_parallel()) {               \
            CODE                                \
        }                                       \
    }

#define ARIES_IF_HAVE_OPENMP(CODE) { CODE }
#define ARIES_IF_NOT_HAVE_OPENMP(CODE)

#else

#define ARIES_omp_version 0
#define ARIES_omp_lock_t int

#define ARIES_omp_init_lock(LOCK_PTR)
#define ARIES_omp_destroy_lock(LOCK_PTR)

#define ARIES_omp_set_lock(LOCK_PTR)
#define ARIES_omp_unset_lock(LOCK_PTR)

#define ARIES_omp_get_num_threads() (1)
#define ARIES_omp_get_max_threads() (1)

#define ARIES_IF_SINGLE_THREAD(CODE) { CODE }
#define ARIES_IF_MULTI_THREAD(CODE)
#define ARIES_IF_IN_PARALLEL_REGION(CODE)
#define ARIES_IF_NOT_IN_PARALLEL_REGION(CODE) { CODE }
#define ARIES_IF_HAVE_OPENMP(CODE)
#define ARIES_IF_NOT_HAVE_OPENMP(CODE) { CODE }

#endif

#endif

