/*
 *================================================================================
 *
 *    Copyright (c) 2016 Vortex Co.,Ltd.
 *    Unpublished - All rights reserved
 *
 *================================================================================
 *    File description:
 *    Utility class for logging.
 *
 *================================================================================
 *    Date            Name                    Description of Change
 *    05-Oct-2016     Jiamin Xu               Creation
 *================================================================================
 */
#ifndef ARIES_CLOCK_HPP
#define ARIES_CLOCK_HPP

/* Aries includes */
#include "AriesMPI.hpp"
/* c++ includes */
#include <ctime>
#include <sys/times.h>
#include <unistd.h>

namespace ARIES
{
    class Clock
    {
        /**
         * Initialize system clock.  Argument must be in the "clock_t" format
         * which is a standard POSIX struct provided on most systems in the
         * <sys/times.h> include file. On Microsoft systems, it is provided in
         * <time.h>.
         */
        static void Initialize(clock_t& clock) { clock = times(&d_tmsBuffer); }

        /**
         * Initialize system clock, where clock is in double format.
         */
        static void Initialize(double& clock) { clock = 0.0; }

        /**
         * Timestamp user, system, and walltime clocks.  Wallclock argument is in
         * double format since it will access wallclock times from
         * SAMRAI_MPI::Wtime() function.
         */
        static void Timestamp(clock_t& user, clock_t& sys, double& wall)
        {
            d_nullClock_t = times(&d_tmsBuffer);
            wall = AriesMPI::Wtime();
            sys = d_tmsBuffer.tms_stime;
            user = d_tmsBuffer.tms_utime;
        }

        /**
         * Returns clock cycle for the system.
         */
        static double GetClockCycle()
        {
            double clock_cycle = 1.0;
            return clock_cycle;
        }

    private:
        static struct tms d_tmsBuffer;
        static clock_t d_nullClock_t;
    };
}

#endif
