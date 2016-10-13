/*
 *================================================================================
 *
 *    Copyright (c) 2016 Vortex Co.,Ltd.
 *    Unpublished - All rights reserved
 *
 *================================================================================
 *    File description:
 *    Timer class to track elapsed time in portions of a program.
 *
 *================================================================================
 *    Date            Name                    Description of Change
 *    23-Sep-2016     Jiamin Xu               Creation
 *================================================================================
 */
#include "Timer.hpp"

#include "SAMRAI_MPI.h"
#include "IOStream.h"
#include "SAMRAIManager.h"
#include "TimerManager.h"
#include "Utilities.h"

namespace AREIS
{
    const int Timer::DEFAULT_NUMBER_OF_TIMERS_INCREMENT = 128;
    const int Timer::TBOX_TIMER_VERSION = 1;

    /*
     * The constructor sets the timer name and initializes timer state.
     */
    Timer::Timer(const std::string& name):
            d_name(name),
            d_is_running(false),
            d_is_active(true),
            d_accesses(0)
    {
        Clock::initialize(d_user_start_exclusive);
        Clock::initialize(d_user_stop_exclusive);
        Clock::initialize(d_system_start_exclusive);
        Clock::initialize(d_system_stop_exclusive);
        Clock::initialize(d_wallclock_start_exclusive);
        Clock::initialize(d_wallclock_stop_exclusive);

        reset();
    }
    
    Timer::~Timer()
    {
        d_concurrent_timers.clear();
    }
    
    /*
     * Start and stop routines for timers.
     *
     * For wallclock time: If we have MPI, we use MPI_Wtime to set the
     *                     start/stop point.  If we don't have MPI but do
     *                     have access to timer utilities in sys/times.h,
     *                     we use the time() utility to set the start/start
     *                     point.  If we have neither, we set the wallclock
     *                     start/stop time to zero.
     *
     * For user time:      If we have access to timer utilities in sys/times.h,
     *                     we use the times() utility to compute user and
     *                     system start/stop point (passing in the tms struct).
     *                     If we don't have these utilities, we simply set the
     *                     user and start/stop times to zero.
     *
     * Note that the stop routine increments the elapsed time information.
     * Also, the timer manager manipulates the exclusive time information
     * the timers when start and stop are called.
     */
    void Timer::start()
    {
        if (d_is_active)
        {
            if (d_is_running == true)
            {
                ARIES_ERROR("Illegal attempt to start timer '"
                            << d_name
                            << "' when it is already started.");
            }
            d_is_running = true;

            ++d_accesses;
            Clock::timestamp(d_user_start_total,
                             d_system_start_total,
                             d_wallclock_start_total);

            TimerManager::GetManager()->startTime(this);
        }
    }

    void Timer::stop()
    {
        if (d_is_active)
        {
            if (d_is_running == false)
            {
                ARIES_ERROR("Illegal attempt to stop timer '"
                            << d_name
                            << "' when it is already stopped.");
            }
            d_is_running = false;

            TimerManager::getManager()->stopTime(this);
            
            Clock::timestamp(d_user_stop_total,
                             d_system_stop_total,
                             d_wallclock_stop_total);

            d_wallclock_total += double(d_wallclock_stop_total - d_wallclock_start_total);
            d_user_total += double(d_user_stop_total - d_user_start_total);
            d_system_total += double(d_system_stop_total - d_system_start_total);
        }
    }

    void Timer::startExclusive()
    {
        if (d_is_active)
        {
            Clock::timestamp(d_user_start_exclusive,
                             d_system_start_exclusive,
                             d_wallclock_start_exclusive);
        }
    }

    void Timer::stopExclusive()
    {
        if (d_is_active)
        {
            Clock::timestamp(d_user_stop_exclusive,
                             d_system_stop_exclusive,
                             d_wallclock_stop_exclusive);

            d_wallclock_exclusive += double(d_wallclock_stop_exclusive - d_wallclock_start_exclusive);
            d_user_exclusive += double(d_user_stop_exclusive - d_user_start_exclusive);
            d_system_exclusive += double(d_system_stop_exclusive - d_system_start_exclusive);
        }
    }
    
    void Timer::barrierAndStart()
    {
        if (d_is_active)
            SAMRAI_MPI::getSAMRAIWorld().Barrier();
        
        start();
    }

    void Timer::barrierAndStop()
    {
        if (d_is_active)
            SAMRAI_MPI::getSAMRAIWorld().Barrier();
        
        stop();
    }

    void Timer::reset()
    {
        d_user_total = 0.0;
        d_system_total = 0.0;
        d_wallclock_total = 0.0;

        d_user_exclusive = 0.0;
        d_system_exclusive = 0.0;
        d_wallclock_exclusive = 0.0;

        d_max_wallclock = 0.0;

        d_concurrent_timers.clear();
    }

    bool Timer::isConcurrentTimer(const Timer& timer) const
    {
        for (std::vector<const Timer *>::const_iterator i = d_concurrent_timers.begin();
             i != d_concurrent_timers.end();
             ++i)
        {
            if (*i == &timer)
            {
                return true;
            }
        }

        return false;
    }

    /*
     * Compute the load balance efficiency based the wallclock time on each
     * processor, using the formula:
     *
     *      eff = (sum(time summed across processors)/#processors) /
     *             max(time across all processors)
     *
     * This formula corresponds to that used to compute load balance
     * efficiency based on the processor distribution of the the number of
     * cells (i.e. in BalanceUtilities::computeLoadBalanceEfficiency).
     */
    double Timer::computeLoadBalanceEfficiency()
    {
        const SAMRAI_MPI& mpi(SAMRAI_MPI::getSAMRAIWorld());
        double wall_time = d_wallclock_total;
        double sum = wall_time;
        if (mpi.getSize() > 1)
        {
            mpi.Allreduce(&wall_time, &sum, 1, MPI_DOUBLE, MPI_SUM);
        }
        computeMaxWallclock();
        int nprocs = mpi.getSize();
        double eff = 100.;
        if (d_max_wallclock > 0.)
        {
            eff = 100. * (sum / (double)nprocs) / d_max_wallclock;
        }
        return eff;
    }

    void Timer::computeMaxWallclock()
    {
        const SAMRAI_MPI& mpi(SAMRAI_MPI::getSAMRAIWorld());
        double wall_time = d_wallclock_total;
        if (mpi.getSize() > 1)
        {
            mpi.Allreduce( &wall_time, &d_max_wallclock, 1, MPI_DOUBLE, MPI_MAX);
        }
    }

    void Timer::putToRestart() const
    {
        //TBOX_ASSERT(restart_db);
        //
        //restart_db->putInteger("TBOX_TIMER_VERSION", TBOX_TIMER_VERSION);
        //
        //restart_db->putString("d_name", d_name);
        //
        //restart_db->putDouble("d_user_total", d_user_total);
        //restart_db->putDouble("d_system_total", d_system_total);
        //restart_db->putDouble("d_wallclock_total", d_wallclock_total);
        //
        //restart_db->putDouble("d_user_exclusive", d_user_exclusive);
        //restart_db->putDouble("d_system_exclusive", d_system_exclusive);
        //restart_db->putDouble("d_wallclock_exclusive", d_wallclock_exclusive);
    }

    void Timer::getFromRestart()
    {
        //TBOX_ASSERT(restart_db);
        //
        //int ver = restart_db->getInteger("TBOX_TIMER_VERSION");
        //if (ver != TBOX_TIMER_VERSION)
        //    TBOX_ERROR("Restart file version different than class version.");
        //
        //
        //d_name = restart_db->getString("d_name");
        //
        //d_user_total = restart_db->getDouble("d_user_total");
        //d_system_total = restart_db->getDouble("d_system_total");
        //d_wallclock_total = restart_db->getDouble("d_wallclock_total");
        //
        //d_user_exclusive = restart_db->getDouble("d_user_exclusive");
        //d_system_exclusive = restart_db->getDouble("d_system_exclusive");
        //d_wallclock_exclusive = restart_db->getDouble("d_wallclock_exclusive");
    }
}
