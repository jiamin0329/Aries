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
 *    16-Oct-2016     Jiamin Xu               Creation
 *================================================================================
 */

#include "Timer.hpp"

/* Aries includes */
#include "AriesMPI.hpp"
#include "IOStream.hpp"
#include "AriesManager.hpp"
#include "TimerManager.hpp"
#include "Logger.hpp"


namespace ARIES
{
    const int Timer::DEFAULT_NUMBER_OF_TIMERS_INCREMENT = 128;
    const int Timer::ARIES_TIMER_VERSION = 1;

    Timer::Timer(const std::string& name):
            d_name(name),
            d_isRunning(false),
            d_isActive(true),
            d_accesses(0)
    {
        Clock::Initialize(d_userStartExclusive);
        Clock::Initialize(d_userStopExclusive);
        Clock::Initialize(d_systemStartExclusive);
        Clock::Initialize(d_systemStopExclusive);
        Clock::Initialize(d_wallclockStartExclusive);
        Clock::Initialize(d_wallclockStopExclusive);

        Reset();
    }
    
    Timer::~Timer()
    {
        d_concurrentTimers.clear();
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
    void Timer::Start()
    {
        if (d_isActive)
        {
            if (d_isRunning == true)
            {
                ARIES_ERROR("Illegal attempt to start timer '" << d_name << "' when it is already started.");
            }
            d_isRunning = true;

            ++d_accesses;
            Clock::Timestamp(d_userStartTotal,
                             d_systemStartTotal,
                             d_wallclockStartTotal);

            TimerManager::GetInstance()->StartTime(this);
        }
    }

    void Timer::Stop()
    {
        if (d_isActive)
        {
            if (d_isRunning == false)
            {
                ARIES_ERROR("Illegal attempt to stop timer '" << d_name << "' when it is already stopped.");
            }
            d_isRunning = false;

            TimerManager::GetInstance()->StopTime(this);
            
            Clock::Timestamp(d_userStopTotal,
                             d_systemStopTotal,
                             d_wallclockStopTotal);

            d_wallclockTotal += double(d_wallclockStopTotal - d_wallclockStartTotal);
            d_userTotal      += double(     d_userStopTotal -      d_userStartTotal);
            d_systemTotal    += double(   d_systemStopTotal -    d_systemStartTotal);
        }
    }

    void Timer::StartExclusive()
    {
        if (d_isActive)
        {
            Clock::Timestamp(d_userStartExclusive,
                             d_systemStartExclusive,
                             d_wallclockStartExclusive);
        }
    }

    void Timer::StopExclusive()
    {
        if (d_isActive)
        {
            Clock::Timestamp(d_userStopExclusive,
                             d_systemStopExclusive,
                             d_wallclockStopExclusive);

            d_wallclockExclusive += double(d_wallclockStopExclusive - d_wallclockStartExclusive);
            d_userExclusive      += double(     d_userStopExclusive -      d_userStartExclusive);
            d_systemExclusive    += double(   d_systemStopExclusive -    d_systemStartExclusive);
        }
    }
    
    void Timer::BarrierAndStart()
    {
        if (d_isActive)
            AriesMPI::GetAriesWorld().Barrier();
        
        Start();
    }

    void Timer::BarrierAndStop()
    {
        if (d_isActive)
            AriesMPI::GetAriesWorld().Barrier();
        
        Stop();
    }

    void Timer::Reset()
    {
        d_userTotal = 0.0;
        d_systemTotal = 0.0;
        d_wallclockTotal = 0.0;

        d_userExclusive = 0.0;
        d_systemExclusive = 0.0;
        d_wallclockExclusive = 0.0;

        d_maxWallclock = 0.0;

        d_concurrentTimers.clear();
    }

    bool Timer::IsConcurrentTimer(const Timer& timer) const
    {
        for (std::vector<const Timer *>::const_iterator i = d_concurrentTimers.begin();
             i != d_concurrentTimers.end();
             ++i)
        {
            if (*i == &timer)
            {
                return true;
            }
        }

        return false;
    }

    double Timer::ComputeLoadBalanceEfficiency()
    {
        const AriesMPI& mpi(AriesMPI::GetAriesWorld());
        double wall_time = d_wallclockTotal;
        double sum = wall_time;
        if (mpi.GetSize() > 1)
        {
            mpi.Allreduce(&wall_time, &sum, 1, MPI_DOUBLE, MPI_SUM);
        }
        ComputeMaxWallclock();
        int nprocs = mpi.GetSize();
        double eff = 100.;
        if (d_maxWallclock > 0.)
        {
            eff = 100. * (sum / (double)nprocs) / d_maxWallclock;
        }
        return eff;
    }

    void Timer::ComputeMaxWallclock()
    {
        const AriesMPI& mpi(AriesMPI::GetAriesWorld());
        double wall_time = d_wallclockTotal;
        if (mpi.GetSize() > 1)
        {
            mpi.Allreduce( &wall_time, &d_maxWallclock, 1, MPI_DOUBLE, MPI_MAX);
        }
    }

    void Timer::PutToRestart() const
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
        //restart_db->putDouble("d_userExclusive", d_userExclusive);
        //restart_db->putDouble("d_systemExclusive", d_systemExclusive);
        //restart_db->putDouble("d_wallclockExclusive", d_wallclockExclusive);
    }

    void Timer::GetFromRestart()
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
        //d_userExclusive = restart_db->getDouble("d_userExclusive");
        //d_systemExclusive = restart_db->getDouble("d_systemExclusive");
        //d_wallclockExclusive = restart_db->getDouble("d_wallclockExclusive");
    }
    
} // end namespace ARIES
